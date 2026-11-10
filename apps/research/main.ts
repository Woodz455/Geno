/** Mode Recherche (semaine 17) : la carte de contacts d'un utilisateur → une structure.
 *
 *  La page ne calcule rien : elle parle au service local (`geno research serve`), qui lit le
 *  `.mcool` avec cooler et reconstruit. Elle importe (fichier déposé, ou chemin sur la machine),
 *  montre ce que le fichier contient, règle les paramètres, lance, affiche la structure et ce
 *  qu'elle vaut, donne le **manifeste de session** — et sait en rejouer un : mêmes octets, ou
 *  elle dit lesquels diffèrent.
 */

import { G3D, HttpSource } from "@geno/g3d";
import { createTargets, distanceFor, ImpostorRenderer, NO_HIT, type Orbit, perspective, Picker, rotate, Ssao, viewMatrix, widthAt, zoom } from "@geno/viewer";

const FOV = (45 * Math.PI) / 180;

declare global {
  interface Window {
    __geno?: Record<string, unknown>;
    __research?: Record<string, unknown>;
  }
}

interface Resolution { resolution: number; bins: number; balanced: boolean; contacts: number }
interface Inspect { name: string; bytes: number; assembly: string | null; chroms: [string, number][]; resolutions: Resolution[]; max_beads: number }
interface Manifest {
  format: string;
  input: { name: string; bytes: number; sha256: string };
  params: { resolution: number; region: string; alpha: number; balance: string; min_coverage: number };
  region: string;
  versions: Record<string, string>;
  output: { file: string; bytes: number; sha256: string };
  stats: { bins: number; beads: number; masked: number; contacts: number; observed_pairs: number; pieces: number; r_log_distance: number; radius_nm: number; gyration_nm: number };
  run: { created: string; seconds: number };
}
interface Job { state: string; step: string; progress: number; session?: string; manifest?: Manifest; error?: string;
                reproduce?: { identical: boolean; expected_sha256: string; sha256: string; versions_differ: Record<string, [string, string]>; note?: string } }

const style = document.createElement("style");
style.textContent = `
  html, body { height: 100%; overflow: hidden; }
  #side { position: fixed; left: 0; top: 0; bottom: 0; width: min(380px, 100vw); overflow-y: auto; box-sizing: border-box;
          background: rgba(14,18,24,.92); padding: 12px 14px; font-size: 12.5px; line-height: 1.45; }
  #side h1 { font-size: 14px; margin: 0 0 4px; }
  #side h2 { font-size: 12px; margin: 14px 0 4px; text-transform: uppercase; letter-spacing: .04em; color: #8b98aa; }
  #side label { display: block; margin-top: 6px; color: #c9d3df; }
  #side input, #side select, #side button { font: inherit; color: inherit; background: #0e131b; border: 1px solid #3a4658;
          border-radius: 6px; padding: 4px 8px; box-sizing: border-box; max-width: 100%; }
  #side input[type=text], #side input[type=number], #side select { width: 100%; }
  #side button { background: #1b2330; cursor: pointer; margin-top: 8px; }
  #side button:disabled { opacity: .5; cursor: default; }
  .k { color: #8b98aa; } .err { color: #ff9c8a; } .ok { color: #8fd694; }
  #bar { height: 4px; background: #1b2330; border-radius: 2px; margin-top: 8px; }
  #bar div { height: 100%; width: 0; background: #6aa7ff; border-radius: 2px; }
  #hover { position: fixed; right: 12px; bottom: 10px; background: rgba(14,18,24,.85); padding: 6px 10px; border-radius: 6px; font-size: 12px; }
  canvas { position: fixed; left: 0; top: 0; touch-action: none; }
  table { border-collapse: collapse; } td { padding: 1px 10px 1px 0; vertical-align: top; }
`;
document.head.appendChild(style);

const canvas = document.createElement("canvas");
document.body.appendChild(canvas);
const side = document.createElement("div");
side.id = "side";
side.innerHTML = `
  <h1>Geno — mode Recherche</h1>
  <div class="k">Votre carte de contacts → une structure 3D. Le calcul tourne sur cette machine
  (service local) ; rien n'en sort.</div>
  <h2>1 · Importer</h2>
  <label>Un fichier .mcool <input id="file" type="file" accept=".mcool,.cool"></label>
  <label>ou son chemin sur cette machine <input id="path" type="text" placeholder="/chemin/vers/carte.mcool" spellcheck="false"></label>
  <button id="open" type="button">Ouvrir ce chemin</button>
  <div id="info" aria-live="polite"></div>
  <h2>2 · Paramètres</h2>
  <label>Résolution <select id="resolution" disabled></select></label>
  <label>Région — un chromosome, chrom:début-fin, ou « génome » <input id="region" type="text" spellcheck="false" disabled></label>
  <label>Équilibrage <select id="balance" disabled>
    <option value="fichier">poids du fichier (ICE déjà calculé)</option>
    <option value="ice">ICE recalculé sur la région</option>
    <option value="aucun">aucun (comptages bruts)</option></select></label>
  <label>alpha — d = f<sup>−alpha</sup> <input id="alpha" type="number" min="0.05" max="2" step="0.01" value="0.35" disabled></label>
  <div class="k">alpha n'est pas une constante : il dépend de la profondeur de la carte (mesuré, semaine 5).</div>
  <button id="run" type="button" disabled>Reconstruire</button>
  <div id="bar"><div></div></div>
  <div id="progress" aria-live="polite"></div>
  <h2>3 · Structure</h2>
  <div id="stats" aria-live="polite" class="k">aucune structure encore</div>
  <a id="manifest-dl" hidden>Télécharger le manifeste de session</a>
  <h2>4 · Rejouer un manifeste</h2>
  <div class="k">Sur le fichier importé ci-dessus : mêmes paramètres, puis comparaison des octets.</div>
  <label>Manifeste (session.json) <input id="replay" type="file" accept=".json,application/json" disabled></label>
  <div id="verdict" aria-live="polite"></div>`;
document.body.appendChild(side);
const hover = document.createElement("div");
hover.id = "hover";
hover.hidden = true;
document.body.appendChild(hover);

const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;
const fmt = (n: number) => n.toLocaleString("fr");
const esc = (s: string) => s.replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]!);

async function api<T>(path: string, body?: BodyInit, json = true): Promise<T> {
  const r = await fetch(path, body === undefined ? {} : { method: "POST", body, headers: json ? { "Content-Type": "application/json" } : {} });
  const data = await r.json();
  if (!r.ok) throw new Error(data.error ?? `erreur ${r.status}`);
  return data as T;
}

// — La vue 3D —
const gl = canvas.getContext("webgl2", { antialias: false })!;
let W = 1, H = 1;
let targets = createTargets(gl, 1, 1);
let ssao = new Ssao(gl, targets);
let picker = new Picker(gl, targets);
const renderer = new ImpostorRenderer(gl);
let orbit: Orbit = { target: [0, 0, 0], yaw: 0.4, pitch: 0.3, distance: 10 };
let reach = 1, radiusUm = 0.1, count = 0, dirty = true;
let beads: { chrom: string; start: number; end: number }[] = [];

function fit(): void {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  W = Math.max(1, Math.round(window.innerWidth * dpr));
  H = Math.max(1, Math.round(window.innerHeight * dpr));
  canvas.width = W;
  canvas.height = H;
  canvas.style.width = `${window.innerWidth}px`;
  canvas.style.height = `${window.innerHeight}px`;
  ssao.dispose();
  targets = createTargets(gl, W, H);
  ssao = new Ssao(gl, targets);
  picker = new Picker(gl, targets);
  dirty = true;
}

function render(): void {
  const d = orbit.distance;
  const proj = perspective(FOV, W / H, Math.max(0.02 * d, d - reach - 1), d + reach + 1);
  gl.bindFramebuffer(gl.FRAMEBUFFER, targets.fbo);
  gl.viewport(0, 0, W, H);
  gl.clearBufferfv(gl.COLOR, 0, [0.055, 0.07, 0.095, 1]);
  gl.clearBufferuiv(gl.COLOR, 1, [NO_HIT, 0, 0, 0]);
  gl.clearBufferfv(gl.COLOR, 2, [1, 0, 0, 1]);
  gl.clearBufferfi(gl.DEPTH_STENCIL, 0, 1.0, 0);
  if (count) renderer.draw(proj, { view: viewMatrix(orbit) });
  ssao.occlusion(proj, { radius: 3 * radiusUm });
  ssao.composite(proj, { radius: 3 * radiusUm, strength: 0.9 });
}

/** Une teinte par chromosome, et du clair au foncé le long de chacun : l'ordre de la chaîne se lit. */
function color(c: number, k: number, t: number): [number, number, number] {
  const h = (c / Math.max(1, k)) * 6.283 + 0.6;
  const l = 1.0 - 0.7 * t;
  return [l * (0.6 + 0.4 * Math.cos(h)), l * (0.6 + 0.4 * Math.cos(h - 2.094)), l * (0.6 + 0.4 * Math.cos(h + 2.094))];
}

async function show(session: string): Promise<void> {
  const g3d = await G3D.open(new HttpSource(`sessions/${session}/structure.g3d`));
  const frame = await g3d.firstFrame();
  const lv = frame.level;
  const n = frame.count;
  const centers = new Float32Array(n * 3), radii = new Float32Array(n), ids = new Uint32Array(n), colors = new Float32Array(n * 3);
  const per = new Array<number>(lv.copies.length).fill(0), seen = new Array<number>(lv.copies.length).fill(0);
  for (let i = 0; i < n; i++) per[frame.copy[i]]++;
  reach = 0;
  beads = [];
  for (let i = 0; i < n; i++) {
    const c = frame.copy[i];
    for (let k = 0; k < 3; k++) centers[3 * i + k] = frame.positions[3 * i + k] / 1000;
    reach = Math.max(reach, Math.hypot(centers[3 * i], centers[3 * i + 1], centers[3 * i + 2]));
    radii[i] = lv.radius_nm[c] / 1000;
    ids[i] = i + 1;
    colors.set(color(c, lv.copies.length, seen[c]++ / Math.max(1, per[c] - 1)), 3 * i);
    const b = await g3d.identify(lv.name, i + 1);
    beads.push({ chrom: lv.copies[c], start: b.start, end: b.end });
  }
  radiusUm = lv.radius_nm[0] / 1000;
  renderer.upload({ centers, radii, ids, colors, count: n });
  count = n;
  // La structure entière dans la partie de l'écran que le panneau ne couvre pas.
  orbit = { target: [0, 0, 0], yaw: 0.4, pitch: 0.3, distance: distanceFor(2.6 * (reach + radiusUm), FOV) };
  dirty = true;
}

// — État et gestes —
let file: { token: string; inspect: Inspect; sha256: string } | null = null;
let last: { session: string; manifest: Manifest } | null = null;
let busy = false;
const info = $<HTMLDivElement>("info"), progress = $<HTMLDivElement>("progress"), stats = $<HTMLDivElement>("stats"), verdict = $<HTMLDivElement>("verdict");
const selRes = $<HTMLSelectElement>("resolution"), inRegion = $<HTMLInputElement>("region"), selBal = $<HTMLSelectElement>("balance"), inAlpha = $<HTMLInputElement>("alpha");
const btnRun = $<HTMLButtonElement>("run"), inReplay = $<HTMLInputElement>("replay"), bar = $<HTMLDivElement>("bar").firstElementChild as HTMLDivElement;
const dl = $<HTMLAnchorElement>("manifest-dl");

function imported(r: { file: string; sha256: string; inspect: Inspect }): void {
  file = { token: r.file, inspect: r.inspect, sha256: r.sha256 };
  const I = r.inspect;
  info.innerHTML = `<b>${esc(I.name)}</b> · ${(I.bytes / 1e6).toFixed(1)} Mo · assemblage ${esc(I.assembly ?? "non déclaré")}<br>` +
    `<span class="k">${I.chroms.length} chromosome(s) : ${I.chroms.slice(0, 6).map(([c, s]) => `${esc(c)} (${(s / 1e6).toFixed(1)} Mb)`).join(", ")}${I.chroms.length > 6 ? "…" : ""}<br>` +
    `sha256 ${r.sha256.slice(0, 16)}…</span>`;
  selRes.innerHTML = I.resolutions.map((x) => `<option value="${x.resolution}">${fmt(x.resolution)} pb · ${fmt(x.bins)} bins${x.balanced ? "" : " · non équilibré"}</option>`).join("");
  // Par défaut : la résolution la plus fine où le premier chromosome tient sous la limite de billes.
  const first = I.chroms[0];
  const okRes = I.resolutions.find((x) => first[1] / x.resolution <= I.max_beads) ?? I.resolutions[I.resolutions.length - 1];
  selRes.value = String(okRes.resolution);
  inRegion.value = first[0];
  if (!okRes.balanced) selBal.value = "ice";
  for (const el of [selRes, inRegion, selBal, inAlpha, btnRun, inReplay]) el.disabled = false;
}

async function follow(job: string): Promise<Job> {
  for (;;) {
    const st = await api<Job>(`api/job/${job}`);
    bar.style.width = `${Math.round(st.progress * 100)}%`;
    progress.textContent = st.state === "en cours" || st.state === "en attente" ? `${st.state} — ${st.step}` : "";
    if (st.state !== "en cours" && st.state !== "en attente") return st;
    await new Promise((r) => setTimeout(r, 120));
  }
}

function describe(m: Manifest): void {
  const s = m.stats;
  stats.className = "";
  stats.innerHTML = `<table>
    <tr><td class="k">région</td><td>${esc(m.region)} à ${fmt(m.params.resolution)} pb</td></tr>
    <tr><td class="k">billes</td><td>${fmt(s.beads)} <span class="k">(${fmt(s.masked)} bins écartés : couverture trop faible ou sans poids)</span></td></tr>
    <tr><td class="k">contacts lus</td><td>${fmt(Math.round(s.contacts))} · ${(s.observed_pairs * 100).toFixed(0)} % des paires observées</td></tr>
    <tr><td class="k">accord</td><td>r = ${s.r_log_distance.toFixed(3)} <span class="k">entre distances demandées et distances de la structure (log) — reproduire sa carte n'est pas être juste</span></td></tr>
    <tr><td class="k">taille</td><td>rayon de giration ${fmt(Math.round(s.gyration_nm))} nm <span class="k">(échelle de convention : deux billes voisines se touchent)</span></td></tr>
    <tr><td class="k">calcul</td><td>${m.run.seconds.toFixed(1)} s · alpha ${m.params.alpha} · équilibrage « ${esc(m.params.balance)} »</td></tr>
    <tr><td class="k">.g3d</td><td>${fmt(m.output.bytes)} octets · sha256 ${m.output.sha256.slice(0, 16)}…</td></tr></table>
    ${s.pieces > 1 ? `<div class="err">${s.pieces} morceaux sans aucun contact entre eux : leur position relative est arbitraire, seule la forme de chacun se lit.</div>` : ""}
    <div class="k">Reconstruite (ShRec3D + MDS) : une structure consensus, à une isométrie près — orientation de convention, chiralité arbitraire.</div>`;
  dl.hidden = false;
  dl.href = URL.createObjectURL(new Blob([JSON.stringify(m, null, 1)], { type: "application/json" }));
  dl.download = "session.json";
}

async function guard(where: HTMLElement, f: () => Promise<void>): Promise<void> {
  if (busy) return;
  busy = true;
  btnRun.disabled = true;
  try {
    await f();
  } catch (e) {
    where.innerHTML = `<span class="err">${esc(e instanceof Error ? e.message : String(e))}</span>`;
  } finally {
    busy = false;
    btnRun.disabled = !file;
  }
}

$<HTMLInputElement>("file").addEventListener("change", (ev) => {
  const f = (ev.target as HTMLInputElement).files?.[0];
  if (!f) return;
  void guard(info, async () => {
    info.textContent = `envoi de ${f.name} au service local…`;
    imported(await api(`api/upload?name=${encodeURIComponent(f.name.replace(/[^A-Za-z0-9._-]/g, "_"))}`, f, false));
  });
});
$<HTMLButtonElement>("open").addEventListener("click", () => {
  void guard(info, async () => {
    info.textContent = "lecture…";
    imported(await api("api/open", JSON.stringify({ path: $<HTMLInputElement>("path").value })));
  });
});
btnRun.addEventListener("click", () => {
  void guard(progress, async () => {
    verdict.textContent = "";
    const params = { resolution: Number(selRes.value), region: inRegion.value, alpha: Number(inAlpha.value), balance: selBal.value };
    const { job } = await api<{ job: string }>("api/run", JSON.stringify({ file: file!.token, params }));
    const st = await follow(job);
    if (st.state !== "fait") throw new Error(st.error ?? st.state);
    last = { session: st.session!, manifest: st.manifest! };
    describe(st.manifest!);
    await show(st.session!);
  });
});
inReplay.addEventListener("change", (ev) => {
  const f = (ev.target as HTMLInputElement).files?.[0];
  if (!f) return;
  void guard(verdict, async () => {
    verdict.textContent = "rejeu…";
    let manifest: unknown;
    try {
      manifest = JSON.parse(await f.text());
    } catch {
      throw new Error("ce fichier n'est pas du JSON");
    }
    const { job } = await api<{ job: string }>("api/reproduce", JSON.stringify({ file: file!.token, manifest }));
    const st = await follow(job);
    if (st.state !== "fait") throw new Error(st.error ?? st.state);
    const r = st.reproduce!;
    last = { session: st.session!, manifest: st.manifest! };
    describe(st.manifest!);
    await show(st.session!);
    const diff = Object.entries(r.versions_differ).map(([k, [a, b]]) => `${k} ${a} → ${b}`).join(", ");
    verdict.innerHTML = r.identical
      ? `<span class="ok">Reproduit à l'identique : le .g3d a les mêmes octets (sha256 ${r.sha256.slice(0, 16)}…).</span>`
      : `<span class="err">Différent : sha256 ${r.sha256.slice(0, 16)}… au lieu de ${r.expected_sha256.slice(0, 16)}…</span><br><span class="k">${esc(r.note ?? "")}${diff ? ` — ${esc(diff)}` : ""}</span>`;
    (ev.target as HTMLInputElement).value = "";
  });
});

let drag: { x: number; y: number } | null = null;
canvas.addEventListener("pointerdown", (e) => {
  drag = { x: e.clientX, y: e.clientY };
  canvas.setPointerCapture(e.pointerId);
});
canvas.addEventListener("pointermove", (e) => {
  if (drag) {
    orbit = rotate(orbit, -(e.clientX - drag.x) * 0.006, (e.clientY - drag.y) * 0.006);
    drag = { x: e.clientX, y: e.clientY };
    dirty = true;
    return;
  }
  if (!count || dirty) return;
  const dpr = W / canvas.clientWidth;
  const id = picker.pick(e.offsetX * dpr, e.offsetY * dpr);
  const b = id !== NO_HIT ? beads[id - 1] : null;
  hover.hidden = !b;
  if (b) hover.textContent = `${b.chrom}:${fmt(b.start + 1)}–${fmt(b.end)} · bille ${id}`;
});
canvas.addEventListener("pointerup", () => (drag = null));
canvas.addEventListener("wheel", (e) => {
  e.preventDefault();
  orbit = zoom(orbit, Math.exp(e.deltaY * 0.0015), 0.05, 400);
  dirty = true;
}, { passive: false });
window.addEventListener("resize", fit);

fit();
const loop = () => {
  if (dirty) {
    render();
    dirty = false;
  }
  requestAnimationFrame(loop);
};
requestAnimationFrame(loop);

// Pour `tools/research.mjs` : l'état, et ce que l'image montre vraiment.
window.__research = {
  busy: () => busy || dirty,
  file: () => file,
  last: () => last,
  beads: () => beads,
  verdict: () => verdict.textContent,
  message: () => [info.textContent, progress.textContent, verdict.textContent],
  /** Billes distinctes et pixels couverts dans le tampon d'identifiants, sous deux orientations. */
  seen() {
    const ids = new Set<number>();
    let covered = 0;
    for (const yaw of [orbit.yaw, orbit.yaw + Math.PI]) {
      const keep = orbit;
      orbit = { ...orbit, yaw };
      render();
      const buf = new Uint32Array(W * H * 4);
      gl.bindFramebuffer(gl.READ_FRAMEBUFFER, targets.fbo);
      gl.readBuffer(gl.COLOR_ATTACHMENT1);
      gl.readPixels(0, 0, W, H, gl.RGBA_INTEGER, gl.UNSIGNED_INT, buf);
      gl.readBuffer(gl.COLOR_ATTACHMENT0);
      for (let i = 0; i < buf.length; i += 4) if (buf[i] !== NO_HIT) { ids.add(buf[i]); if (yaw === keep.yaw) covered++; }
      orbit = keep;
    }
    render();
    return { distinct: ids.size, covered, max: Math.max(0, ...ids), W, H, widthUm: widthAt(orbit.distance, FOV) };
  },
};
window.__geno = { ok: true, research: true };
