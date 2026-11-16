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
  params: { resolution: number; region: string; alpha: number; balance: string; min_coverage: number; ensemble?: number; seed?: number };
  region: string;
  versions: Record<string, string>;
  output: { file: string; bytes: number; sha256: string };
  files?: Record<string, { file: string; bytes: number; sha256: string }>;
  ensemble?: { replicates: number; beads: number };
  stats: { bins: number; beads: number; masked: number; contacts: number; observed_pairs: number; pieces: number; r_log_distance: number; radius_nm: number; gyration_nm: number };
  run: { created: string; seconds: number };
}
interface Job { state: string; step: string; progress: number; session?: string; manifest?: Manifest; error?: string;
                reproduce?: { identical: boolean; files_compared: string[]; files_differ: string[]; expected_sha256: string; sha256: string; versions_differ: Record<string, [string, string]>; note?: string } }

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
  #side a { color: #9cc3ff; cursor: pointer; }
  #charts svg { display: block; margin-top: 6px; background: #0e131b; border-radius: 6px; }
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
  <label>Un fichier .mcool ou .hic <input id="file" type="file" accept=".mcool,.cool,.hic"></label>
  <label>ou son chemin sur cette machine <input id="path" type="text" placeholder="/chemin/vers/carte.mcool ou .hic" spellcheck="false"></label>
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
  <label>Ensemble — réplicats par rééchantillonnage de la carte (0 : aucun) <input id="ensemble" type="number" min="0" max="64" step="1" value="0" disabled></label>
  <label>Graine <input id="seed" type="number" min="0" step="1" value="0" disabled></label>
  <button id="run" type="button" disabled>Reconstruire</button>
  <div id="bar"><div></div></div>
  <div id="progress" aria-live="polite"></div>
  <h2>3 · Structure</h2>
  <div id="stats" aria-live="polite" class="k">aucune structure encore</div>
  <a id="manifest-dl" hidden>Télécharger le manifeste de session</a>
  <div id="explore" hidden>
    <h2>Ensemble</h2>
    <label><span id="replica-label">consensus</span> <input id="replica" type="range" min="0" max="0" step="1" value="0"></label>
    <label>Couleur <select id="colour"><option value="chain">le long de la chaîne</option><option value="var">variabilité (ensemble)</option></select></label>
    <div class="k">La variabilité est ce que le bruit de comptage de la carte laisse indéterminé — pas la différence d'une cellule à l'autre.</div>
  </div>
  <div id="charts"></div>
  <div id="exports" hidden>
    <h2>Exporter</h2>
    <a id="ex-pdb">PDB</a> · <a id="ex-3dg">3DG</a> · <a id="ex-gltf">glTF</a> ·
    <button id="ex-png" type="button">PNG</button> <button id="ex-svg" type="button">SVG</button>
    <div class="k">PDB, 3DG, glTF : la structure du fichier .g3d. PNG, SVG : la vue courante.</div>
  </div>
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
let beads: { chrom: string; start: number; end: number; variability: number }[] = [];
/** Semaine 18 : ce que la vue dessine — le consensus ou un réplicat —, et de quoi la recolorer. */
let cur: { session: string; consensus: Float32Array; shown: Float32Array; radii: Float32Array; ids: Uint32Array; chain: Float32Array;
           copy: Uint8Array | Uint16Array | Uint32Array | number[]; copies: string[]; replica: number } | null = null;
const replicas = new Map<string, Float32Array>();
let lastProj: Float32Array = new Float32Array(16);

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
  lastProj = proj;
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
    beads.push({ chrom: lv.copies[c], start: b.start, end: b.end, variability: b.variabilityNm });
  }
  radiusUm = lv.radius_nm[0] / 1000;
  cur = { session, consensus: centers, shown: centers, radii, ids, chain: colors, copy: Array.from(frame.copy), copies: lv.copies, replica: 0 };
  count = n;
  upload();
  // La structure entière dans la partie de l'écran que le panneau ne couvre pas.
  orbit = { target: [0, 0, 0], yaw: 0.4, pitch: 0.3, distance: distanceFor(2.6 * (reach + radiusUm), FOV) };
  dirty = true;
}

/** Du clair au foncé : la variabilité de la bille, rapportée au 95e centile de la structure. */
function varColors(): Float32Array {
  const v = beads.map((b) => b.variability);
  const top = v.slice().sort((a, b) => a - b)[Math.floor(0.95 * (v.length - 1))] || 1;
  const out = new Float32Array(v.length * 3);
  v.forEach((x, i) => {
    const t = Math.min(1, x / top);
    out.set([0.93 - 0.45 * t, 0.78 - 0.68 * t, 0.95 - 0.35 * t], 3 * i);
  });
  return out;
}

function upload(): void {
  if (!cur) return;
  const colors = $<HTMLSelectElement>("colour").value === "var" ? varColors() : cur.chain;
  renderer.upload({ centers: cur.shown, radii: cur.radii, ids: cur.ids, colors, count: cur.ids.length });
  dirty = true;
}

/** Un réplicat de l'ensemble, lu **à la demande** par une requête Range dans ensemble.f32. */
async function replica(k: number): Promise<void> {
  if (!cur) return;
  const n = cur.ids.length;
  if (k === 0) cur.shown = cur.consensus;
  else {
    const key = `${cur.session}/${k}`;
    let r = replicas.get(key);
    if (!r) {
      const at = (k - 1) * n * 12;
      const res = await fetch(`sessions/${cur.session}/ensemble.f32`, { headers: { Range: `bytes=${at}-${at + n * 12 - 1}` } });
      if (res.status !== 206) throw new Error(`réplicat ${k} : réponse ${res.status}`);
      const nm = new Float32Array(await res.arrayBuffer());
      r = Float32Array.from(nm, (x) => x / 1000);
      replicas.set(key, r);
      ensembleBytes += nm.byteLength;
    }
    cur.shown = r;
  }
  cur.replica = k;
  $<HTMLSpanElement>("replica-label").textContent = k === 0 ? "consensus" : `réplicat ${k} / ${$<HTMLInputElement>("replica").max}`;
  upload();
}
let ensembleBytes = 0;

/** Les deux courbes de statistics.json, en log-log : P(s) de la carte, R(s) de la structure. */
function chart(title: string, pts: number[][], unit: string, slope: number | null, name: string): string {
  const ok = pts.filter((p) => p[0] > 0 && p[1] > 0);
  if (ok.length < 2) return "";
  const lx = ok.map((p) => Math.log10(p[0])), ly = ok.map((p) => Math.log10(p[1]));
  const [x0, x1, y0, y1] = [Math.min(...lx), Math.max(...lx), Math.min(...ly), Math.max(...ly)];
  const X = (v: number) => 34 + ((v - x0) / (x1 - x0 || 1)) * 290, Y = (v: number) => 104 - ((v - y0) / (y1 - y0 || 1)) * 84;
  const d = ok.map((_, i) => `${i ? "L" : "M"}${X(lx[i]).toFixed(1)} ${Y(ly[i]).toFixed(1)}`).join(" ");
  return `<svg width="340" height="128" role="img" aria-label="${esc(title)}"><text x="8" y="14" fill="#c9d3df" font-size="11">${esc(title)}${slope === null ? "" : ` — ${name} ${slope.toFixed(2)}`}</text>
    <path d="${d}" fill="none" stroke="#6aa7ff" stroke-width="1.6"/>
    <text x="34" y="122" fill="#8b98aa" font-size="10">${(10 ** x0).toFixed(0)} bin</text><text x="324" y="122" fill="#8b98aa" font-size="10" text-anchor="end">${(10 ** x1).toFixed(0)} bins</text>
    <text x="4" y="106" fill="#8b98aa" font-size="10">${(10 ** y0).toPrecision(2)}</text><text x="4" y="26" fill="#8b98aa" font-size="10">${(10 ** y1).toPrecision(2)} ${unit}</text></svg>`;
}

async function extras(session: string, m: Manifest): Promise<void> {
  const charts = $<HTMLDivElement>("charts");
  charts.innerHTML = "";
  const ex = $<HTMLDivElement>("exports"), explore = $<HTMLDivElement>("explore");
  ex.hidden = !m.files;
  explore.hidden = !m.ensemble;
  if (!m.files) return;                       // un manifeste de la semaine 17 : la structure seule
  for (const k of ["pdb", "3dg", "gltf"]) {
    const a = $<HTMLAnchorElement>(`ex-${k}`);
    a.href = `sessions/${session}/${m.files[k].file}`;
    a.download = m.files[k].file;
  }
  const st = await (await fetch(`sessions/${session}/statistics.json`)).json();
  stats18 = st;
  charts.innerHTML = `<h2>Statistiques</h2>` +
    chart("contacts lus selon la séparation, P(s)", st.contacts_vs_separation, "contacts", st.ps_slope, "pente") +
    chart("distance dans la structure, R(s)", st.distance_vs_separation, "nm", st.rs_nu, "exposant") +
    (st.variability_nm ? `<div>variabilité : médiane ${fmt(Math.round(st.variability_nm.median))} nm, 90e centile ${fmt(Math.round(st.variability_nm.p90))} nm <span class="k">(rayon d'une bille : ${fmt(Math.round(m.stats.radius_nm))} nm)</span></div>` : "") +
    `<div class="k">Les deux pentes sont lues entre ${st.fit_range_bins[0]} et ${Math.round(st.fit_range_bins[1])} bins. La seconde décrit la structure reconstruite, pas le génome.</div>`;
  const slider = $<HTMLInputElement>("replica");
  slider.max = String(m.ensemble?.replicates ?? 0);
  slider.value = "0";
  $<HTMLSelectElement>("colour").value = "chain";
  $<HTMLSpanElement>("replica-label").textContent = "consensus";
}
let stats18: Record<string, unknown> | null = null;

/** La vue courante en SVG : la silhouette exacte de chaque bille, peinte du fond vers l'avant, et un
 *  masque là où une bille sécante peinte avant elle est devant — la bille visible en chaque point
 *  est celle de l'écran. La légende dit d'où vient la structure. Rend aussi les centres, pour la
 *  mesure. */
/** Une bille projetée : centre projeté (x, y) ; silhouette exacte de la sphère — une ellipse de centre
 *  (ex, ey), demi-axes rx (radial) et ry, tournée de rot degrés. Les masques sont dans le SVG. */
type SvgBead = { id: number; x: number; y: number; ex: number; ey: number; rx: number; ry: number; rot: number };

function svgExport(): { svg: string; circles: SvgBead[] } {
  const view = viewMatrix(orbit);
  const w = canvas.clientWidth, h = canvas.clientHeight;
  // Pixel = s·(u, −v) + décalage, où (u, v) = (x, y)/(−z) en espace vue ; s isotrope pour une perspective.
  const s = lastProj[5] * h / 2;
  const toPx = (u: number, v: number) => [((lastProj[0] * u - lastProj[8]) * 0.5 + 0.5) * w, (1 - ((lastProj[5] * v - lastProj[9]) * 0.5 + 0.5)) * h];
  const proj = (x: number[]) => toPx(x[0] / -x[2], x[1] / -x[2]);
  type Item = SvgBead & { power: number; c: number[]; r: number; hide: number[][][] };
  const items: Item[] = [];
  if (cur) for (let i = 0; i < cur.ids.length; i++) {
    const p = [cur.shown[3 * i], cur.shown[3 * i + 1], cur.shown[3 * i + 2]];
    const e = [0, 1, 2].map((k) => view[k] * p[0] + view[4 + k] * p[1] + view[8 + k] * p[2] + view[12 + k]);
    const r = cur.radii[i], D = -e[2], rho = Math.hypot(e[0], e[1]);
    if (D <= r) continue;                                     // sphère coupée par le plan de l'œil
    // Silhouette exacte : le cône tangent depuis l'œil coupe le plan image selon une ellipse.
    // Avec k = d² − r² (puissance de l'œil) : demi-axe transverse r/√(D²−r²), radial r√k/(D²−r²),
    // centre décalé vers l'extérieur à D·ρ/(D²−r²). Le rendu (impostors.ts) lance le même rayon.
    const k = rho * rho + D * D - r * r, q = D * D - r * r;
    const [x, y] = toPx(e[0] / D, e[1] / D);
    const off = rho > 0 ? D * rho / q : 0;
    const [ex, ey] = rho > 0 ? toPx(off * e[0] / rho, off * e[1] / rho) : [x, y];
    items.push({ id: i + 1, x, y, ex, ey, rx: s * r * Math.sqrt(k) / q, ry: s * r / Math.sqrt(q), rot: Math.atan2(-e[1], e[0]) * 180 / Math.PI, hide: [], power: k, c: e, r });
  }
  // Ordre du peintre par puissance de l'œil (d² − r²), pas par distance du centre : pour des sphères
  // disjointes, c'est un ordre de visibilité exact (théorème d'acyclicité d'Edelsbrunner, 1989).
  const back = items.slice().sort((a, b) => b.power - a.power);
  // Sphères sécantes (ici jusqu'à une centaine de voisines par bille) : aucun ordre de disques n'est
  // exact. Pour X peinte après V, V est devant exactement là où le point d'entrée du rayon dans X est à
  // l'intérieur de V (du côté de l'œil du plan radical, V ⊂ X ; de l'autre, X ⊂ V). Ce morceau de la
  // face avant de X est borné par l'arc avant du cercle d'intersection et par l'arc de son contour
  // apparent intérieur à V : projeté, il devient un masque sur X.
  // Masquer toutes les paires sécantes coûterait des centaines de Mo ; seules comptent celles qui se
  // voient. On lance donc, pour chaque pixel, le même rayon que le rendu (impostors.ts) : V est la
  // bille d'entrée la plus proche, et toute bille X touchée par ce rayon et peinte après V doit être
  // masquée par V (elle la coupe forcément : disjointe, l'ordre par puissance la mettrait derrière).
  const sub = (u: number[], v: number[]) => [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
  const dot = (u: number[], v: number[]) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  const basis = (nv: number[]) => {
    const t = Math.abs(nv[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
    const u = [nv[1] * t[2] - nv[2] * t[1], nv[2] * t[0] - nv[0] * t[2], nv[0] * t[1] - nv[1] * t[0]];
    const lu = Math.hypot(u[0], u[1], u[2]);
    const uu = u.map((x) => x / lu);
    return [uu, [nv[1] * uu[2] - nv[2] * uu[1], nv[2] * uu[0] - nv[0] * uu[2], nv[0] * uu[1] - nv[1] * uu[0]]];
  };
  const circle = (m: number[], nv: number[], rad: number, count: number) => {
    const [u, v] = basis(nv);
    return Array.from({ length: count }, (_, j) => {
      const t = (2 * Math.PI * j) / count, ct = Math.cos(t) * rad, st = Math.sin(t) * rad;
      return [m[0] + ct * u[0] + st * v[0], m[1] + ct * u[1] + st * v[1], m[2] + ct * u[2] + st * v[2]];
    });
  };
  /** La suite circulaire de points qui vérifient `keep`, dans l'ordre (un seul arc attendu). */
  const run = (pts: number[][], keep: boolean[]) => {
    const len = pts.length;
    if (keep.every(Boolean)) return { pts, closed: true };
    if (!keep.some(Boolean)) return { pts: [] as number[][], closed: false };
    const start = keep.findIndex((k, j) => k && !keep[(j + len - 1) % len]);
    const out: number[][] = [];
    for (let j = 0; j < len && keep[(start + j) % len]; j++) out.push(pts[(start + j) % len]);
    return { pts: out, closed: false };
  };
  const n = back.length;
  const CX = new Float64Array(n), CY = new Float64Array(n), CZ = new Float64Array(n), CC = new Float64Array(n);
  back.forEach((it, j) => { CX[j] = it.c[0]; CY[j] = it.c[1]; CZ[j] = it.c[2]; CC[j] = it.power; });
  const CELL = 8, gw = Math.ceil(w / CELL), gh = Math.ceil(h / CELL);
  const cells: number[][] = Array.from({ length: gw * gh }, () => []);
  back.forEach((it, j) => {
    const t = (it.rot * Math.PI) / 180, c = Math.cos(t), sn = Math.sin(t);
    const hx = Math.hypot(it.rx * c, it.ry * sn) + 1, hy = Math.hypot(it.rx * sn, it.ry * c) + 1;
    const x0 = Math.max(0, Math.floor((it.ex - hx) / CELL)), x1 = Math.min(gw - 1, Math.floor((it.ex + hx) / CELL));
    const y0 = Math.max(0, Math.floor((it.ey - hy) / CELL)), y1 = Math.min(gh - 1, Math.floor((it.ey + hy) / CELL));
    for (let gy = y0; gy <= y1; gy++) for (let gx = x0; gx <= x1; gx++) cells[gy * gw + gx].push(j);
  });
  const pairs = new Set<number>();
  const P0 = lastProj[0], P5 = lastProj[5], P8 = lastProj[8], P9 = lastProj[9];
  for (let gy = 0; gy < gh; gy++) for (let gx = 0; gx < gw; gx++) {
    const list = cells[gy * gw + gx];
    if (list.length < 2) continue;
    for (let py = gy * CELL; py < Math.min(h, (gy + 1) * CELL); py++) for (let px = gx * CELL; px < Math.min(w, (gx + 1) * CELL); px++) {
      const u = ((2 * (px + 0.5)) / w - 1 + P8) / P0, v = (1 - (2 * (py + 0.5)) / h + P9) / P5, a = u * u + v * v + 1;
      let best = Infinity, vis = -1;
      for (const j of list) {
        const bb = u * CX[j] + v * CY[j] - CZ[j], disc = bb * bb - a * CC[j];
        if (disc < 0) continue;
        const t = (bb - Math.sqrt(disc)) / a;
        if (t < best) { best = t; vis = j; }
      }
      if (vis < 0) continue;
      for (const j of list) {
        if (j <= vis) continue;
        const bb = u * CX[j] + v * CY[j] - CZ[j];
        if (bb * bb - a * CC[j] >= 0) pairs.add(j * n + vis);
      }
    }
  }
  for (const key of pairs) {
    const A = back[Math.floor(key / n)], B = back[key % n];
    const ab = sub(B.c, A.c), L = Math.hypot(ab[0], ab[1], ab[2]);
    if (L >= A.r + B.r) continue;                             // disjointes à l'arrondi près
    if (L + A.r <= B.r) { A.hide.push(rimPts(A, A).map(proj)); continue; }   // X dans V : jamais vue
    if (L + B.r <= A.r) continue;
    const nn = ab.map((x) => x / L), t = (L * L + A.r * A.r - B.r * B.r) / (2 * L), rad = Math.sqrt(Math.max(0, A.r * A.r - t * t));
    const C = circle([A.c[0] + t * nn[0], A.c[1] + t * nn[1], A.c[2] + t * nn[2]], nn, rad, samples(A, Math.max(rad, A.r)));
    const front = C.map((x) => dot(sub(x, A.c), x) < 0);
    const rim = rimPts(A, A);
    const inB = rim.map((x) => dot(sub(x, B.c), sub(x, B.c)) < B.r * B.r);
    const arc = run(C, front), rarc = run(rim, inB);
    let poly: number[][];
    if (arc.closed) poly = arc.pts;
    else if (rarc.closed) poly = rim;
    else if (!arc.pts.length) continue;
    else {
      const end = arc.pts[arc.pts.length - 1], r0 = rarc.pts[0], r1 = rarc.pts[rarc.pts.length - 1];
      const d2 = (p: number[], q: number[]) => dot(sub(p, q), sub(p, q));
      poly = arc.pts.concat(rarc.pts.length && d2(end, r1) < d2(end, r0) ? rarc.pts.slice().reverse() : rarc.pts);
    }
    A.hide.push(poly.map(proj));
  }
  /** Points d'échantillonnage : environ un par pixel du contour projeté, de 16 à 128 ; la corde
   *  s'écarte alors de l'arc vrai de moins de 0,05 px dès 2,5 px de rayon. */
  function samples(A: Item, rad: number): number {
    return Math.max(16, Math.min(128, Math.ceil((2 * Math.PI * s * rad) / -A.c[2])));
  }
  function rimPts(A: Item, _: Item): number[][] {
    const d2 = dot(A.c, A.c), nn = A.c.map((x) => x / Math.sqrt(d2)), f = 1 - (A.r * A.r) / d2;
    return circle([A.c[0] * f, A.c[1] * f, A.c[2] * f], nn, A.r * Math.sqrt(d2 - A.r * A.r) / Math.sqrt(d2), samples(A, A.r));
  }

  const colors = cur ? ($<HTMLSelectElement>("colour").value === "var" ? varColors() : cur.chain) : new Float32Array();
  const rgb = (i: number) => `rgb(${[0, 1, 2].map((k) => Math.round(255 * Math.min(1, colors[3 * (i - 1) + k]))).join(",")})`;
  const m = last?.manifest;
  const caption = m ? `${m.input.name} · ${m.region} à ${m.params.resolution} pb · alpha ${m.params.alpha} · reconstruite (ShRec3D + MDS), chiralité arbitraire${cur && cur.replica ? ` · réplicat ${cur.replica}` : ""}` : "";
  const f2 = (v: number) => v.toFixed(2);
  const masks = back.filter((b) => b.hide.length).map((b) =>
    `<mask id="m${b.id}" maskUnits="userSpaceOnUse" x="0" y="0" width="${w}" height="${h}"><rect width="${w}" height="${h}" fill="#fff"/>` +
    b.hide.map((pl) => `<path d="M${pl.map((q) => `${f2(q[0])} ${f2(q[1])}`).join("L")}Z" fill="#000"/>`).join("") + `</mask>`).join("");
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}">` + (masks ? `<defs>${masks}</defs>` : "") + `<rect width="100%" height="100%" fill="#ffffff"/>` +
    back.map((b) => `<ellipse cx="${f2(b.ex)}" cy="${f2(b.ey)}" rx="${f2(b.rx)}" ry="${f2(b.ry)}" transform="rotate(${f2(b.rot)} ${f2(b.ex)} ${f2(b.ey)})" fill="${rgb(b.id)}" stroke="#222" stroke-width="0.4"${b.hide.length ? ` mask="url(#m${b.id})"` : ""} data-bead="${b.id}"/>`).join("") +
    `<text x="${w - 8}" y="${h - 10}" text-anchor="end" font-family="sans-serif" font-size="11" fill="#333">${esc(caption)}</text></svg>`;
  return { svg, circles: items.map(({ power: _p, c: _c, r: _r, hide: _h, ...b }) => b) };
}

function save(name: string, blob: Blob): void {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = name;
  a.click();
  exported.push({ name, bytes: blob.size, type: blob.type });
}
const exported: { name: string; bytes: number; type: string }[] = [];

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
  for (const el of [selRes, inRegion, selBal, inAlpha, btnRun, inReplay, $<HTMLInputElement>("ensemble"), $<HTMLInputElement>("seed")]) el.disabled = false;
  if ((r.inspect as Inspect & { converted_from_hic?: boolean }).converted_from_hic) info.innerHTML += `<br><span class="k">.hic converti en .mcool par hic2cool, sur cette machine</span>`;
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
    ${s.observed_pairs < 0.05 ? `<div class="err">Carte trop creuse à cette résolution : ${(s.observed_pairs * 100).toFixed(1)} % des paires de billes ont un contact. La structure tient presque toute à la complétion par les plus courts chemins — choisir une résolution plus grossière.</div>` : ""}
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
    const params = { resolution: Number(selRes.value), region: inRegion.value, alpha: Number(inAlpha.value), balance: selBal.value,
                     ensemble: Number($<HTMLInputElement>("ensemble").value), seed: Number($<HTMLInputElement>("seed").value) };
    const { job } = await api<{ job: string }>("api/run", JSON.stringify({ file: file!.token, params }));
    const st = await follow(job);
    if (st.state !== "fait") throw new Error(st.error ?? st.state);
    last = { session: st.session!, manifest: st.manifest! };
    describe(st.manifest!);
    await show(st.session!);
    await extras(st.session!, st.manifest!);
  });
});
$<HTMLInputElement>("replica").addEventListener("input", (e) => void replica(Number((e.target as HTMLInputElement).value)));
$<HTMLSelectElement>("colour").addEventListener("change", upload);
$<HTMLButtonElement>("ex-svg").addEventListener("click", () => save("structure.svg", new Blob([svgExport().svg], { type: "image/svg+xml" })));
$<HTMLButtonElement>("ex-png").addEventListener("click", () => {
  // Le tampon WebGL n'est garanti qu'au sein de la tâche qui l'a dessiné : on dessine, puis on copie.
  render();
  gl.bindFramebuffer(gl.READ_FRAMEBUFFER, null);
  const px = new Uint8Array(W * H * 4);
  gl.readPixels(0, 0, W, H, gl.RGBA, gl.UNSIGNED_BYTE, px);
  const c2 = document.createElement("canvas");
  c2.width = W;
  c2.height = H;
  const ctx = c2.getContext("2d")!;
  const img = ctx.createImageData(W, H);
  for (let y = 0; y < H; y++) img.data.set(px.subarray((H - 1 - y) * W * 4, (H - y) * W * 4), y * W * 4);
  ctx.putImageData(img, 0, 0);
  c2.toBlob((b) => b && save("structure.png", b), "image/png");
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
    await extras(st.session!, st.manifest!);
    const diff = Object.entries(r.versions_differ).map(([k, [a, b]]) => `${k} ${a} → ${b}`).join(", ");
    verdict.innerHTML = r.identical
      ? `<span class="ok">Reproduit à l'identique : le .g3d a les mêmes octets (sha256 ${r.sha256.slice(0, 16)}…)${r.files_compared.length ? `, et ${r.files_compared.length} fichiers d'accompagnement aussi (${r.files_compared.join(", ")})` : ""}.</span>`
      : `<span class="err">Différent : sha256 ${r.sha256.slice(0, 16)}… au lieu de ${r.expected_sha256.slice(0, 16)}…</span><br><span class="k">${r.files_differ.length ? `fichiers différents : ${esc(r.files_differ.join(", "))} — ` : ""}${esc(r.note ?? "")}${diff ? ` — ${esc(diff)}` : ""}</span>`;
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
  if (b) hover.textContent = `${b.chrom}:${fmt(b.start + 1)}–${fmt(b.end)} · bille ${id}${b.variability > 0 ? ` · variabilité ${fmt(Math.round(b.variability))} nm` : ""}`;
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
  replica: () => cur?.replica ?? 0,
  shown: () => (cur ? Array.from(cur.shown) : []),
  ensembleBytes: () => ensembleBytes,
  statistics: () => stats18,
  charts: () => document.querySelectorAll("#charts svg path").length,
  svg: () => svgExport(),
  exported: () => exported,
  /** Pour la mesure : l'identifiant lu au pixel CSS (x, y) du dernier rendu. */
  /** Les billes que traverse le rayon du pixel (x, y), de la plus proche à la plus lointaine, avec la
   *  distance d'entrée en µm — le même test que le rendu et que l'export SVG. */
  rayAt: (x: number, y: number) => {
    if (!cur) return [];
    const view = viewMatrix(orbit), w = canvas.clientWidth, h = canvas.clientHeight;
    const u = ((2 * x) / w - 1 + lastProj[8]) / lastProj[0], v = (1 - (2 * y) / h + lastProj[9]) / lastProj[5], a = u * u + v * v + 1;
    const out: { id: number; t: number; power: number }[] = [];
    for (let i = 0; i < cur.ids.length; i++) {
      const p = [cur.shown[3 * i], cur.shown[3 * i + 1], cur.shown[3 * i + 2]];
      const e = [0, 1, 2].map((k) => view[k] * p[0] + view[4 + k] * p[1] + view[8 + k] * p[2] + view[12 + k]);
      const cc = e[0] * e[0] + e[1] * e[1] + e[2] * e[2] - cur.radii[i] * cur.radii[i], b = u * e[0] + v * e[1] - e[2], disc = b * b - a * cc;
      if (disc >= 0) out.push({ id: i + 1, t: ((b - Math.sqrt(disc)) / a) * Math.sqrt(a), power: cc });
    }
    return out.sort((p, q) => p.t - q.t);
  },
  idAt: (x: number, y: number) => { render(); const dpr = W / canvas.clientWidth; return picker.pick(x * dpr, y * dpr); },
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
