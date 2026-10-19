/** La plongée : du noyau entier jusqu'à une base de TP53, sans coupure (semaines 10 à 12).
 *
 *  Deux niveaux du même `.g3d` — l'aperçu à 3 Mb, puis le noyau à 750 kb —, fondus l'un
 *  dans l'autre par tramage complémentaire ; les 45 autres copies qui s'effacent autour de
 *  chr7:a ; les couleurs de territoire qui glissent vers celles de compartiment ; une
 *  occlusion ambiante qui rend le volume lisible. Tout dépend de la caméra (`scene.ts`).
 *
 *  Paramètres d'URL :
 *    f=…        fichier .g3d (défaut data/gm12878.g3d)
 *    bench=N    plongée scriptée de N images, mesurée, puis résultat publié dans window.__geno
 *    hard       témoin : transitions sèches au lieu de fondus (la mesure doit y voir des sauts)
 *    w, h       taille du canevas en mode mesure (défaut 960 × 540)
 *    view=…     une seule image, à une étape de la plongée : final, loop, anchor ; et, en
 *               semaine 12, fibre, nuc, helix, base
 *    legibility mesure la lisibilité de la boucle visée ; noisolate : le témoin sans isolement
 *    junction   le raccord S8 seul, sur la région de la cible (chr7 ou chr17)
 *    target=…   tp53 (défaut, semaine 12) ou actb (la plongée des semaines 10 et 11)
 *    stream=S   la plongée jouée en temps réel en S secondes, caméra réglée sur l'horloge, et
 *               le compte des images qui montrent un chargement ; noprefetch : le témoin
 *    fidelity   les lettres des barreaux visibles, relues sur les pixels de la vue finale
 *    labels     dessine la surimpression des lettres même en mode mesure (pour l'enregistrement)
 *    record=N   N images de la plongée, rendues une à une sur demande (`tools/record.mjs`)
 *    fiche=…    fichier de fiches (défaut : celui du .g3d, en .fiche) — semaine 13
 *    pick       mode interactif instrumenté pour `tools/pick.mjs` : chaque clic est chronométré
 *               (pixel, fiche, DOM, image peinte) ; nofichepre : le témoin sans préchargement
 */

import { FicheReader, HttpSource } from "@geno/g3d";
import { BASE_COLORS, baseOf, dive, type Orbit, rotate, zoom } from "@geno/viewer";
import { divePath, runBench } from "./bench.js";
import { type DiveScene, FOV, openDive, type Selection, type TargetName } from "./dive.js";
import { compose, EVIDENCE_COLORS, FichePanel } from "./fiche.js";
import { IDS } from "./fiber.js";
import { measureJunction, measureLoop } from "./legibility.js";

const q = new URLSearchParams(location.search);
const FILE = q.get("f") ?? "data/gm12878.g3d";
const BENCH = Number(q.get("bench") ?? 0);
const HARD = q.has("hard");
const VIEW = q.get("view") ?? (q.has("legibility") ? "loop" : q.has("fidelity") ? "base" : null);
const TARGET = (q.get("target") ?? "tp53") as TargetName;
const STREAM = Number(q.get("stream") ?? 0);
const RECORD = Number(q.get("record") ?? 0);
const FICHE = q.get("fiche") ?? FILE.replace(/\.g3d$/, ".fiche");
const PICK = q.has("pick");

declare global {
  interface Window {
    __geno?: Record<string, unknown>;
  }
}

const style = document.createElement("style");
style.textContent = `
  html, body { height: 100%; overflow: hidden; }
  #hud { position: fixed; left: 12px; top: 10px; max-width: min(460px, calc(100vw - 24px));
         background: rgba(14,18,24,.82); padding: 10px 12px; border-radius: 8px; font-size: 12px; }
  #hud h1 { font-size: 13px; margin: 0 0 6px; font-weight: 600; }
  #hud button { font: inherit; margin: 6px 6px 0 0; padding: 4px 10px; border-radius: 6px;
                border: 1px solid #3a4658; background: #1b2330; color: inherit; cursor: pointer; }
  #hud .k { color: #8b98aa; }
  .sw { display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin: 0 4px -1px 8px; }
  #info { margin-top: 6px; white-space: pre-wrap; }
  #evidence { margin-top: 6px; color: #c9d3df; }
  canvas { touch-action: none; }
  #overlay { position: fixed; left: 0; top: 0; pointer-events: none; }
`;
document.head.appendChild(style);

const canvas = document.createElement("canvas");
document.body.appendChild(canvas);
const overlay = document.createElement("canvas");
overlay.id = "overlay";
document.body.appendChild(overlay);

function fit(): [number, number] {
  if (BENCH || VIEW || STREAM || RECORD) {
    const w = Number(q.get("w") ?? 960);
    const h = Number(q.get("h") ?? 540);
    canvas.width = w;
    canvas.height = h;
    canvas.style.width = `${w}px`;
    canvas.style.height = `${h}px`;
    sizeOverlay(w, h, w, h);
    return [w, h];
  }
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  canvas.width = Math.round(window.innerWidth * dpr);
  canvas.height = Math.round(window.innerHeight * dpr);
  canvas.style.width = `${window.innerWidth}px`;
  canvas.style.height = `${window.innerHeight}px`;
  sizeOverlay(canvas.width, canvas.height, window.innerWidth, window.innerHeight);
  return [canvas.width, canvas.height];
}

function sizeOverlay(w: number, h: number, cssW: number, cssH: number): void {
  overlay.width = w;
  overlay.height = h;
  overlay.style.width = `${cssW}px`;
  overlay.style.height = `${cssH}px`;
}

/** Les lettres des bases et l'étiquette de la base visée, par-dessus le rendu WebGL. */
function drawLabels(scene: DiveScene): void {
  const ctx = overlay.getContext("2d");
  if (!ctx) return;
  ctx.clearRect(0, 0, overlay.width, overlay.height);
  const scale = overlay.width / Math.max(1, overlay.clientWidth || overlay.width);
  // L'étiquette de la base visée en dernier : au-dessus des lettres.
  for (const l of scene.labels().sort((a, b) => Number(!!a.target) - Number(!!b.target))) {
    ctx.globalAlpha = l.alpha;
    if (l.target) {
      const lines = l.text.split("\n");
      const fs = Math.round(13 * scale);
      ctx.font = `600 ${fs}px system-ui, sans-serif`;
      const w = Math.max(...lines.map((t) => ctx.measureText(t).width)) + 16 * scale;
      const h = lines.length * fs * 1.35 + 10 * scale;
      // Dans le coin haut droit, reliée à la base par un trait : posée à côté de la base, elle
      // couvrait les barreaux voisins à 4 nm de champ.
      const x = overlay.width - w - 16 * scale;
      const y = 16 * scale;
      ctx.strokeStyle = "rgba(255,255,255,0.9)";
      ctx.lineWidth = 1.5 * scale;
      ctx.beginPath();
      ctx.moveTo(l.x, l.y);
      ctx.lineTo(x + w / 2, y + h);
      ctx.stroke();
      ctx.beginPath();
      ctx.arc(l.x, l.y, 4 * scale, 0, 2 * Math.PI);
      ctx.stroke();
      ctx.fillStyle = "rgba(14,18,24,0.88)";
      ctx.fillRect(x, y, w, h);
      ctx.fillStyle = "#fff";
      lines.forEach((t, i) => {
        ctx.font = `${i === 0 ? 600 : 400} ${fs}px system-ui, sans-serif`;
        ctx.fillText(t, x + 8 * scale, y + 6 * scale + (i + 0.85) * fs * 1.35);
      });
      continue;
    }
    const [r, g, b] = l.color.map((c) => Math.round(c * 255));
    ctx.font = `700 ${Math.round(l.size)}px ui-monospace, monospace`;
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.lineWidth = Math.max(2, l.size / 6);
    ctx.strokeStyle = "rgba(10,12,16,0.85)";
    ctx.strokeText(l.text, l.x, l.y);
    ctx.fillStyle = `rgb(${r},${g},${b})`;
    ctx.fillText(l.text, l.x, l.y);
    ctx.textAlign = "start";
    ctx.textBaseline = "alphabetic";
  }
  ctx.globalAlpha = 1;
}

/** La légende de l'enregistrement : ce qu'on regarde, à quelle échelle, et d'où ça vient. */
function drawCaption(scene: DiveScene): void {
  const ctx = overlay.getContext("2d");
  if (!ctx) return;
  const w = scene.last.widthUm;
  const field = w >= 1 ? `${w.toLocaleString("fr", { maximumFractionDigits: 1 })} µm` : `${(w * 1000).toLocaleString("fr", { maximumFractionDigits: w < 0.01 ? 1 : 0 })} nm`;
  const k = overlay.height / 720;
  const lines = [`Geno · GM12878 · du noyau à une base de TP53 · champ ${field}`, scene.evidence()];
  ctx.font = `${Math.round(14 * k)}px system-ui, sans-serif`;
  const bw = Math.max(...lines.map((t) => ctx.measureText(t).width)) + 20 * k;
  const bh = lines.length * 20 * k + 12 * k;
  ctx.fillStyle = "rgba(14,18,24,0.78)";
  ctx.fillRect(14 * k, overlay.height - bh - 14 * k, bw, bh);
  lines.forEach((t, i) => {
    ctx.fillStyle = i ? "#c9d3df" : "#ffffff";
    ctx.font = `${i ? 400 : 600} ${Math.round(14 * k)}px system-ui, sans-serif`;
    ctx.fillText(t, 24 * k, overlay.height - bh - 14 * k + (i + 1) * 20 * k);
  });
}

/** Fidélité : chaque pixel de barreau de la vue finale, classé par teinte, puis un vote par
 *  base. Rend { rang: lettre } pour chaque brin ; l'outil compare à la séquence lue par
 *  Python dans le .2bit, par un autre chemin que celui du navigateur. */
function readLetters(scene: DiveScene): Record<string, string> {
  const gl = scene.gl;
  const [W, H] = scene.size();
  const ids = scene.ids();
  const rgb = new Uint8Array(W * H * 4);
  gl.bindFramebuffer(gl.READ_FRAMEBUFFER, null);
  gl.readPixels(0, 0, W, H, gl.RGBA, gl.UNSIGNED_BYTE, rgb);
  const hue = (r: number, g: number, b: number) => {
    const mx = Math.max(r, g, b), mn = Math.min(r, g, b), d = mx - mn;
    if (!d || d / mx < 0.35) return -1;
    let h = mx === r ? ((g - b) / d) % 6 : mx === g ? (b - r) / d + 2 : (r - g) / d + 4;
    h *= 60;
    return h < 0 ? h + 360 : h;
  };
  const ref = Object.fromEntries((["A", "C", "G", "T"] as const).map((k) => [k, hue(...(BASE_COLORS[k].map((x) => x * 255) as [number, number, number]))]));
  const votes = new Map<string, Record<string, number>>();
  for (let i = 0; i < ids.length; i++) {
    const id = ids[i];
    if (id < IDS.base || id === 0xffffffff) continue;
    const h = hue(rgb[4 * i], rgb[4 * i + 1], rgb[4 * i + 2]);
    if (h < 0) continue;
    let best = "N", gap = Infinity;
    for (const [k, v] of Object.entries(ref)) {
      const dh = Math.min(Math.abs(h - v), 360 - Math.abs(h - v));
      if (dh < gap) [best, gap] = [k, dh];
    }
    const { rank, strand } = baseOf(id - IDS.base);
    const key = `${rank}:${strand}`;
    const v = votes.get(key) ?? {};
    v[best] = (v[best] ?? 0) + 1;
    votes.set(key, v);
  }
  const out: Record<string, string> = {};
  for (const [k, v] of votes) out[k] = Object.entries(v).sort((a, b) => b[1] - a[1])[0][0];
  return out;
}

function publish(r: Record<string, unknown>): void {
  window.__geno = r;
}

async function main(): Promise<void> {
  fit();
  const scene = await openDive(canvas, FILE, { hard: HARD, noIsolate: q.has("noisolate"), target: TARGET, noPrefetch: q.has("noprefetch") });

  if (q.has("junction")) {
    publish({ ok: true, junction: await measureJunction(scene) });
    return;
  }

  if (q.has("legibility")) {
    const loop = measureLoop(scene);
    publish({ ...loop, isolate: !q.has("noisolate"), junction: q.has("noisolate") ? null : await measureJunction(scene) });
    return;
  }

  if (VIEW) {
    const st = scene.stops;
    const o = ({ final: scene.final, loop: scene.loopView, anchor: scene.anchorView, fibre: st[2], nuc: st[3], helix: st[4], base: st[5] } as Record<string, Orbit | null | undefined>)[VIEW];
    scene.arrival = 1;
    await scene.preload();
    scene.growBudgetMs = Infinity;
    if (o) {
      scene.render(o);
      scene.render(o);      // la seconde image a sa fenêtre d'hélice
    }
    if (q.has("labels")) drawLabels(scene);
    publish({ ok: !!o, view: VIEW, state: scene.last, focus: scene.fiber?.focus ?? null,
              arcs: scene.fiber?.arcs.length ?? 0, held: scene.fiber?.held.length ?? 0,
              nucleosomes: scene.fiber?.nucleosomeCount() ?? 0, target: scene.target,
              region_start: scene.sequence?.regionStart ?? null, target_rank: scene.sequence?.targetRank ?? null,
              letters: q.has("fidelity") ? readLetters(scene) : null, evidence: scene.evidence(),
              loading: scene.sequence?.stats ?? null });
    return;
  }

  if (RECORD) {
    // Image par image, sur demande de l'outil : `render(k)` dessine l'image k et sa
    // surimpression ; l'outil capture l'écran et passe à la suivante.
    const path = divePath(scene);
    scene.arrival = 1;
    await scene.preload();
    scene.growBudgetMs = Infinity;
    (window as unknown as { __record: unknown }).__record = {
      frames: RECORD,
      render(k: number): Record<string, unknown> {
        const o = path.at((Math.min(RECORD - 1, k) / (RECORD - 1)) * path.length);
        scene.render(o);
        drawLabels(scene);
        drawCaption(scene);
        return { ...scene.last };
      },
    };
    publish({ ok: true, record: RECORD, stops: scene.stops.length });
    return;
  }

  if (STREAM) {
    // Temps réel : la caméra suit l'horloge, le réseau fait ce qu'il peut, et chaque image
    // dit si elle montre un chargement. Le chronomètre part au premier rendu de la plongée.
    const path = divePath(scene);
    scene.arrival = 1;
    const frames: { t: number; w: number; helix: number; loading: boolean }[] = [];
    const t0 = performance.now();
    await new Promise<void>((done) => {
      const tick = () => {
        const t = (performance.now() - t0) / 1000;
        const k = Math.min(1, t / STREAM);
        const before = scene.sequence?.stats.loading ?? 0;
        scene.render(path.at(k * path.length));
        frames.push({ t, w: scene.last.widthUm, helix: scene.last.helix, loading: (scene.sequence?.stats.loading ?? 0) > before });
        if (k < 1) requestAnimationFrame(tick);
        else done();
      };
      requestAnimationFrame(tick);
    });
    publish({ ok: true, seconds: STREAM, prefetch: !q.has("noprefetch"), frames, stats: scene.sequence?.stats ?? null,
              sequence_bytes: scene.sequence?.bytes ?? 0, requests: scene.g3d.requests, first_render_ms: scene.firstRenderMs });
    return;
  }

  if (BENCH) {
    const r = await runBench(scene, BENCH, { continuity: true });
    publish({
      ok: r.diffs.length === BENCH && Number.isFinite(r.diffs[0]),
      hard: HARD,
      frames: r.frames,
      size: r.size,
      gpu: r.gpu,
      path_length: r.pathLength,
      first_render_ms: scene.firstRenderMs,
      instances: r.instances,
      diffs: r.diffs,
      frame_ms: r.frameMs,
      raf_ms: r.rafMs,
      states: r.states,
    });
    return;
  }

  // Mode interactif : glisser tourne, la molette zoome, un clic sans glisser identifie.
  let orbit: Orbit = { ...scene.overview };
  let anim: { path: ReturnType<typeof dive>; t0: number; duration: number } | null = null;
  let dirty = true;
  const arrivedAt = performance.now();
  const frameTimes: number[] = [];
  let lastFrame = performance.now();

  const hud = document.createElement("div");
  hud.id = "hud";
  hud.innerHTML = `<h1>Geno — du noyau à une base de TP53</h1>
    <div class="k">GM12878 · modèle simulé (médoïde de 200 repliements) · compartiments appelés sur la
    carte de contacts du modèle, pas sur du Hi-C mesuré · séquence hg38 de référence</div>
    <div><span class="sw" style="background:rgb(224,107,51)"></span>A
         <span class="sw" style="background:rgb(56,107,219)"></span>B
         <span class="sw" style="background:rgb(128,133,140)"></span>sans compartiment (Xi)</div>`;
  document.body.appendChild(hud);
  const go = (label: string, to: () => Orbit) => {
    const b = document.createElement("button");
    b.textContent = label;
    b.onclick = () => {
      const path = dive(orbit, to(), FOV);
      anim = { path, t0: performance.now(), duration: 900 * Math.max(0.5, path.length) };
    };
    hud.appendChild(b);
  };
  if (scene.stops.length > 4) {
    // La plongée entière, étape par étape, comme la mesure la joue.
    const b = document.createElement("button");
    b.textContent = `Plonger vers ${scene.target?.name ?? "la base"}`;
    b.onclick = () => {
      const path = divePath({ ...scene, stops: [orbit, ...scene.stops.slice(1)] });
      anim = { path, t0: performance.now(), duration: 900 * Math.max(0.5, path.length) };
    };
    hud.appendChild(b);
  } else {
    go("Plonger vers ACTB", () => scene.final);
  }
  go("Remonter", () => scene.overview);
  // Semaine 13 : la fiche. Le fichier s'ouvre après le premier rendu et se précharge en
  // arrière-plan, le chromosome visé d'abord ; un clic n'attend jamais ce préchargement.
  const panel = new FichePanel(document);
  let reader: FicheReader | null = null;
  const picks: Record<string, unknown>[] = [];
  const prefetch: Record<string, number> = { started: performance.now() };
  let marker: number[] | null = null;
  const focusChrom = scene.region?.chrom ?? "chr17";
  const opened = FicheReader.open(new HttpSource(FICHE)).then(async (r) => {
    reader = r;
    if (q.has("nofichepre")) return r;
    const order = [
      { level: "noyau", first: [focusChrom] },
      { level: scene.g3d.header.first, first: [focusChrom] },
      ...r.header.levels.filter((l) => l.fine).map((l) => ({ level: l.name })),
    ];
    await r.prefetch(order, (step) => (prefetch[step] = performance.now() - prefetch.started));
    prefetch.all = performance.now() - prefetch.started;
    prefetch.bytes = r.bytes;
    return r;
  }).catch((e) => {
    info.textContent = `fiches indisponibles : ${e}`;
    return null;
  });

  /** Un clic : le pixel, la fiche, le DOM — chronométrés séparément, puis l'image peinte. */
  const onClick = (x: number, y: number) => {
    const t0 = performance.now();
    const sel: Selection | null = scene.select(x, y);
    const t1 = performance.now();
    if (!sel) {
      panel.hide();
      marker = null;
      picks.push({ hit: false, pickMs: t1 - t0 });
      return;
    }
    marker = sel.world;
    const r = reader;
    const card = r ? r.card(sel.level, sel.bead) : null;
    const t2 = performance.now();
    const c = compose(sel, card, (r ?? { header: null }).header as never, r, scene.sequence?.reader ?? null);
    if (r) panel.show(c, r.header, !card);
    else info.textContent = `${c.title}\n${c.subtitle}\n(fiches en cours d'ouverture)`;
    const t3 = performance.now();
    const entry: Record<string, unknown> = {
      hit: true, what: sel.what, level: sel.level, bead: sel.bead, pos: sel.pos, ready: !!card,
      pickMs: t1 - t0, cardMs: t2 - t1, domMs: t3 - t2, expected: c.expected, valued: c.valued, missing: c.missing,
    };
    picks.push(entry);
    requestAnimationFrame(() => {
      entry.frameMs = performance.now() - t0;
    });
    if (!card) {
      // À froid : les coordonnées sont déjà là, le reste arrive avec son bloc.
      opened.then(async (rr) => {
        if (!rr) return;
        const full = await rr.cardAsync(sel.level, sel.bead);
        const c2 = compose(sel, full, rr.header, rr, scene.sequence?.reader ?? null);
        panel.show(c2, rr.header, false);
        entry.waitMs = performance.now() - t0;
        entry.expected = c2.expected;
        entry.valued = c2.valued;
        entry.missing = c2.missing;
      });
    }
  };

  const status = document.createElement("div");
  const evidence = document.createElement("div");
  evidence.id = "evidence";
  const info = document.createElement("div");
  info.id = "info";
  hud.appendChild(status);
  hud.appendChild(evidence);
  hud.appendChild(info);

  let drag: { x: number; y: number; moved: number } | null = null;
  canvas.addEventListener("pointerdown", (e) => {
    drag = { x: e.clientX, y: e.clientY, moved: 0 };
    canvas.setPointerCapture(e.pointerId);
  });
  canvas.addEventListener("pointermove", (e) => {
    if (!drag) return;
    const dx = e.clientX - drag.x;
    const dy = e.clientY - drag.y;
    drag = { x: e.clientX, y: e.clientY, moved: drag.moved + Math.abs(dx) + Math.abs(dy) };
    orbit = rotate(orbit, -dx * 0.006, dy * 0.006);
    anim = null;
    dirty = true;
  });
  canvas.addEventListener("pointerup", async (e) => {
    const wasClick = drag && drag.moved < 4;
    drag = null;
    if (!wasClick) return;
    const dpr = canvas.width / canvas.clientWidth;
    onClick(e.offsetX * dpr, e.offsetY * dpr);
  });
  canvas.addEventListener(
    "wheel",
    (e) => {
      e.preventDefault();
      orbit = zoom(orbit, Math.exp(e.deltaY * 0.0015), 0.003, 150);
      anim = null;
      dirty = true;
    },
    { passive: false },
  );
  window.addEventListener("resize", () => {
    scene.resize(...fit());
    dirty = true;
  });

  const loop = (now: number) => {
    if (scene.fine && scene.arrival < 1) {
      scene.arrival = Math.min(1, (now - arrivedAt) / 400);   // seul fondu réglé par le temps
      dirty = true;
    }
    if (anim) {
      const k = Math.min(1, (now - anim.t0) / anim.duration);
      orbit = anim.path.at(k * anim.path.length);
      if (k >= 1) anim = null;
      dirty = true;
    }
    // Les nucléosomes se génèrent par tranches et la séquence arrive par le réseau : tant que
    // la fenêtre n'est pas complète, l'image suivante peut montrer plus.
    if (scene.last.loop > 0 || scene.last.helix > 0) dirty = true;
    if (dirty || marker) {
      scene.render(orbit);
      drawLabels(scene);
      if (marker) drawMarker(scene, marker);
      evidence.textContent = scene.evidence();
      dirty = false;
      frameTimes.push(now - lastFrame);
      if (frameTimes.length > 60) frameTimes.shift();
      const sorted = frameTimes.slice().sort((a, b) => a - b);
      const s = scene.last;
      status.innerHTML =
        `<span class="k">champ</span> ${s.widthUm >= 1 ? `${s.widthUm.toFixed(1)} µm` : `${(s.widthUm * 1000).toFixed(s.widthUm < 0.01 ? 1 : 0)} nm`} · ` +
        `<span class="k">750 kb</span> ${(s.lod * 100).toFixed(0)} % · ` +
        `<span class="k">autres copies</span> ${((1 - s.focus) * 100).toFixed(0)} % · ` +
        `<span class="k">compartiments</span> ${(s.mix * 100).toFixed(0)} % · ` +
        `<span class="k">hélice</span> ${(s.helix * 100).toFixed(0)} % · ` +
        `<span class="k">image</span> ${sorted[sorted.length >> 1]?.toFixed(1) ?? "–"} ms<br>` +
        `<span class="k">${scene.gpu}</span>`;
    }
    lastFrame = now;
    requestAnimationFrame(loop);
  };
  requestAnimationFrame(loop);

  if (PICK) {
    // Pour `tools/pick.mjs` : aller à une étape de la plongée, tirer des pixels qui portent
    // un objet, relire le journal des clics et l'état du préchargement.
    scene.arrival = 1;
    await scene.preload();
    scene.growBudgetMs = Infinity;
    (window as unknown as { __pick: unknown }).__pick = {
      stops: scene.stops.length,
      goto(i: number) {
        orbit = scene.stops[Math.min(i, scene.stops.length - 1)];
        anim = null;
        scene.render(orbit);
        scene.render(orbit);
        dirty = true;
        return { ...scene.last };
      },
      /** `n` pixels CSS au hasard (graine fixe) parmi ceux qui portent un identifiant. */
      hits(n: number, seed: number) {
        const ids = scene.ids();
        const [W, H] = scene.size();
        const dpr = canvas.width / canvas.clientWidth;
        let st = seed >>> 0;
        const rnd = () => ((st = (st * 1664525 + 1013904223) >>> 0) / 2 ** 32);
        const out: [number, number][] = [];
        for (let tries = 0; out.length < n && tries < 200 * n; tries++) {
          const x = Math.floor(rnd() * W), y = Math.floor(rnd() * H);
          const cx = (x + 0.5) / dpr, cy = (y + 0.5) / dpr;
          // Un pixel sous le HUD ou sous la fiche ouverte n'est pas cliquable : on l'écarte.
          if (ids[(H - 1 - y) * W + x] && document.elementFromPoint(cx, cy) === canvas) out.push([cx, cy]);
        }
        return out;
      },
      log: picks,
      prefetch,
      reader: () => (reader ? { bytes: reader.bytes, requests: reader.requests, tables: reader.tablesBytes(),
        levels: Object.fromEntries(reader.header.levels.map((l) => [l.name, reader!.levelBytes(l.name)])) } : null),
    };
    publish({ ok: true, pick: true, stops: scene.stops.length, gpu: scene.gpu });
  }
}

/** Le point cliqué, repéré sur la surimpression : un cercle et une croix fine. */
function drawMarker(scene: DiveScene, p: number[]): void {
  const ctx = overlay.getContext("2d")!;
  const [x, y] = scene.project(p);
  const r = 9 * (overlay.width / overlay.clientWidth || 1);
  ctx.save();
  ctx.strokeStyle = EVIDENCE_COLORS.autre;
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.arc(x, y, r, 0, 2 * Math.PI);
  ctx.moveTo(x - 1.8 * r, y);
  ctx.lineTo(x - 1.2 * r, y);
  ctx.moveTo(x + 1.2 * r, y);
  ctx.lineTo(x + 1.8 * r, y);
  ctx.stroke();
  ctx.restore();
}

main().catch((e) => publish({ ok: false, error: String(e) }));
