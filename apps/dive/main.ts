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
 *    search=…   index de recherche (défaut : celui du .g3d, en .search) — semaine 16
 *    nav        mode interactif instrumenté pour `tools/search.mjs` (recherche, URL, signets)
 *
 *  Semaine 16 : le fragment (#v=1&o=…) porte la vue — orbite, sélection, surimpression, fiche —
 *  et la page le rejoue au chargement (`url.ts`).
 */

import { FicheReader, HttpSource, SearchIndex, type SearchHit as Hit } from "@geno/g3d";
import { BASE_COLORS, baseOf, dive, type Orbit, rotate, zoom } from "@geno/viewer";
import { divePath, runBench } from "./bench.js";
import { type DiveScene, FOV, openDive, type Region, type Selection, type TargetName } from "./dive.js";
import { compose, EVIDENCE_COLORS, FichePanel } from "./fiche.js";
import { IDS } from "./fiber.js";
import { measureJunction, measureLoop } from "./legibility.js";
import { depthText, type Destination, destination } from "./navigate.js";
import { type Bookmark, type CardState, decode, encode, loadBookmarks, saveBookmarks, type ViewState } from "./url.js";

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
const NAV = q.has("nav");
const SEARCH = q.get("search") ?? FILE.replace(/\.g3d$/, ".search");

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
  #search { display: flex; gap: 6px; margin-top: 8px; }
  #search input { flex: 1; min-width: 0; font: inherit; padding: 4px 8px; border-radius: 6px;
                  border: 1px solid #3a4658; background: #0e131b; color: inherit; }
  #search button { margin: 0; }
  #found { margin-top: 4px; color: #c9d3df; white-space: pre-wrap; }
  #found button, #marks button { margin: 4px 4px 0 0; padding: 2px 8px; }
  #marks a { color: #9cc3ff; margin-right: 4px; }
`;
document.head.appendChild(style);

const canvas = document.createElement("canvas");
document.body.appendChild(canvas);
const overlay = document.createElement("canvas");
overlay.id = "overlay";
document.body.appendChild(overlay);
// Semaine 16 : la surimpression est rastérisée par le processeur dès le départ. Chrome bascule
// sinon un canevas 2D du GPU au processeur après la première relecture de ses pixels, et ses
// arcs et ses lettres ne se dessinent plus au même pixel près : une vue relue depuis son URL
// différait de l'originale pour cette seule raison (8 vues sur 30).
overlay.getContext("2d", { willReadFrequently: true });

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
/** La fiche est ouverte : l'étiquette de la base visée, qui en dirait une partie, se tait —
 *  dans le coin haut droit, elle passait sous la fiche, et à sa gauche, sous le HUD. */
let ficheOpen = false;

function drawLabels(scene: DiveScene): void {
  const ctx = overlay.getContext("2d");
  if (!ctx) return;
  ctx.clearRect(0, 0, overlay.width, overlay.height);
  const scale = overlay.width / Math.max(1, overlay.clientWidth || overlay.width);
  // L'étiquette de la base visée en dernier : au-dessus des lettres.
  for (const l of scene.labels().sort((a, b) => Number(!!a.target) - Number(!!b.target))) {
    ctx.globalAlpha = l.alpha;
    if (l.target) {
      if (ficheOpen) continue;
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
    if (id < IDS.base || id >= IDS.linker || id === 0xffffffff) continue;   // l'ADN de liaison est au-dessus des bases
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
  panel.onHide = () => {
    ficheOpen = false;
    marker = null;
    cardState = null;
    drawLabels(scene);
  };
  let reader: FicheReader | null = null;
  const picks: Record<string, unknown>[] = [];
  const prefetch: Record<string, number> = { started: performance.now() };
  let marker: number[] | null = null;
  /** Semaine 16 : la fiche ouverte (pour l'URL), et ce qui est surligné. */
  let cardState: ViewState["card"] = null;
  let selRegions: Region[] | null = null;
  const focusChrom = scene.region?.chrom ?? "chr17";
  // Deux promesses : le fichier ouvert (en-tête lu) — ce qu'un clic à froid attend —, et le
  // préchargement, qui dure, lui, jusqu'à deux minutes et demie sous Slow 4G.
  const opened = FicheReader.open(new HttpSource(FICHE)).then((r) => {
    reader = r;
    return r;
  }).catch((e) => {
    info.textContent = `fiches indisponibles : ${e}`;
    return null;
  });
  void opened.then(async (r) => {
    if (!r || q.has("nofichepre")) return;
    const order = [
      { level: "noyau", first: [focusChrom] },
      { level: scene.g3d.header.first, first: [focusChrom] },
      ...r.header.levels.filter((l) => l.fine).map((l) => ({ level: l.name })),
    ];
    await r.prefetch(order, (step) => (prefetch[step] = performance.now() - prefetch.started));
    prefetch.all = performance.now() - prefetch.started;
    prefetch.bytes = r.bytes;
  });

  /** Un clic : le pixel, la fiche, le DOM — chronométrés séparément, puis l'image peinte. */
  const onClick = (x: number, y: number) => {
    const t0 = performance.now();
    const sel: Selection | null = scene.select(x, y);
    const t1 = performance.now();
    if (!sel) {
      panel.hide();
      marker = null;
      cardState = null;
      ficheOpen = false;
      drawLabels(scene);
      picks.push({ hit: false, pickMs: t1 - t0 });
      return;
    }
    marker = sel.world;
    cardState = { level: sel.level, bead: sel.bead, pos: sel.pos, what: sel.what, strand: sel.strand,
                  precision: sel.precision, detail: sel.detail ?? "" };
    const r = reader;
    const card = r ? r.card(sel.level, sel.bead) : null;
    const t2 = performance.now();
    const c = compose(sel, card, (r ?? { header: null }).header as never, r, scene.sequence?.reader ?? null);
    if (r) panel.show(c, r.header, !card);
    ficheOpen = !!r;
    // Le repère du point cliqué : la surimpression seule, pas une image WebGL de plus.
    drawLabels(scene);
    if (marker) drawMarker(scene, marker);
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

  // Semaine 15 : la surimpression clinique, bille par bille, depuis les fiches — ClinVar P/LP,
  // associations GWAS, contrainte gnomAD (LOEUF le plus bas des gènes de la bille).
  const layer = document.createElement("select");
  layer.innerHTML = `<option value="">surimpression : aucune</option><option value="clinvar">ClinVar P/LP par bille</option>
    <option value="gwas">associations GWAS par bille</option><option value="loeuf">contrainte gnomAD (LOEUF minimal)</option>`;
  layer.setAttribute("aria-label", "Surimpression clinique sur la structure");
  layer.style.cssText = "font: inherit; margin: 6px 6px 0 0; padding: 3px 6px; border-radius: 6px; border: 1px solid #3a4658; background: #1b2330; color: inherit;";
  hud.appendChild(layer);
  const legend = document.createElement("div");
  legend.className = "k";
  hud.appendChild(legend);
  const overlay2: Record<string, Record<string, Float64Array>> = {};
  const clinicalRaw = async (level: string, kind: string): Promise<Float64Array> => {
    overlay2[level] ??= {};
    if (overlay2[level][kind]) return overlay2[level][kind];
    const r = (await opened)!;
    const n = r.level(level)!.n_beads;
    const out = new Float64Array(n);
    for (let k = 1; k <= n; k++) {
      const c = r.card(level, k) ?? (await r.cardAsync(level, k));
      if (kind === "clinvar") out[k - 1] = c.clinvar?.v?.plp_total ?? NaN;
      else if (kind === "gwas") out[k - 1] = c.gwas?.v?.n ?? NaN;
      else {
        const v = Object.values((c.contrainte?.v ?? {}) as Record<string, { loeuf: number | null }>)
          .map((x) => x.loeuf).filter((x): x is number => typeof x === "number");
        out[k - 1] = v.length ? Math.min(...v) : NaN;
      }
    }
    overlay2[level][kind] = out;
    return out;
  };
  const toUnit = (kind: string, raw: Float64Array): Float32Array => {
    const fin = Array.from(raw).filter(Number.isFinite).sort((a, b) => a - b);
    const q99 = fin.length ? Math.max(1, fin[Math.floor(0.99 * (fin.length - 1))]) : 1;
    return Float32Array.from(raw, (x) => !Number.isFinite(x) ? NaN
      : kind === "loeuf" ? Math.max(0, Math.min(1, 1 - x / 1.2)) : Math.sqrt(Math.min(1, x / q99)));
  };
  layer.onchange = async () => {
    const kind = layer.value;
    const t0 = performance.now();
    for (const lv of [scene.g3d.header.first, "noyau"]) {
      scene.paint(lv, kind ? toUnit(kind, await clinicalRaw(lv, kind)) : null);
    }
    legend.textContent = !kind ? "" : kind === "loeuf"
      ? "foncé : gènes contraints (LOEUF bas) · gris : pas de gène codant noté — gnomAD v4.1, annotation"
      : `foncé : beaucoup ${kind === "clinvar" ? "de variants ClinVar pathogènes ou probablement pathogènes" : "d'associations GWAS"} dans l'intervalle de la bille (racine, plafond au 99e centile) — annotation`;
    paintMs = performance.now() - t0;
    dirty = true;
  };
  let paintMs = 0;

  // Semaine 16 : la recherche — un gène, une position, un rsID —, le vol vers la cible, la vue
  // dans l'URL et les signets.
  const form = document.createElement("form");
  form.id = "search";
  form.setAttribute("role", "search");
  const input = document.createElement("input");
  input.type = "search";
  input.placeholder = "BRCA1, chr17:43,044,295, rs334";
  input.setAttribute("aria-label", "Rechercher un gène, une position hg38 ou un rsID");
  const submit = document.createElement("button");
  submit.type = "submit";
  submit.textContent = "Aller";
  form.append(input, submit);
  hud.appendChild(form);
  const found = document.createElement("div");
  found.id = "found";
  found.setAttribute("aria-live", "polite");
  hud.appendChild(found);
  const marks = document.createElement("div");
  marks.id = "marks";
  hud.appendChild(marks);

  const index = SearchIndex.open(new HttpSource(SEARCH)).catch((e) => {
    found.textContent = `recherche indisponible : ${e}`;
    return null;
  });
  let query = "";
  const flights: Record<string, unknown>[] = [];
  let flight: Record<string, unknown> | null = null;
  const fmt = (n: number) => n.toLocaleString("en-US");
  const ERRORS: Record<string, string> = {
    "malformed": "requête illisible — un gène (BRCA1), une position hg38 (chr17:43,044,295) ou un rsID (rs334)",
    "unknown-gene": "aucun gène de ce nom ni de cet identifiant dans GENCODE v50",
    "unknown-rs": "rsID absent de ClinVar et du catalogue GWAS — l'index ne couvre que ces deux sources, pas tout dbSNP",
    "unknown-chrom": "chromosome inconnu de hg38",
    "out-of-range": "position au-delà de la fin du chromosome",
  };

  /** La fiche d'une bille ou d'un monomère, ouverte comme par un clic (sans repère de pixel). */
  const showCard = async (st: CardState): Promise<void> => {
    const b = await scene.g3d.identify(st.level, st.bead);
    const sel: Selection = {
      what: st.what, level: st.level, bead: st.bead, copy: b.copy, chrom: b.copy.split(":")[0], start: b.start, end: b.end,
      pos: st.pos, strand: st.strand, precision: st.precision, detail: st.detail || undefined, world: null, pickMs: 0,
    };
    const r = await opened;
    if (!r) return;
    const card = r.card(st.level, st.bead) ?? (await r.cardAsync(st.level, st.bead));
    panel.show(compose(sel, card, r.header, r, scene.sequence?.reader ?? null), r.header, false);
    ficheOpen = true;
    cardState = st;
    dirty = true;
  };
  const precisionOf = (level: string): string => {
    const bp = scene.g3d.level(level).bp_per_bead;
    return level === "noyau" || level === scene.g3d.header.first ? `bille de ${(bp / 1e3).toFixed(0)} kb` : `monomère de ${(bp / 1e3).toFixed(0)} kb`;
  };

  const fly = async (hit: Hit, kind: string, t0: number): Promise<Destination | null> => {
    const regions: Region[] = [{ chrom: hit.chrom, start: hit.start, end: hit.end }];
    selRegions = regions;
    scene.highlight(regions);
    // Un gène : la bille de son TSS ; une position, un rsID : la leur.
    const isGene = kind === "gene" || kind === "gene-id";
    const pos = isGene ? (hit.strand === "-" ? hit.end - 1 : hit.start) : hit.start;
    const dest = await destination(scene, { chrom: hit.chrom, start: hit.start, end: hit.end, pos }, orbit);
    if (!dest) {
      found.textContent = `${hit.label} : aucune bille du modèle ne porte cette position`;
      return null;
    }
    const path = dive(orbit, dest.orbit, FOV);
    anim = { path, t0: performance.now(), duration: 900 * Math.max(0.5, path.length) };
    marker = null;
    flight = { query, kind, hit, pos, depth: dest.depth, level: dest.level, bead: dest.bead, copy: dest.copy,
               start: dest.start, end: dest.end,
               centre: dest.centre, orbit: dest.orbit, resolveMs: performance.now() - t0, durationMs: anim.duration,
               startedAt: performance.now() };
    flights.push(flight);
    found.textContent = `${hit.label} → ${hit.chrom}:${fmt(hit.start + 1)}${hit.end - hit.start > 1 ? `-${fmt(hit.end)}` : ""}` +
      `${hit.sources ? ` (${hit.sources.join(" + ")})` : ""}${hit.type ? ` · ${hit.type}` : ""}\n${depthText(dest, scene)}`;
    await showCard({ level: dest.level, bead: dest.bead, pos: isGene ? -1 : pos, strand: -1, detail: "",
                     what: dest.depth === "base" ? "base" : dest.depth === "fibre" ? "fibre" : "bille du noyau",
                     precision: dest.depth === "base" ? "la base recherchée" : precisionOf(dest.level) });
    return dest;
  };

  const search = async (text: string): Promise<Record<string, unknown>> => {
    const t0 = performance.now();
    input.value = text.trim();
    const ix = await index;
    if (!ix) return { ok: false, error: "index" };
    const res = await ix.resolve(text);
    const inModel = res.hits.filter((h) => h.in_model);
    const out: Record<string, unknown> = { resolution: res, resolveMs: performance.now() - t0 };
    if (res.error) {
      found.textContent = ERRORS[res.error] ?? res.error;
      return { ...out, flew: false };
    }
    if (!inModel.length) {
      found.textContent = `${res.hits[0].label} : ${res.hits[0].chrom} n'est pas dans le modèle (46,XX — ni chrY, ni l'ADN mitochondrial)`;
      return { ...out, flew: false };
    }
    query = text.trim();
    const dest = await fly(inModel[0], res.kind, t0);
    if (inModel.length > 1) {
      // Plusieurs loci (homonymes, rsID placé deux fois) : le premier, et les autres au choix.
      const more = document.createElement("div");
      more.textContent = `${inModel.length} cibles :`;
      inModel.forEach((h, k) => {
        const b = document.createElement("button");
        b.type = "button";
        b.textContent = `${k + 1}. ${h.chrom}:${fmt(h.start + 1)}`;
        b.onclick = () => void fly(h, res.kind, performance.now());
        more.appendChild(b);
      });
      found.appendChild(more);
    }
    return { ...out, flew: !!dest, depth: dest?.depth ?? null };
  };
  form.addEventListener("submit", (e) => {
    e.preventDefault();
    void search(input.value);
  });

  // — La vue dans l'URL —
  const stateNow = (): ViewState => ({ orbit, sel: selRegions, query, layer: layer.value, card: cardState, marker });
  let lastHash = "";
  const writeUrl = () => {
    const h = encode(stateNow());
    if (h === lastHash) return;
    lastHash = h;
    history.replaceState(null, "", `${location.pathname}${location.search}#${h}`);
  };
  const applyState = async (st: ViewState): Promise<void> => {
    anim = null;
    orbit = st.orbit;
    selRegions = st.sel;
    scene.highlight(st.sel);
    query = st.query;
    input.value = st.query;
    found.textContent = "";
    if (layer.value !== st.layer || st.layer) {
      layer.value = st.layer;
      await layer.onchange!(new Event("change"));
    }
    if (st.card) await showCard(st.card);
    else panel.hide();
    marker = st.marker;
    lastHash = encode(st);
    dirty = true;
  };
  window.addEventListener("hashchange", () => {
    const st = decode(location.hash);
    if (st && encode(st) !== lastHash) void applyState(st);
  });

  // — Signets : dans ce navigateur seulement —
  let bookmarks: Bookmark[] = loadBookmarks();
  const drawMarks = () => {
    marks.textContent = "";
    const add = document.createElement("button");
    add.type = "button";
    add.textContent = "+ signet";
    add.onclick = () => {
      writeUrl();
      bookmarks = [...bookmarks, { name: query || `vue ${bookmarks.length + 1}`, hash: lastHash, at: Date.now() }];
      if (!saveBookmarks(bookmarks)) found.textContent = "signet non enregistré : le stockage de ce navigateur est indisponible";
      drawMarks();
    };
    marks.appendChild(add);
    bookmarks.forEach((b, k) => {
      const a = document.createElement("a");
      a.href = `#${b.hash}`;
      a.textContent = b.name;
      const x = document.createElement("button");
      x.type = "button";
      x.textContent = "×";
      x.setAttribute("aria-label", `retirer le signet ${b.name}`);
      x.onclick = () => {
        bookmarks = bookmarks.filter((_, i) => i !== k);
        saveBookmarks(bookmarks);
        drawMarks();
      };
      marks.append(a, x);
    });
  };
  drawMarks();

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
    // la fenêtre n'est pas complète, l'image suivante peut montrer plus. Sinon on ne redessine
    // pas : une image en vol ferait attendre la lecture du pixel du prochain clic.
    if (scene.busy) dirty = true;
    if (dirty) {
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
    if (flight && !anim && !flight.arrivedMs) flight.arrivedMs = now - (flight.startedAt as number);
    // L'URL suit la vue au repos (replaceState : ni historique à chaque image, ni hashchange).
    if (!anim && !dirty && !drag && urlReady) writeUrl();
    requestAnimationFrame(loop);
  };
  // Une URL qui porte une vue la rejoue d'abord ; l'URL ne s'écrit qu'ensuite.
  let urlReady = false;
  const initial = decode(location.hash);
  if (initial) await applyState(initial);
  urlReady = true;
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
      /** Coût de chaque lecture, sur la vue courante : identifiant seul, profondeur seule. */
      probe(x: number, y: number, n = 20) {
        const dpr = canvas.width / canvas.clientWidth;
        const t = (f: () => unknown) => {
          const a = performance.now();
          for (let k = 0; k < n; k++) f();
          return (performance.now() - a) / n;
        };
        return { select: t(() => scene.select(x * dpr, y * dpr)), ...scene.probe(x * dpr, y * dpr, n) };
      },
      /** Plus rien à dessiner : ni image demandée, ni génération ou lecture en cours. */
      idle: () => !dirty && !anim && !scene.busy,
      /** Semaine 15 : les valeurs brutes de la surimpression, par bille, et le temps de la peindre. */
      clinical: async (level: string, kind: string) => Array.from(await clinicalRaw(level, kind), (x) => (Number.isFinite(x) ? x : null)),
      paint: async (kind: string) => { layer.value = kind; await layer.onchange!(new Event("change")); return paintMs; },
      reader: () => (reader ? { bytes: reader.bytes, requests: reader.requests, tables: reader.tablesBytes(),
        levels: Object.fromEntries(reader.header.levels.map((l) => [l.name, reader!.levelBytes(l.name)])) } : null),
    };
    publish({ ok: true, pick: true, stops: scene.stops.length, gpu: scene.gpu });
  }

  if (NAV) {
    // Pour `tools/search.mjs` : ce que la page a résolu, où le vol s'est posé, ce que l'image
    // montre à l'arrivée, et l'empreinte de l'image (WebGL et surimpression).
    const digest = async (b: Uint8Array<ArrayBuffer> | Uint8ClampedArray<ArrayBuffer>) =>
      Array.from(new Uint8Array(await crypto.subtle.digest("SHA-256", b)), (x) => x.toString(16).padStart(2, "0")).join("");
    (window as unknown as { __nav: unknown }).__nav = {
      search,
      resolve: async (text: string) => (await index)?.resolve(text) ?? null,
      /** Le texte de la fiche ouverte, ou null si elle est fermée. */
      panel: () => {
        const el = document.getElementById("fiche");
        return el?.classList.contains("open") ? el.textContent : null;
      },
      message: () => found.textContent,
      /** `n` pixels CSS au hasard (graine fixe) parmi ceux qui portent une bille du noyau et
       *  que rien ne couvre (ni le HUD, ni la fiche). */
      hits(n: number, seed: number) {
        const ids = scene.ids();
        const [W, H] = scene.size();
        const dpr = canvas.width / canvas.clientWidth;
        let st = seed >>> 0;
        const rnd = () => ((st = (st * 1664525 + 1013904223) >>> 0) / 2 ** 32);
        const out: [number, number][] = [];
        for (let tries = 0; out.length < n && tries < 400 * n; tries++) {
          const x = Math.floor(rnd() * W), y = Math.floor(rnd() * H);
          const id = ids[(H - 1 - y) * W + x];
          const cx = (x + 0.5) / dpr, cy = (y + 0.5) / dpr;
          if (id && id < 1_000_000 && document.elementFromPoint(cx, cy) === canvas) out.push([cx, cy]);
        }
        return out;
      },
      flight: () => flight,
      flights: () => flights,
      idle: () => !dirty && !anim && !scene.busy && !!(scene.arrival >= 1),
      state: () => encode(stateNow()),
      url: () => location.href,
      card: () => cardState,
      selection: () => selRegions,
      bookmarks: () => loadBookmarks(),
      /** À l'arrivée : pixels de la cible dans le tampon d'identifiants, son centre projeté, les
       *  copies surlignées, et la fiche ouverte. */
      arrived() {
        const f = flight as (Record<string, unknown> & { depth: string; bead: number; level: string; centre: number[];
                                                         pos: number }) | null;
        if (!f) return null;
        scene.render(orbit);
        const ids = scene.ids();
        const [W, H] = scene.size();
        const rank = scene.sequence ? f.pos - scene.sequence.regionStart : -1;
        let px = 0;
        for (let i = 0; i < ids.length; i++) {
          const id = ids[i];
          if (f.depth === "bille" ? id === f.bead
            : f.depth === "fibre" ? id === IDS.fiber + f.bead
            : id >= IDS.base && id < IDS.linker && Math.floor((id - IDS.base - 1) / 2) === rank) px++;
        }
        const [x, y] = scene.project(f.centre);
        const hl = scene.highlight(selRegions);
        // Ce qui est sous le centre projeté : la cible, ou ce qui la couvre.
        const cx = Math.round(x), cy = Math.round(y);
        const centreId = cx >= 0 && cx < W && cy >= 0 && cy < H ? ids[(H - 1 - cy) * W + cx] : -1;
        return { pixels: px, x, y, W, H, centreId, central: x >= W / 3 && x <= (2 * W) / 3 && y >= H / 3 && y <= (2 * H) / 3,
                 copies: hl.copies, card: cardState, depth: f.depth, widthUm: scene.last.widthUm };
      },
      async snapshot() {
        scene.render(orbit);
        const gl = scene.gl;
        const [W, H] = scene.size();
        const rgb = new Uint8Array(W * H * 4);
        gl.bindFramebuffer(gl.READ_FRAMEBUFFER, null);
        gl.readPixels(0, 0, W, H, gl.RGBA, gl.UNSIGNED_BYTE, rgb);
        drawLabels(scene);
        if (marker) drawMarker(scene, marker);
        const ov = overlay.getContext("2d")!.getImageData(0, 0, overlay.width, overlay.height).data;
        return { gl: await digest(rgb), overlay: await digest(ov), W, H, state: encode(stateNow()) };
      },
      /** Le témoin : la même vue, le lacet décalé, comparée pixel à pixel à l'image courante. */
      async perturbed(dyaw: number) {
        const gl = scene.gl;
        const [W, H] = scene.size();
        const read = () => {
          const b = new Uint8Array(W * H * 4);
          gl.bindFramebuffer(gl.READ_FRAMEBUFFER, null);
          gl.readPixels(0, 0, W, H, gl.RGBA, gl.UNSIGNED_BYTE, b);
          return b;
        };
        scene.render(orbit);
        const a = read();
        scene.render({ ...orbit, yaw: orbit.yaw + dyaw });
        const b = read();
        scene.render(orbit);
        let diff = 0;
        for (let i = 0; i < a.length; i += 4) if (a[i] !== b[i] || a[i + 1] !== b[i + 1] || a[i + 2] !== b[i + 2]) diff++;
        return { differing: diff, of: W * H };
      },
    };
    publish({ ok: true, nav: true, gpu: scene.gpu, region: scene.region });
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
