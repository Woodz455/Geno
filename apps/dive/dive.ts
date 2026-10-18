/** La scène de la plongée : deux niveaux, une cible, un rendu par orbite.
 *
 *  Partagée par la page interactive, la mesure headless et la page de mesure sur matériel
 *  réel : les trois dessinent exactement la même chose, au pixel près.
 */

import { type Frame, G3D, HttpSource, SequenceReader, type Source, type TargetSpec } from "@geno/g3d";
import {
  createTargets,
  distanceFor,
  eye,
  ImpostorRenderer,
  NO_HIT,
  type Orbit,
  perspective,
  Picker,
  Ssao,
  type Targets,
  viewMatrix,
  widthAt,
} from "@geno/viewer";
import { type Fiber, IDS, openFiber } from "./fiber.js";
import { compartmentColors, type State, stateAt, territoryPalette } from "./scene.js";
import { type Label, SequenceLayer } from "./sequence.js";

export const FOV = (45 * Math.PI) / 180;
const PREVIEW_IDS = 1_000_000;           // les identifiants de l'aperçu ne recouvrent pas ceux du noyau
const ACTB: [string, number, number] = ["chr7", 5_527_151, 5_530_601];
const AO_STRENGTH = 0.9;

/** Les cibles de la plongée : le niveau fin, la copie, la base désignée par l'annotation. */
const TARGETS = {
  actb: { level: "fin", copy: "chr7:a", target: "ACTB c.1", near: 763 },
  tp53: { level: "fin-tp53", copy: "chr17:a", target: "TP53 c.818", near: -1 },
} as const;
export type TargetName = keyof typeof TARGETS;

export interface Layer {
  name: string;
  renderer: ImpostorRenderer;
  frame: Frame;
  /** Rayon médian d'une bille, µm. */
  radiusUm: number;
  idBase: number;
  e1: Float32Array | null;
}

function layer(gl: WebGL2RenderingContext, frame: Frame, e1: Float32Array | null, idBase: number): Layer {
  const n = frame.count;
  const lv = frame.level;
  const centers = new Float32Array(n * 3);
  for (let i = 0; i < n * 3; i++) centers[i] = frame.positions[i] / 1000;   // nm → µm
  const radii = new Float32Array(n);
  const ids = new Uint32Array(n);
  const colors = new Float32Array(n * 3);
  const palette = territoryPalette(lv.copies);
  for (let i = 0; i < n; i++) {
    const c = frame.copy[i];
    radii[i] = lv.radius_nm[c] / 1000;
    ids[i] = idBase + i + 1;
    colors.set(palette[c], i * 3);
  }
  const renderer = new ImpostorRenderer(gl);
  renderer.upload({ centers, radii, ids, colors, colors2: compartmentColors(e1, n), groups: frame.copy, count: n });
  const sorted = Array.from(lv.radius_nm).sort((a, b) => a - b);
  return { name: lv.name, renderer, frame, radiusUm: sorted[sorted.length >> 1] / 1000, idBase, e1 };
}

export interface DiveScene {
  gl: WebGL2RenderingContext;
  gpu: string;
  g3d: G3D;
  preview: Layer;
  fine: Layer | null;
  /** Vue d'ensemble interactive : le noyau remplit ~40 % du champ. */
  overview: Orbit;
  /** Départ de la plongée mesurée : assez loin pour traverser le changement de niveau. */
  start: Orbit;
  /** Arrivée de la semaine 10 : la bille d'ACTB sur chr7:a, 2,6 µm de champ. */
  final: Orbit;
  /** Semaine 11 : le niveau fin dessiné, s'il est dans le fichier. */
  fiber: Fiber | null;
  /** La boucle CTCF visée, cadrée entière. */
  loopView: Orbit | null;
  /** Son ancre gauche, à l'échelle des nucléosomes (100 nm de champ). */
  anchorView: Orbit | null;
  /** La plongée entière, étape par étape (semaine 12 : jusqu'à la base). */
  stops: Orbit[];
  /** Semaine 12 : la base visée, et la couche de séquence qui la dessine. */
  target: TargetSpec | null;
  /** La région fine dessinée : son niveau, la copie du noyau où elle est posée, ses bornes. */
  region: { level: string; copy: string; chrom: string; start: number; end: number } | null;
  sequence: SequenceLayer | null;
  /** Budget de génération des nucléosomes par image, ms : un quart de l'intervalle entre les
   *  deux dernières images, entre 4 et 50 ms — un appareil lent génère moins par image, pas
   *  moins par seconde. `Infinity` pour la mesure image par image, qui le lève. */
  growBudgetMs: number;
  /** Ms entre l'ouverture et le premier rendu (aperçu seul). */
  firstRenderMs: number;
  /** État visuel du dernier rendu, avec la largeur de champ. */
  last: State & { widthUm: number };
  /** Fondu d'arrivée du noyau (0 → 1), réglé par l'appelant. */
  arrival: number;
  /** Pour la mesure seulement : impose des grandeurs d'état au lieu de les lire sur la caméra. */
  override: Partial<State> | null;
  render(o: Orbit): void;
  resize(width: number, height: number): void;
  size(): [number, number];
  /** Bille sous un pixel du canevas : texte prêt à afficher, ou null. */
  describe(x: number, y: number): Promise<string | null>;
  /** Le tampon d'identifiants du dernier rendu, ligne 0 en bas (convention WebGL). */
  ids(): Uint32Array;
  /** Pixel d'un point monde au dernier rendu, origine en haut à gauche. */
  project(p: number[]): [number, number];
  /** Lettres et étiquette de la base visée, pour la surimpression 2D du dernier rendu. */
  labels(): Label[];
  /** Génère toute la fenêtre et lit la séquence autour de la base visée — pour une mesure
   *  image par image, qui doit être déterministe et ne rend pas la main au réseau. */
  preload(): Promise<void>;
  /** Ce que l'on regarde, et d'où ça vient (principe n° 2), à l'état du dernier rendu. */
  evidence(): string;
}

export interface OpenOptions {
  /** Témoin de la mesure : transitions sèches au lieu de fondus. */
  hard?: boolean;
  prefix?: number;
  /** Témoin de lisibilité : la boucle visée n'est pas isolée du reste de la région. */
  noIsolate?: boolean;
  /** Semaine 12 : où plonge-t-on. `actb` reproduit les semaines 10 et 11. */
  target?: TargetName;
  /** Témoin de chargement : la séquence n'est lue qu'au moment où l'hélice en a besoin. */
  noPrefetch?: boolean;
}

export async function openDive(canvas: HTMLCanvasElement, file: string | Source, opts: OpenOptions = {}): Promise<DiveScene> {
  const t0 = performance.now();
  const gl = canvas.getContext("webgl2", { antialias: false, powerPreference: "high-performance" });
  if (!gl) throw new Error("WebGL2 indisponible sur cet appareil");
  const dbg = gl.getExtension("WEBGL_debug_renderer_info");
  const gpu = dbg ? String(gl.getParameter(dbg.UNMASKED_RENDERER_WEBGL)) : String(gl.getParameter(gl.RENDERER));

  let W = canvas.width;
  let H = canvas.height;
  let targets: Targets = createTargets(gl, W, H);
  let ssao = new Ssao(gl, targets);
  let picker = new Picker(gl, targets);

  const source = typeof file === "string" ? new HttpSource(file) : file;
  const g3d = await G3D.open(source, { prefix: opts.prefix ?? 16 * 1024 });
  const firstName = g3d.header.first;
  const preview = layer(gl, await g3d.firstFrame(firstName), await g3d.attribute(firstName, "compartment"), PREVIEW_IDS);
  const copies = preview.frame.level.copies;
  const groupFade = new Float32Array(copies.length);
  let focusCopy = copies.indexOf("chr7:a");

  const overview: Orbit = { target: [0, 0, 0], yaw: 0.35, pitch: 0.28, distance: distanceFor(24, FOV) };
  let lastView: Float32Array = new Float32Array(16);
  let lastProj: Float32Array = new Float32Array(16);
  let lastOrbit: Orbit = overview;
  let lastFrameAt = performance.now();
  const which = TARGETS[opts.target ?? "tp53"];

  const scene: DiveScene = {
    gl,
    gpu,
    g3d,
    preview,
    fine: null,
    overview,
    start: { ...overview, distance: distanceFor(80, FOV) },
    final: { ...overview },
    fiber: null,
    loopView: null,
    anchorView: null,
    stops: [],
    target: null,
    region: null,
    sequence: null,
    growBudgetMs: NaN,
    firstRenderMs: 0,
    last: { lod: 0, focus: 0, mix: 0, fiber: 0, loop: 0, nuc: 0, helix: 0, bases: 0, widthUm: 0 },
    arrival: 0,
    override: null,

    render(o: Orbit): void {
      const view = viewMatrix(o);
      const d = o.distance;
      // Le « champ » est la hauteur visible au point visé : la plus petite dimension d'un
      // écran paysage, celle qui dit ce qui tient dans la vue.
      const widthUm = widthAt(d, FOV);
      const fine = scene.fine;
      const fiber = scene.fiber;
      const st = stateAt({ widthUm, beadPx: fine ? (fine.radiusUm / widthUm) * H : 0 }, opts.hard);
      if (!fiber) st.fiber = st.loop = st.nuc = st.helix = st.bases = 0;
      if (!scene.sequence) st.helix = st.bases = 0;
      if (scene.override) Object.assign(st, scene.override);
      lastOrbit = o;
      // Semaine 12 : nucléosomes générés par tranches dès la boucle, séquence préchargée dès la
      // fibre, et la part d'hélice que la séquence permet vraiment de dessiner.
      const now = performance.now();
      const budget = Number.isNaN(scene.growBudgetMs) ? Math.max(4, Math.min(50, 0.25 * (now - lastFrameAt))) : scene.growBudgetMs;
      lastFrameAt = now;
      if (fiber && st.fiber > 0) fiber.grow(o.target, budget);
      let helix = 0;
      if (scene.sequence) {
        scene.sequence.prefetch(o.target, widthUm);
        helix = scene.sequence.update(o.target, widthUm, st);
      }
      // Profondeur de champ utile : tout le noyau tant qu'on le voit, la région fine ensuite.
      const reach = Math.exp((1 - st.loop) * Math.log(Math.hypot(...o.target) + 6) + st.loop * Math.log(1.2));
      const proj = perspective(FOV, W / H, Math.max(0.02 * d, d - reach), d + reach);
      lastView = view;
      lastProj = proj;
      const t = fine ? st.lod * scene.arrival : 0;
      // La copie visée cède à la fibre par tramage complémentaire : elle garde h ≥ fiber, la
      // fibre h < fiber.
      for (let c = 0; c < copies.length; c++) groupFade[c] = c === focusCopy ? 1 - st.fiber : 1 - st.focus;
      scene.last = { ...st, lod: t, widthUm };

      gl.bindFramebuffer(gl.FRAMEBUFFER, targets.fbo);
      gl.viewport(0, 0, W, H);
      gl.clearBufferfv(gl.COLOR, 0, [0.055, 0.07, 0.095, 1]);
      gl.clearBufferuiv(gl.COLOR, 1, [NO_HIT, 0, 0, 0]);
      gl.clearBufferfv(gl.COLOR, 2, [1, 0, 0, 1]);
      gl.clearBufferfi(gl.DEPTH_STENCIL, 0, 1.0, 0);
      if (t < 1) preview.renderer.draw(proj, { view, fade: 1 - t, groupFade, mix: st.mix });
      if (fine && t > 0 && st.fiber < 1) fine.renderer.draw(proj, { view, fade: t, invert: true, groupFade, mix: st.mix });
      if (fiber) fiber.draw(proj, view, H, { ...st, helix });
      if (scene.sequence) scene.sequence.draw(proj, view, H, helix);
      // Rayon d'occlusion : l'échelle de ce qu'on regarde, interpolée en logarithme — deux
      // billes de 750 kb, puis deux fois le rayon de la fibre, puis un nucléosome, puis le
      // diamètre de la double hélice.
      const r750 = 2 * (fine ?? preview).radiusUm;
      const rFib = fiber ? 2 * fiber.radiusUm : r750;
      const lr = (1 - st.nuc) * ((1 - st.fiber) * Math.log(r750) + st.fiber * Math.log(rFib)) + st.nuc * Math.log(0.012);
      const radius = Math.exp((1 - helix) * lr + helix * Math.log(0.0015));
      ssao.occlusion(proj, { radius });
      ssao.composite(proj, { radius, strength: AO_STRENGTH });
    },

    resize(width: number, height: number): void {
      W = width;
      H = height;
      ssao.dispose();
      targets = createTargets(gl, W, H);
      ssao = new Ssao(gl, targets);
      picker = new Picker(gl, targets);
    },

    size: () => [W, H],

    ids(): Uint32Array {
      const buf = new Uint32Array(W * H * 4);
      gl.bindFramebuffer(gl.READ_FRAMEBUFFER, targets.fbo);
      gl.readBuffer(gl.COLOR_ATTACHMENT1);
      gl.readPixels(0, 0, W, H, gl.RGBA_INTEGER, gl.UNSIGNED_INT, buf);
      gl.readBuffer(gl.COLOR_ATTACHMENT0);
      return buf.filter((_, i) => i % 4 === 0);
    },

    project(p: number[]): [number, number] {
      const v = lastView;
      const m = lastProj;
      const e = [0, 1, 2].map((k) => v[k] * p[0] + v[4 + k] * p[1] + v[8 + k] * p[2] + v[12 + k]);
      const cx = m[0] * e[0] + m[4] * e[1] + m[8] * e[2] + m[12];
      const cy = m[1] * e[0] + m[5] * e[1] + m[9] * e[2] + m[13];
      const cw = m[3] * e[0] + m[7] * e[1] + m[11] * e[2] + m[15];
      return [((cx / cw) * 0.5 + 0.5) * W, (1 - ((cy / cw) * 0.5 + 0.5)) * H];
    },

    async preload(): Promise<void> {
      const F = scene.fiber;
      const last = scene.stops[scene.stops.length - 1];
      if (!F || !last) return;
      F.grow(last.target, Infinity);
      if (scene.sequence) await scene.sequence.preload(last.target);
    },

    labels(): Label[] {
      return scene.sequence ? scene.sequence.labels(scene.last, (p) => scene.project(p), eye(lastOrbit), W, H) : [];
    },

    evidence(): string {
      const s = scene.last;
      if (s.helix > 0.5 && scene.sequence) {
        return "séquence : référence hg38, mesurée · hélice : B idéale, déterministe · position : nucléosomes illustratifs sur une fibre simulée (CTCF plantés)";
      }
      if (s.nuc > 0.5) return "nucléosomes : illustratifs (aucun positionnement lu) · fibre : simulée, extrusion sur CTCF plantés";
      if (s.fiber > 0.5) return "fibre de 2 kb : simulée (polymère + extrusion, CTCF plantés) · posée sur le noyau, raccord en échec (S8)";
      return "noyau : simulé (contraintes de la semaine 6) · compartiments : prédits par le modèle, pas mesurés";
    },

    async describe(x: number, y: number): Promise<string | null> {
      const id = picker.pick(x, y);
      if (id === NO_HIT) return null;
      if (id >= IDS.base && scene.sequence) return scene.sequence.describe(id);
      if (id >= IDS.fiber && scene.fiber) return scene.fiber.describe(id);
      const lay = id > PREVIEW_IDS ? preview : scene.fine!;
      const b = await g3d.identify(lay.name, id - lay.idBase);
      const e1 = lay.e1 ? lay.e1[id - lay.idBase - 1] : NaN;
      const comp = Number.isFinite(e1) ? `${e1 > 0 ? "A" : "B"} (E1 ${e1.toFixed(3)})` : "aucun";
      return (
        `${b.copy}:${b.start.toLocaleString("fr")}–${b.end.toLocaleString("fr")}  (niveau ${lay.name})\n` +
        `compartiment ${comp} · dispersion de profondeur ${b.variabilityNm.toFixed(0)} nm sur l'ensemble`
      );
    },
  };

  // Premier rendu : l'aperçu seul, dès qu'il est là.
  scene.render(overview);
  scene.firstRenderMs = performance.now() - t0;

  // Puis le noyau, ses compartiments, et la cible de la plongée.
  const fineName = g3d.header.levels.find((l) => l.name === "noyau")?.name;
  const fiberLevel = g3d.header.levels.find((l) => l.name === which.level);
  const tgt = fiberLevel?.targets?.find((t) => t.name === which.target) ?? null;
  scene.target = tgt;
  let target: [number, number, number] = [0, 0, 0];
  if (fineName) {
    scene.fine = layer(gl, await g3d.firstFrame(fineName), await g3d.attribute(fineName, "compartment"), 0);
    const locus: [string, number, number] = tgt ? [tgt.chrom, tgt.position - 1, tgt.position] : ACTB;
    const hits = await g3d.locate(fineName, ...locus);
    const hit = hits.find((h) => h.copy === which.copy) ?? hits[0];
    if (hit) {
      focusCopy = copies.indexOf(hit.copy);
      const i = hit.ids[0] - 1;
      target = [0, 1, 2].map((k) => scene.fine!.frame.positions[i * 3 + k] / 1000) as [number, number, number];
    }
  }
  scene.final = { target, yaw: overview.yaw + 1.1, pitch: 0.12, distance: distanceFor(2.6, FOV) };

  // Puis la fibre, sa boucle visée et son ancre.
  const seqSpec = tgt ? (g3d.header.sequences ?? []).find((q) => q.chrom === tgt.chrom && q.start < tgt.position && tgt.position <= q.end) : null;
  const regionStart = seqSpec?.start ?? 0;
  if (fiberLevel) {
    const chrom = which.copy.split(":")[0];
    const first = await g3d.identify(which.level, 1);
    const last = await g3d.identify(which.level, fiberLevel.n_beads);
    scene.region = { level: which.level, copy: which.copy, chrom, start: first.start, end: last.end };
  }
  const targetMonomer = tgt && fiberLevel ? Math.floor((tgt.position - 1 - regionStart) / fiberLevel.bp_per_bead) : -1;
  // ACTB rejoue les semaines 10 et 11 telles quelles : sa boucle visée est la plus proche
  // d'ACTB, pas celle qui contient une base — on ne lui passe pas de base visée.
  const actb = opts.target === "actb";
  scene.fiber = await openFiber(gl, g3d, which.level, which.near >= 0 ? which.near : targetMonomer, !opts.noIsolate, actb ? -1 : targetMonomer);
  if (scene.fiber && seqSpec && tgt) {
    scene.sequence = new SequenceLayer(gl, new SequenceReader(g3d), scene.fiber, tgt.chrom, regionStart, tgt, !opts.noPrefetch);
  }
  // Les nucléosomes autour de la base visée, dès l'ouverture : ce sont les premiers qu'une
  // plongée verra, et les générer en route laissait la descente les rattraper (S12).
  if (scene.fiber && targetMonomer >= 0) {
    const P = scene.fiber.points;
    for (const k of scene.fiber.windowNear([P[3 * targetMonomer], P[3 * targetMonomer + 1], P[3 * targetMonomer + 2]]).slice(0, 24)) {
      scene.fiber.grow([P[3 * k], P[3 * k + 1], P[3 * k + 2]], 0);
    }
  }
  const focus = scene.fiber?.focus;
  if (scene.fiber && focus) {
    // Caméra **de côté** : perpendiculaire à l'axe base → boucle et à l'écart entre les deux
    // ancres. Les deux flèches CTCF pointent alors dans le plan de l'image, vers le corps de
    // la boucle, et la convergence se lit. Deux essais l'ont précédée : l'orientation de la
    // plongée regardait la base à travers la pelote ; l'axe base → boucle mettait les deux
    // flèches dans l'axe de visée, fondues en une seule tache.
    const F = scene.fiber;
    const at = (k: number): [number, number, number] => [F.points[3 * k], F.points[3 * k + 1], F.points[3 * k + 2]];
    const a = at(focus.i);
    const b = at(focus.j);
    const base = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2];
    const ax = [focus.centroid[0] - base[0], focus.centroid[1] - base[1], focus.centroid[2] - base[2]];
    const sep = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let v = [ax[1] * sep[2] - ax[2] * sep[1], ax[2] * sep[0] - ax[0] * sep[2], ax[0] * sep[1] - ax[1] * sep[0]];
    if (Math.hypot(v[0], v[1], v[2]) < 1e-9) v = [ax[2], 0, -ax[0]];
    const l = Math.hypot(v[0], v[1], v[2]) || 1;
    const pitch = Math.max(-1.2, Math.min(1.2, Math.asin(v[1] / l)));
    const yaw = Math.atan2(v[0], v[2]);
    const mid: [number, number, number] = [
      (base[0] + focus.centroid[0]) / 2, (base[1] + focus.centroid[1]) / 2, (base[2] + focus.centroid[2]) / 2,
    ];
    const w = Math.max(0.3, Math.min(0.58, 2.4 * focus.extent));
    scene.loopView = { target: mid, yaw, pitch, distance: distanceFor(w, FOV) };
    scene.anchorView = { target: a, yaw: yaw + 0.35, pitch, distance: distanceFor(0.12, FOV) };
  }

  // Semaine 12 : la plongée jusqu'à la base. De la bille du noyau à la boucle si la base y
  // est, sinon au monomère ; puis les nucléosomes, l'hélice et la base, vue **de côté** :
  // l'œil sur la dyade de la paire, son barreau couché dans le plan de l'image.
  scene.stops = [scene.start, scene.final, scene.loopView, scene.anchorView].filter((o): o is Orbit => !!o);
  const F = scene.fiber;
  if (F && scene.sequence && !actb) {
    const rank = scene.sequence.targetRank;
    const bp = F.path.basePair(rank);
    const orient = (v: number[]): { yaw: number; pitch: number } => {
      const l = Math.hypot(v[0], v[1], v[2]) || 1;
      return { yaw: Math.atan2(v[0], v[2]), pitch: Math.asin(Math.max(-1, Math.min(1, v[1] / l))) };
    };
    let side = orient(bp.d);
    if (Math.abs(side.pitch) > 1.1) {
      // Dyade presque verticale : on tourne autour de l'axe de l'ADN jusqu'à un œil plus bas.
      const t = bp.t, d = bp.d;
      const e = [t[1] * d[2] - t[2] * d[1], t[2] * d[0] - t[0] * d[2], t[0] * d[1] - t[1] * d[0]];
      side = orient([0.6 * d[0] + 0.8 * e[0], 0.6 * d[1] + 0.8 * e[1], 0.6 * d[2] + 0.8 * e[2]]);
    }
    const P = F.points;
    const k = targetMonomer;
    const inLoop = focus && k >= focus.i && k <= focus.j;
    const fibreView: Orbit = inLoop && scene.loopView
      ? scene.loopView
      : { target: [P[3 * k], P[3 * k + 1], P[3 * k + 2]], yaw: scene.final.yaw + 0.3, pitch: scene.final.pitch, distance: distanceFor(0.6, FOV) };
    const at = (w: number, dyaw = 0): Orbit => ({ target: [bp.p[0], bp.p[1], bp.p[2]], yaw: side.yaw + dyaw, pitch: side.pitch, distance: distanceFor(w, FOV) });
    scene.stops = [scene.start, scene.final, fibreView, at(0.1, 0.5), at(0.02, 0.15), at(0.004)];
  }
  return scene;
}
