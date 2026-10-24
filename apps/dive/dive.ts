/** La scène de la plongée : deux niveaux, une cible, un rendu par orbite.
 *
 *  Partagée par la page interactive, la mesure headless et la page de mesure sur matériel
 *  réel : les trois dessinent exactement la même chose, au pixel près.
 */

import { BeadIndex, type Frame, G3D, HttpSource, SequenceReader, type Source, type TargetSpec } from "@geno/g3d";
import {
  baseOf,
  createTargets,
  distanceFor,
  eye,
  ImpostorRenderer,
  invert,
  KIND,
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

/** a × b, matrices 4 × 4 en colonnes. */
function mat4mul(a: Float32Array, b: Float32Array): Float32Array {
  const o = new Float32Array(16);
  for (let c = 0; c < 4; c++) for (let r = 0; r < 4; r++) {
    let v = 0;
    for (let k = 0; k < 4; k++) v += a[k * 4 + r] * b[c * 4 + k];
    o[c * 4 + r] = v;
  }
  return o;
}
const PREVIEW_IDS = 1_000_000;           // les identifiants de l'aperçu ne recouvrent pas ceux du noyau
const ACTB: [string, number, number] = ["chr7", 5_527_151, 5_530_601];
const AO_STRENGTH = 0.9;

/** Les cibles de la plongée : le niveau fin, la copie, la base désignée par l'annotation. */
const TARGETS = {
  actb: { level: "fin", copy: "chr7:a", target: "ACTB c.1", near: 763 },
  tp53: { level: "fin-tp53", copy: "chr17:a", target: "TP53 c.818", near: -1 },
} as const;
export type TargetName = keyof typeof TARGETS;

/** Ce qu'un clic désigne (semaine 13) : le niveau et la bille dont on lira la fiche, et
 *  la position génomique aussi fine que l'objet cliqué le permet. */
export interface Selection {
  /** L'objet : « aperçu », « bille », « fibre », « nucléosome », « ADN », « base »… */
  what: string;
  /** Niveau du `.g3d` et du fichier de fiches, et identifiant de la bille (≥ 1). */
  level: string;
  bead: number;
  copy: string;
  chrom: string;
  /** Intervalle de la bille ou du monomère, 0-based demi-ouvert. */
  start: number;
  end: number;
  /** La base, 0-based hg38, quand l'objet en désigne une ; sinon −1. */
  pos: number;
  /** Brin du barreau cliqué (0 : +, 1 : −), pour une base de l'hélice. */
  strand: 0 | 1 | -1;
  precision: string;
  detail?: string | undefined;
  /** Point monde sous le pixel, µm (profondeur relue). */
  world: number[] | null;
  /** Coût de la lecture du pixel et du mapping inverse, ms. */
  pickMs: number;
}

export interface Layer {
  name: string;
  renderer: ImpostorRenderer;
  frame: Frame;
  /** Rayon médian d'une bille, µm. */
  radiusUm: number;
  idBase: number;
  e1: Float32Array | null;
  /** Les deux palettes téléversées (territoires, compartiments), pour revenir d'une surbrillance. */
  colors: Float32Array;
  colors2: Float32Array;
  /** Semaine 15 : une couche clinique peinte à la place des deux palettes, ou null. */
  painted: Float32Array | null;
}

/** Semaine 14 : un intervalle génomique, 0-based demi-ouvert, sans copie — la sélection
 *  des vues liées est haploïde, comme le navigateur linéaire et la carte de contacts. */
export interface Region {
  chrom: string;
  start: number;
  end: number;
}

/** Ce qu'une sélection a surligné dans la 3D. */
export interface Highlighted {
  /** Identifiants locaux (≥ 1) des billes de l'aperçu et du noyau qui la recouvrent. */
  preview: number[];
  nucleus: number[];
  /** Monomères de la fibre (rangs 0-based, bornes incluses), ou null. */
  fiber: [number, number] | null;
  /** Copies qui portent une bille surlignée. */
  copies: string[];
}

/** Semaine 15 : couleur d'une valeur clinique dans [0, 1] — une seule teinte (pourpre), du clair
 *  au foncé ; gris neutre sans donnée. Le jaune de la surbrillance n'y figure pas. */
export function paintColor(v: number): [number, number, number] {
  if (!Number.isFinite(v) || v <= 0) return [0.36, 0.39, 0.44];
  const t = Math.min(1, v);
  return [0.93 - 0.45 * t, 0.78 - 0.68 * t, 0.95 - 0.35 * t];
}

/** Okabe–Ito jaune : absent des palettes de territoire et de compartiment. */
export const HIGHLIGHT: [number, number, number] = [0.94, 0.89, 0.26];
const HL_DIM = 0.4;          // les billes non sélectionnées gardent 40 % de leur couleur
const HL_FADE = 0.15;        // les copies sans sélection, 15 % de leurs pixels
const HL_OWN = 0.3;          // le reste d'une copie qui porte la sélection, 30 % : la bille
                             // surlignée ne reste pas enfouie dans son propre territoire
/** Groupes des billes surlignées : celles de la copie visée (qui cèdent à la fibre), les autres. */
const SEL_FOCUS = 62;
const SEL_OTHER = 63;

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
  const colors2 = compartmentColors(e1, n);
  return { name: lv.name, renderer, frame, radiusUm: sorted[sorted.length >> 1] / 1000, idBase, e1, colors, colors2, painted: null };
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
  /** Le dernier rendu a généré des nucléosomes, attend de la séquence ou en lit : l'image
   *  suivante peut montrer plus. Sinon, inutile de redessiner (semaine 13 : une image en
   *  vol fait attendre la lecture du pixel d'un clic). */
  busy: boolean;
  render(o: Orbit): void;
  resize(width: number, height: number): void;
  size(): [number, number];
  /** Bille sous un pixel du canevas : texte prêt à afficher, ou null. */
  describe(x: number, y: number): Promise<string | null>;
  /** Semaine 13 : ce que désigne un pixel, **synchrone** — une lecture de pixel et des
   *  tableaux déjà en mémoire. `null` sur le fond, ou si le niveau n'est pas encore lu. */
  select(x: number, y: number): Selection | null;
  /** Pour la mesure : coût moyen d'une lecture d'identifiant et d'une lecture de profondeur. */
  probe(x: number, y: number, n: number): { id: number; depth: number };
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
  /** Semaine 14 : surligne ce qui recouvre `regions` — billes des deux niveaux en jaune, les
   *  autres assombries, les copies sans sélection presque effacées, monomères de la fibre ;
   *  `null` rend la scène telle qu'elle était. Synchrone : les intervalles sont en mémoire. */
  highlight(regions: Region[] | null): Highlighted;
  /** Semaine 15 : peint une couche par bille (`values[i]` dans [0, 1], NaN : pas de donnée)
   *  sur le niveau nommé, à la place des territoires et des compartiments ; `null` les rend. */
  paint(level: string, values: Float32Array | null): void;
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

/** Une paire de bases vue **de côté** (semaine 12) : l'œil sur la dyade de la paire `rank`
 *  (rang dans la région fine), son barreau couché dans le plan de l'image, `width` µm de champ.
 *  Partagé par la plongée et, en semaine 16, par la recherche d'une position. */
export function sideOrbit(F: Fiber, rank: number, width: number, dyaw = 0): Orbit {
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
  return { target: [bp.p[0], bp.p[1], bp.p[2]], yaw: side.yaw + dyaw, pitch: side.pitch, distance: distanceFor(width, FOV) };
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
  const beads = new BeadIndex(g3d);
  const firstName = g3d.header.first;
  const preview = layer(gl, await g3d.firstFrame(firstName), await g3d.attribute(firstName, "compartment"), PREVIEW_IDS);
  const copies = preview.frame.level.copies;
  const groupFade = new Float32Array(64);       // MAX_GROUPS : les copies, puis SEL_FOCUS et SEL_OTHER
  let focusCopy = copies.indexOf("chr7:a");

  const overview: Orbit = { target: [0, 0, 0], yaw: 0.35, pitch: 0.28, distance: distanceFor(24, FOV) };
  let lastView: Float32Array = new Float32Array(16);
  let lastProj: Float32Array = new Float32Array(16);
  let lastOrbit: Orbit = overview;
  let lastFrameAt = performance.now();
  const which = TARGETS[opts.target ?? "tp53"];
  /** Semaine 14 : les copies qui portent une bille surlignée ; `null` sans sélection. */
  let hlCopy: Uint8Array | null = null;
  let lastRegions: Region[] | null = null;

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
    busy: false,

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
      const grown = fiber && st.fiber > 0 ? fiber.grow(o.target, budget) : 0;
      let helix = 0;
      let waiting = false;
      if (scene.sequence) {
        scene.sequence.prefetch(o.target, widthUm);
        const before = scene.sequence.stats.loading;
        helix = scene.sequence.update(o.target, widthUm, st);
        waiting = scene.sequence.stats.loading > before || scene.sequence.pending;
      }
      scene.busy = grown > 0 || waiting;
      // Profondeur de champ utile : tout le noyau tant qu'on le voit, la région fine ensuite.
      const reach = Math.exp((1 - st.loop) * Math.log(Math.hypot(...o.target) + 6) + st.loop * Math.log(1.2));
      const proj = perspective(FOV, W / H, Math.max(0.02 * d, d - reach), d + reach);
      lastView = view;
      lastProj = proj;
      const t = fine ? st.lod * scene.arrival : 0;
      // La copie visée cède à la fibre par tramage complémentaire : elle garde h ≥ fiber, la
      // fibre h < fiber.
      // Une sélection (semaine 14) l'emporte sur l'effacement de la plongée : une copie qui
      // porte une bille surlignée reste entière à toutes les échelles, les autres tombent à
      // HL_FADE. La copie visée cède toujours à sa fibre.
      for (let c = 0; c < copies.length; c++) {
        const base = c === focusCopy ? 1 - st.fiber : 1 - st.focus;
        groupFade[c] = !hlCopy ? base : (hlCopy[c] ? (c === focusCopy ? 1 - st.fiber : 1) : base) * (hlCopy[c] ? HL_OWN : HL_FADE);
      }
      groupFade[SEL_FOCUS] = 1 - st.fiber;
      groupFade[SEL_OTHER] = 1;
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

    highlight(regions: Region[] | null): Highlighted {
      lastRegions = regions;
      const out: Highlighted = { preview: [], nucleus: [], fiber: null, copies: [] };
      const on = !!regions && regions.length > 0;
      hlCopy = on ? new Uint8Array(copies.length) : null;
      const hit = (chrom: string, s: number, e: number) => !!regions && regions.some((r) => r.chrom === chrom && r.start < e && r.end > s);
      for (const [lay, ids] of [[preview, out.preview], [scene.fine, out.nucleus]] as const) {
        if (!lay) continue;
        const b = beads.now(lay.name);
        const n = lay.frame.count;
        const c1 = (lay.painted ?? lay.colors).slice();
        const c2 = (lay.painted ?? lay.colors2).slice();
        const groups = lay.frame.copy.slice();
        if (on && b) {
          const labels = lay.frame.level.copies;
          for (let i = 0; i < n; i++) {
            const copy = lay.frame.copy[i];
            if (hit(labels[copy].split(":")[0], b.start[i], b.end[i])) {
              ids.push(i + 1);
              c1.set(HIGHLIGHT, 3 * i);
              c2.set(HIGHLIGHT, 3 * i);
              groups[i] = copy === focusCopy ? SEL_FOCUS : SEL_OTHER;
              if (lay === scene.fine || !scene.fine) hlCopy![copy] = 1;
            } else {
              for (let k = 0; k < 3; k++) {
                c1[3 * i + k] *= HL_DIM;
                c2[3 * i + k] *= HL_DIM;
              }
            }
          }
        }
        lay.renderer.recolor(c1, c2, groups);
      }
      out.copies = copies.filter((_, c) => hlCopy?.[c]);
      const R = scene.region;
      const F = scene.fiber;
      if (R && F) {
        const b = beads.now(R.level);
        let lo = -1, hi = -1;
        if (on && b) {
          for (let k = 0; k < b.start.length; k++) {
            if (hit(R.chrom, b.start[k], b.end[k])) {
              if (lo < 0) lo = k;
              hi = k;
            }
          }
        }
        F.highlight(lo, hi);
        out.fiber = lo >= 0 ? [lo, hi] : null;
      }
      return out;
    },

    paint(level: string, values: Float32Array | null): void {
      const lay = [preview, scene.fine].find((l) => l?.name === level);
      if (!lay) return;
      if (!values) lay.painted = null;
      else {
        const c = new Float32Array(lay.frame.count * 3);
        for (let i = 0; i < lay.frame.count; i++) c.set(paintColor(values[i]), 3 * i);
        lay.painted = c;
      }
      // Même chemin que la surbrillance, qui part de la peinture si elle existe.
      scene.highlight(lastRegions);
    },

    probe(x: number, y: number, n: number) {
      const t = (f: () => unknown) => {
        const a = performance.now();
        for (let k = 0; k < n; k++) f();
        return (performance.now() - a) / n;
      };
      return { id: t(() => picker.pick(x, y)), depth: t(() => picker.depthAt(x, y)) };
    },

    select(x: number, y: number): Selection | null {
      const t0 = performance.now();
      const hit = picker.hit(x, y);
      if (!hit) return null;
      const { id, kind } = hit;
      // Le point monde (profondeur relue) et la ligne de visée du pixel : l'inverse de
      // projection × vue, au centre du pixel que le picking a lu.
      const m = invert(mat4mul(lastProj, lastView));
      const nx = (2 * (Math.floor(x) + 0.5)) / W - 1;
      const ny = 1 - (2 * (Math.floor(y) + 0.5)) / H;
      const unproject = (nz: number) => {
        const v = [0, 1, 2, 3].map((r) => m[r] * nx + m[4 + r] * ny + m[8 + r] * nz + m[12 + r]);
        return [v[0] / v[3], v[1] / v[3], v[2] / v[3]];
      };
      const world = unproject(2 * hit.depth - 1);
      const near = unproject(-1);
      const far = unproject(1);
      const dl = Math.hypot(far[0] - near[0], far[1] - near[1], far[2] - near[2]) || 1;
      const ray = { o: near, d: [(far[0] - near[0]) / dl, (far[1] - near[1]) / dl, (far[2] - near[2]) / dl] };
      const done = (s: Omit<Selection, "pickMs" | "world">): Selection => ({ ...s, world, pickMs: performance.now() - t0 });
      const bead = (lay: Layer, local: number, what: string, precision: string): Selection | null => {
        const b = beads.now(lay.name);
        if (!b) return null;
        const copy = lay.frame.level.copies[lay.frame.copy[local - 1]];
        return done({ what, level: lay.name, bead: local, copy, chrom: copy.split(":")[0], start: b.start[local - 1],
                      end: b.end[local - 1], pos: -1, strand: -1, precision });
      };
      const F = scene.fiber;
      const R = scene.region;
      const fibre = (monomer: number, pos: number, strand: 0 | 1 | -1, what: string, precision: string, detail?: string) => {
        const b = beads.now(R!.level);
        if (!b || monomer < 0 || monomer >= b.start.length) return null;
        return done({ what, level: R!.level, bead: monomer + 1, copy: R!.copy, chrom: R!.chrom, start: b.start[monomer],
                      end: b.end[monomer], pos, strand, precision, detail });
      };
      if (kind === KIND.helix && id >= IDS.base && R && F) {
        const { rank, strand } = baseOf(id - IDS.base);
        return fibre(Math.floor(rank / g3d.level(R.level).bp_per_bead), R.start + rank, strand, "base", "exacte : un barreau = une base");
      }
      if (id >= IDS.fiber && F && R) {
        const h = F.locate(id, ray);
        if (!h) return null;
        const what = { fibre: "fibre", ctcf: "site CTCF", anneau: "cohésine", arc: "paire d'ancres", nucleosome: "nucléosome", adn: "ADN" }[h.what];
        return fibre(h.monomer, h.rank >= 0 ? R.start + h.rank : -1, -1, what, h.precision, h.detail);
      }
      if (id > PREVIEW_IDS) return bead(preview, id - PREVIEW_IDS, "bille de l'aperçu", `bille de ${(preview.frame.level.bp_per_bead / 1e6).toFixed(0)} Mb`);
      if (scene.fine) return bead(scene.fine, id, "bille du noyau", `bille de ${(scene.fine.frame.level.bp_per_bead / 1e3).toFixed(0)} kb`);
      return null;
    },

    async describe(x: number, y: number): Promise<string | null> {
      const id = picker.pick(x, y);
      if (id === NO_HIT) return null;
      if (id >= IDS.linker && scene.fiber) return scene.fiber.describe(id);
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
    // Semaine 13 : début et fin de chaque bille, une fois, pour qu'un clic n'attende rien.
    await Promise.all([beads.load(firstName), beads.load(fineName)]);
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
    await beads.load(which.level);
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
    const P = F.points;
    const k = targetMonomer;
    const inLoop = focus && k >= focus.i && k <= focus.j;
    const fibreView: Orbit = inLoop && scene.loopView
      ? scene.loopView
      : { target: [P[3 * k], P[3 * k + 1], P[3 * k + 2]], yaw: scene.final.yaw + 0.3, pitch: scene.final.pitch, distance: distanceFor(0.6, FOV) };
    const at = (w: number, dyaw = 0): Orbit => sideOrbit(F, rank, w, dyaw);
    scene.stops = [scene.start, scene.final, fibreView, at(0.1, 0.5), at(0.02, 0.15), at(0.004)];
  }
  return scene;
}
