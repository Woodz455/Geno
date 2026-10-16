/** La scène de la plongée : deux niveaux, une cible, un rendu par orbite.
 *
 *  Partagée par la page interactive, la mesure headless et la page de mesure sur matériel
 *  réel : les trois dessinent exactement la même chose, au pixel près.
 */

import { type Frame, G3D, HttpSource } from "@geno/g3d";
import {
  createTargets,
  distanceFor,
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
import { compartmentColors, type State, stateAt, territoryPalette } from "./scene.js";

export const FOV = (45 * Math.PI) / 180;
const PREVIEW_IDS = 1_000_000;           // les identifiants de l'aperçu ne recouvrent pas ceux du noyau
const ACTB: [string, number, number] = ["chr7", 5_527_151, 5_530_601];
const AO_STRENGTH = 0.9;

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
  /** Arrivée : la bille d'ACTB sur chr7:a, 2,6 µm de champ. */
  final: Orbit;
  /** Ms entre l'ouverture et le premier rendu (aperçu seul). */
  firstRenderMs: number;
  /** État visuel du dernier rendu, avec la largeur de champ. */
  last: State & { widthUm: number };
  /** Fondu d'arrivée du noyau (0 → 1), réglé par l'appelant. */
  arrival: number;
  render(o: Orbit): void;
  resize(width: number, height: number): void;
  size(): [number, number];
  /** Bille sous un pixel du canevas : texte prêt à afficher, ou null. */
  describe(x: number, y: number): Promise<string | null>;
}

export interface OpenOptions {
  /** Témoin de la mesure : transitions sèches au lieu de fondus. */
  hard?: boolean;
  prefix?: number;
}

export async function openDive(canvas: HTMLCanvasElement, file: string, opts: OpenOptions = {}): Promise<DiveScene> {
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

  const g3d = await G3D.open(new HttpSource(file), { prefix: opts.prefix ?? 16 * 1024 });
  const firstName = g3d.header.first;
  const preview = layer(gl, await g3d.firstFrame(firstName), await g3d.attribute(firstName, "compartment"), PREVIEW_IDS);
  const copies = preview.frame.level.copies;
  const groupFade = new Float32Array(copies.length);
  let focusCopy = copies.indexOf("chr7:a");

  const overview: Orbit = { target: [0, 0, 0], yaw: 0.35, pitch: 0.28, distance: distanceFor(24, FOV) };

  const scene: DiveScene = {
    gl,
    gpu,
    g3d,
    preview,
    fine: null,
    overview,
    start: { ...overview, distance: distanceFor(80, FOV) },
    final: { ...overview },
    firstRenderMs: 0,
    last: { lod: 0, focus: 0, mix: 0, widthUm: 0 },
    arrival: 0,

    render(o: Orbit): void {
      const view = viewMatrix(o);
      const d = o.distance;
      const reach = Math.hypot(...o.target) + 6;          // rayon du noyau 5 µm, et de la marge
      const proj = perspective(FOV, W / H, Math.max(0.02 * d, d - reach), d + reach);
      // Le « champ » est la hauteur visible au point visé : la plus petite dimension d'un
      // écran paysage, celle qui dit ce qui tient dans la vue.
      const widthUm = widthAt(d, FOV);
      const fine = scene.fine;
      const st = stateAt({ widthUm, beadPx: fine ? (fine.radiusUm / widthUm) * H : 0 }, opts.hard);
      const t = fine ? st.lod * scene.arrival : 0;
      for (let c = 0; c < copies.length; c++) groupFade[c] = c === focusCopy ? 1 : 1 - st.focus;
      scene.last = { ...st, lod: t, widthUm };

      gl.bindFramebuffer(gl.FRAMEBUFFER, targets.fbo);
      gl.viewport(0, 0, W, H);
      gl.clearBufferfv(gl.COLOR, 0, [0.055, 0.07, 0.095, 1]);
      gl.clearBufferuiv(gl.COLOR, 1, [NO_HIT, 0, 0, 0]);
      gl.clearBufferfv(gl.COLOR, 2, [1, 0, 0, 1]);
      gl.clearBufferfi(gl.DEPTH_STENCIL, 0, 1.0, 0);
      if (t < 1) preview.renderer.draw(proj, { view, fade: 1 - t, groupFade, mix: st.mix });
      if (fine && t > 0) fine.renderer.draw(proj, { view, fade: t, invert: true, groupFade, mix: st.mix });
      const radius = 2 * (fine ?? preview).radiusUm;
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

    async describe(x: number, y: number): Promise<string | null> {
      const id = picker.pick(x, y);
      if (id === NO_HIT) return null;
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
  let target: [number, number, number] = [0, 0, 0];
  if (fineName) {
    scene.fine = layer(gl, await g3d.firstFrame(fineName), await g3d.attribute(fineName, "compartment"), 0);
    const hits = await g3d.locate(fineName, ...ACTB);
    const hit = hits.find((h) => h.copy === "chr7:a") ?? hits[0];
    if (hit) {
      focusCopy = copies.indexOf(hit.copy);
      const i = hit.ids[0] - 1;
      target = [0, 1, 2].map((k) => scene.fine!.frame.positions[i * 3 + k] / 1000) as [number, number, number];
    }
  }
  scene.final = { target, yaw: overview.yaw + 1.1, pitch: 0.12, distance: distanceFor(2.6, FOV) };
  return scene;
}
