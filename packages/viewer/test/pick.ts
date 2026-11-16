/** Suite de correction de la semaine 13 : picking et mapping inverse, dans un vrai navigateur.
 *
 *  Comme les suites précédentes, chaque propriété a son témoin, construit pour la violer :
 *  - chaque type de rendu écrit son identifiant **et** son type (canal G de la cible
 *    RGBA32UI), y compris un identifiant de 32 bits ;
 *  - l'ADN de liaison et les cœurs de nucléosome ne partagent plus aucun identifiant —
 *    témoin : l'ancien schéma (identifiants par défaut du tube, plage des cœurs) ;
 *  - aller-retour base → point → pixel → identifiant + profondeur → point monde → base,
 *    sur 500 bases ; témoin : sans la profondeur (identifiant seul), la base n'est connue
 *    qu'au segment près.
 */

import {
  createTargets,
  DnaPath,
  distanceFor,
  EllipsoidRenderer,
  HelixRenderer,
  ImpostorRenderer,
  invert,
  KIND,
  NO_HIT,
  OCTAMER_NM,
  perspective,
  Picker,
  type Targets,
  TubeRenderer,
  viewMatrix,
} from "../src/index.js";

const W = 256;
const H = 256;
const FOV = (45 * Math.PI) / 180;
const IDENT = new Float32Array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]);
// Les plages de la plongée (apps/dive/fiber.ts), recopiées : la suite ne dépend pas de l'app.
const NUC = 4_000_000;
const LINKER = 14_000_000;
const AXIS_EVERY = 6;

interface Result {
  name: string;
  ok: boolean;
  detail: string;
}
const results: Result[] = [];
const check = (name: string, ok: boolean, detail = "") => results.push({ name, ok, detail });

type V = [number, number, number];
const len = (a: ArrayLike<number>) => Math.hypot(a[0], a[1], a[2]);
const sub = (a: ArrayLike<number>, b: ArrayLike<number>): V => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];

function mulberry32(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** Fibre synthétique, comme la suite de l'hélice : pas 46,3 nm, en µm, volume exclu. */
function fibre(n: number, seed = 5): Float32Array {
  const r = mulberry32(seed);
  const pts = new Float32Array(n * 3);
  let p: V = [0, 0, 0];
  let d: V = [1, 0, 0];
  for (let k = 0; k < n; k++) {
    pts.set(p, 3 * k);
    for (let attempt = 0; attempt < 200; attempt++) {
      const q: V = [d[0] + (r() - 0.5) * 1.6, d[1] + (r() - 0.5) * 1.6, d[2] + (r() - 0.5) * 1.6];
      const l = len(q);
      const nd: V = [q[0] / l, q[1] / l, q[2] / l];
      const np: V = [p[0] + nd[0] * 0.0463, p[1] + nd[1] * 0.0463, p[2] + nd[2] * 0.0463];
      let free = true;
      for (let j = 0; j < k; j++) free &&= len(sub(np, pts.subarray(3 * j, 3 * j + 3))) >= 0.85 * 0.0463;
      if (free || attempt === 199) {
        d = nd;
        p = np;
        break;
      }
    }
  }
  return pts;
}

function clear(gl: WebGL2RenderingContext, t: Targets): void {
  gl.bindFramebuffer(gl.FRAMEBUFFER, t.fbo);
  gl.viewport(0, 0, W, H);
  gl.clearBufferfv(gl.COLOR, 0, [0, 0, 0, 1]);
  gl.clearBufferuiv(gl.COLOR, 1, [NO_HIT, 0, 0, 0]);
  gl.clearBufferfv(gl.COLOR, 2, [1, 0, 0, 1]);
  gl.clearBufferfi(gl.DEPTH_STENCIL, 0, 1.0, 0);
}

/** Identifiant et type de chaque pixel, ligne 0 en bas. */
function readIds(gl: WebGL2RenderingContext, t: Targets): { id: Uint32Array; kind: Uint32Array } {
  const buf = new Uint32Array(W * H * 4);
  gl.bindFramebuffer(gl.READ_FRAMEBUFFER, t.fbo);
  gl.readBuffer(gl.COLOR_ATTACHMENT1);
  gl.readPixels(0, 0, W, H, gl.RGBA_INTEGER, gl.UNSIGNED_INT, buf);
  gl.readBuffer(gl.COLOR_ATTACHMENT0);
  return { id: buf.filter((_, i) => i % 4 === 0), kind: buf.filter((_, i) => i % 4 === 1) };
}

function mat4mul(a: Float32Array, b: Float32Array): Float32Array {
  const o = new Float32Array(16);
  for (let c = 0; c < 4; c++) for (let r = 0; r < 4; r++) {
    let v = 0;
    for (let k = 0; k < 4; k++) v += a[k * 4 + r] * b[c * 4 + k];
    o[c * 4 + r] = v;
  }
  return o;
}

/** Pixel d'un point monde (origine en haut à gauche), et l'inverse : pixel + profondeur → monde. */
function camera(proj: Float32Array, view: Float32Array) {
  const m = mat4mul(proj, view);
  const inv = invert(m);
  const toPixel = (p: ArrayLike<number>): [number, number] => {
    const c = [0, 1, 2, 3].map((r) => m[r] * p[0] + m[4 + r] * p[1] + m[8 + r] * p[2] + m[12 + r]);
    return [((c[0] / c[3]) * 0.5 + 0.5) * W, (1 - ((c[1] / c[3]) * 0.5 + 0.5)) * H];
  };
  const toWorld = (x: number, y: number, depth: number): V => {
    const nx = (2 * (Math.floor(x) + 0.5)) / W - 1;
    const ny = 1 - (2 * (Math.floor(y) + 0.5)) / H;
    const nz = 2 * depth - 1;
    const v = [0, 1, 2, 3].map((r) => inv[r] * nx + inv[4 + r] * ny + inv[8 + r] * nz + inv[12 + r]);
    return [v[0] / v[3], v[1] / v[3], v[2] / v[3]];
  };
  return { toPixel, toWorld };
}

/** 1. Chaque type de rendu écrit son identifiant et son type. */
function kindsSuite(gl: WebGL2RenderingContext, t: Targets, picker: Picker): void {
  const proj = perspective(FOV, 1, 0.1, 50);
  const project = (p: number[]): [number, number] => {
    const f = 1 / Math.tan(FOV / 2);
    const w = -p[2];
    return [((f * p[0]) / w * 0.5 + 0.5) * W, (1 - ((f * p[1]) / w * 0.5 + 0.5)) * H];
  };
  clear(gl, t);
  const BIG = 4_294_967_294;          // le dernier identifiant avant 0xFFFFFFFF
  const sph = new ImpostorRenderer(gl);
  sph.upload({ centers: new Float32Array([-0.9, 0.6, -4]), radii: new Float32Array([0.3]), ids: new Uint32Array([BIG]),
               colors: new Float32Array([0.8, 0.3, 0.3]), count: 1 });
  sph.draw(proj, { view: IDENT });
  const tube = new TubeRenderer(gl);
  tube.upload({ points: new Float32Array([-0.4, -0.6, -4, 0.4, -0.6, -4]), radii: new Float32Array([0.12, 0.12]),
                colors: new Float32Array([0.3, 0.8, 0.3, 0.3, 0.8, 0.3]), ids: new Uint32Array([3_999_999, 3_999_999]) });
  tube.draw(proj, { view: IDENT, height: H, idBase: LINKER });
  const disc = new EllipsoidRenderer(gl);
  disc.upload({ centers: new Float32Array([0.9, 0.6, -4]), axes: new Float32Array([0, 0, 1]), radii: new Float32Array([0.3, 0.2]),
                colors: new Float32Array([0.3, 0.3, 0.8]), count: 1, ids: new Uint32Array([NUC + 40]) });
  disc.draw(proj, { view: IDENT });
  const helix = new HelixRenderer(gl, 0.01);
  const n = 21;
  const frames = new Float32Array(n * 9);
  for (let b = 0; b < n; b++) {
    const a = (2 * Math.PI * b) / 10.5;
    frames.set([0.34 * 0.01 * (b - 10), 0.6, -4, 1, 0, 0, 0, Math.cos(a), Math.sin(a)], 9 * b);
  }
  helix.upload({ ranges: [[1000, n]], frames, bases: "ACGTACGTACGTACGTACGTA" });
  helix.draw(proj, { view: IDENT, height: H, idBase: 5_000_000 });

  const at = (p: number[]) => picker.hit(...project(p));
  const s = at([-0.9, 0.6, -4]);
  check("imposteur : identifiant 32 bits et type sphère", s?.id === BIG && s.kind === KIND.sphere, JSON.stringify(s));
  const u = at([0, -0.6, -4]);
  check("tube : plage + identifiant local, type tube", u?.id === LINKER + 3_999_999 && u.kind === KIND.tube, JSON.stringify(u));
  const d = at([0.9, 0.6, -4]);
  check("ellipsoïde : identifiant et type disque", d?.id === NUC + 40 && d.kind === KIND.disc, JSON.stringify(d));
  // L'hélice : un pixel de barreau près du centre de la fenêtre.
  const { id, kind } = readIds(gl, t);
  let hel = 0, helOk = 0;
  for (let i = 0; i < id.length; i++) if (kind[i] === KIND.helix) {
    hel++;
    const local = id[i] - 5_000_000 - 1;
    if (local >= 0 && local >> 1 >= 1000 && local >> 1 < 1000 + n) helOk++;
  }
  check("hélice : type hélice, et chaque identifiant nomme une base de sa fenêtre", hel > 0 && helOk === hel, `${helOk}/${hel} pixels`);
  const bg = picker.hit(2, 2);
  check("fond : aucun objet", bg === null, JSON.stringify(bg));
  const depth = s ? s.depth : NaN;
  check("la profondeur relue avec l'identifiant est celle du fragment", depth > 0 && depth < 1, depth.toFixed(4));
}

/** La scène de nucléosomes de la plongée, à l'un ou l'autre schéma d'identifiants. */
function nucleosomes(gl: WebGL2RenderingContext, path: DnaPath, ks: number[], scheme: "neuf" | "ancien") {
  const cores = path.cores(ks);
  const m = cores.length;
  const disc = new EllipsoidRenderer(gl);
  const centers = new Float32Array(m * 3), axes = new Float32Array(m * 3), radii = new Float32Array(m * 2), colors = new Float32Array(m * 3);
  cores.forEach((c, q) => {
    centers.set(c.center, 3 * q);
    axes.set(c.axis, 3 * q);
    radii.set([OCTAMER_NM[0] / 1000, OCTAMER_NM[1] / 1000], 2 * q);
    colors.set([0.7, 0.5, 0.4], 3 * q);
  });
  disc.upload({ centers, axes, radii, colors, count: m, ids: Uint32Array.from(cores, (c) => NUC + c.index + 1) });
  const pts: number[] = [];
  const ranks: number[] = [];
  for (const k of ks) {
    const b0 = k * path.bpPerMonomer, b1 = (k + 1) * path.bpPerMonomer;
    const ax = path.axis(b0, b1, AXIS_EVERY);
    for (let i = 0; i < ax.length; i += 3) {
      pts.push(ax[i], ax[i + 1], ax[i + 2]);
      ranks.push(Math.min(b1 - 1, b0 + (i / 3) * AXIS_EVERY));
    }
  }
  const np = pts.length / 3;
  const dna = new TubeRenderer(gl);
  dna.upload({
    points: Float32Array.from(pts), radii: new Float32Array(np).fill(0.001),
    colors: Float32Array.from(Array.from({ length: np }, () => [0.85, 0.85, 0.9]).flat()),
    ...(scheme === "neuf" ? { ids: Uint32Array.from(ranks, (b) => b + 1) } : {}),
  });
  return {
    cores,
    draw(proj: Float32Array, view: Float32Array) {
      disc.draw(proj, { view });
      dna.draw(proj, { view, height: H, idBase: scheme === "neuf" ? LINKER : NUC });
    },
  };
}

/** 2. L'ADN et les cœurs ne partagent plus d'identifiant ; 3. l'aller-retour à la base. */
function dnaSuite(gl: WebGL2RenderingContext, t: Targets, picker: Picker): void {
  const path = new DnaPath(fibre(6), 2000, 23.1);
  const ks = [2, 3];
  const c = path.monomer(2).cores[5].center;
  const orbit = { target: [c[0], c[1], c[2]] as V, yaw: 0.4, pitch: 0.3, distance: distanceFor(0.06, FOV) };
  const view = viewMatrix(orbit);
  const proj = perspective(FOV, 1, orbit.distance * 0.2, orbit.distance * 5);

  const collisions: Record<string, { both: number; tube: number; disc: number }> = {};
  for (const scheme of ["neuf", "ancien"] as const) {
    clear(gl, t);
    nucleosomes(gl, path, ks, scheme).draw(proj, view);
    const { id, kind } = readIds(gl, t);
    const discIds = new Set<number>(), tubeIds = new Set<number>();
    for (let i = 0; i < id.length; i++) {
      if (kind[i] === KIND.disc) discIds.add(id[i]);
      else if (kind[i] === KIND.tube) tubeIds.add(id[i]);
    }
    // Ce que fait le décodage de la plongée : il ne voit que l'identifiant. Un pixel d'ADN
    // dont l'identifiant tombe dans la plage des cœurs est lu comme un nucléosome.
    const coreRange = (v: number) => v > NUC && v <= NUC + path.n * 10;
    collisions[scheme] = { both: [...tubeIds].filter(coreRange).length, tube: tubeIds.size, disc: discIds.size };
  }
  check("aucun identifiant d'ADN ne tombe dans la plage des nucléosomes", collisions.neuf.both === 0 && collisions.neuf.tube > 0,
        `${collisions.neuf.both} sur ${collisions.neuf.tube} identifiants d'ADN visibles`);
  check("témoin : l'ancien schéma confondait l'ADN avec des nucléosomes", collisions.ancien.both > 0,
        `${collisions.ancien.both} identifiants d'ADN sur ${collisions.ancien.tube} lus comme des cœurs`);

  // Aller-retour : bases de lien visibles, tirées au hasard.
  clear(gl, t);
  const scene = nucleosomes(gl, path, ks, "neuf");
  scene.draw(proj, view);
  const cam = camera(proj, view);
  const r = mulberry32(13);
  const wrapped = new Set<number>();
  for (const cc of scene.cores) for (let j = 0; j < 147; j++) wrapped.add(cc.firstBp + j);
  let tried = 0, visible = 0, exact = 0, withinOne = 0, withinPixel = 0, segOnly = 0, surfaceOk = 0;
  const b0 = ks[0] * 2000, b1 = (ks[ks.length - 1] + 1) * 2000;
  while (visible < 500 && tried < 20000) {
    tried++;
    const b = b0 + Math.floor(r() * (b1 - b0));
    const [fx, fy] = cam.toPixel(path.basePoint(b));
    if (fx < 1 || fy < 1 || fx > W - 2 || fy > H - 2) continue;
    // Le pixel qui contient le point (son indice), comme un clic : le picking lit ce
    // pixel-là, et la ligne de visée passe par son centre.
    const x = Math.floor(fx), y = Math.floor(fy);
    // Résolution du pixel le long de l'axe : bases qu'il couvre. Vue de biais, l'ADN
    // raccourcit et un pixel en contient plusieurs — personne ne peut faire mieux.
    const [gx, gy] = cam.toPixel(path.basePoint(b + 1));
    const perPixel = 1 / Math.max(1e-6, Math.hypot(gx - fx, gy - fy));
    const tol = Math.max(1, Math.ceil(perPixel));
    const h = picker.hit(x, y);
    if (!h || h.kind !== KIND.tube) continue;
    const seg = h.id - LINKER - 1;
    if (Math.abs(seg - b) > 2 * AXIS_EVERY) continue;     // autre brin d'ADN devant : pas cette base
    visible++;
    // Comme la plongée : la base dont l'axe passe au plus près de la ligne de visée.
    const o = cam.toWorld(x, y, 0);
    const f = cam.toWorld(x, y, 1);
    const dl = len(sub(f, o));
    const dir = sub(f, o).map((v) => v / dl);
    const toRay = (p: ArrayLike<number>) => {
      const v = sub(p, o);
      const t = v[0] * dir[0] + v[1] * dir[1] + v[2] * dir[2];
      return len([v[0] - t * dir[0], v[1] - t * dir[1], v[2] - t * dir[2]]);
    };
    let best = seg, bd = Infinity;
    for (let q = Math.max(0, seg - AXIS_EVERY); q < seg + 2 * AXIS_EVERY; q++) {
      const d = toRay(path.basePoint(q));
      if (d < bd) [bd, best] = [d, q];
    }
    // Témoin du témoin : le point de surface, ce que la première version utilisait.
    const w = cam.toWorld(x, y, h.depth);
    let sbest = seg, sd = Infinity;
    for (let q = Math.max(0, seg - AXIS_EVERY); q < seg + 2 * AXIS_EVERY; q++) {
      const d = len(sub(path.basePoint(q), w));
      if (d < sd) [sd, sbest] = [d, q];
    }
    if (Math.abs(sbest - b) <= 1) surfaceOk++;
    if (best === b) exact++;
    if (Math.abs(best - b) <= 1) withinOne++;
    if (Math.abs(best - b) <= tol) withinPixel++;
    if (seg === b) segOnly++;
  }
  check("aller-retour base → pixel → base, par l'identifiant et la ligne de visée, à la résolution du pixel",
        visible >= 300 && withinPixel / visible >= 0.95,
        `${withinPixel}/${visible} dans la résolution du pixel (≥ 1 pb), ${withinOne} à ± 1 pb, ${exact} exactes · ${visible} bases visibles sur ${tried} tirées`);
  check("témoin : l'identifiant seul ne donne la base qu'au segment près", segOnly / visible < 0.4,
        `${segOnly}/${visible} exactes sans la ligne de visée`);
  check("témoin : le point de surface glisse le long de l'axe (r / tan θ)", surfaceOk / visible < 0.8,
        `${surfaceOk}/${visible} à ± 1 pb par le point de surface`);
  void wrapped;
}

function run(): void {
  const canvas = document.createElement("canvas");
  canvas.width = W;
  canvas.height = H;
  document.body.appendChild(canvas);
  const gl = canvas.getContext("webgl2", { antialias: false });
  if (!gl) {
    check("WebGL2 disponible", false);
    return;
  }
  const targets = createTargets(gl, W, H);
  const picker = new Picker(gl, targets);
  // Coût d'une lecture : 200 `hit` sur une image déjà rendue. SwiftShader, borne haute.
  kindsSuite(gl, targets, picker);
  const t0 = performance.now();
  for (let k = 0; k < 200; k++) picker.hit((k * 37) % W, (k * 53) % H);
  const per = (performance.now() - t0) / 200;
  check("une lecture de pixel (identifiant + profondeur) coûte moins d'une milliseconde", per < 1, `${per.toFixed(3)} ms`);
  dnaSuite(gl, targets, picker);
}

try {
  run();
} catch (e) {
  results.push({ name: "exécution de la suite", ok: false, detail: String((e as Error).stack ?? e) });
}
const pre = document.createElement("pre");
pre.id = "results";
pre.textContent = JSON.stringify({ passed: results.filter((r) => r.ok).length, failed: results.filter((r) => !r.ok).length, results });
document.body.appendChild(pre);
