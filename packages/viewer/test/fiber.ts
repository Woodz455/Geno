/** Suite de correction de la semaine 11 : tubes et nucléosomes, dans un vrai navigateur.
 *
 *  Comme les suites des semaines 4 et 10 : chaque propriété qui compte a son **témoin**,
 *  construit pour la violer. Correction seulement ; SwiftShader ne dit rien des fps.
 */

import {
  createTargets,
  EllipsoidRenderer,
  ImpostorRenderer,
  NO_HIT,
  perspective,
  Picker,
  type Targets,
  transport,
  TubeRenderer,
} from "../src/index.js";

const W = 256;
const H = 256;
const FOV = (45 * Math.PI) / 180;
const IDENT = new Float32Array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]);

interface Result {
  name: string;
  ok: boolean;
  detail: string;
}
const results: Result[] = [];
const check = (name: string, ok: boolean, detail = "") => results.push({ name, ok, detail });

function clear(gl: WebGL2RenderingContext, t: Targets): void {
  gl.bindFramebuffer(gl.FRAMEBUFFER, t.fbo);
  gl.viewport(0, 0, W, H);
  gl.clearBufferfv(gl.COLOR, 0, [0, 0, 0, 1]);
  gl.clearBufferuiv(gl.COLOR, 1, [NO_HIT, 0, 0, 0]);
  gl.clearBufferfv(gl.COLOR, 2, [1, 0, 0, 1]);
  gl.clearBufferfi(gl.DEPTH_STENCIL, 0, 1.0, 0);
}

function readIds(gl: WebGL2RenderingContext, t: Targets): Uint32Array {
  const buf = new Uint32Array(W * H * 4);
  gl.bindFramebuffer(gl.READ_FRAMEBUFFER, t.fbo);
  gl.readBuffer(gl.COLOR_ATTACHMENT1);
  gl.readPixels(0, 0, W, H, gl.RGBA_INTEGER, gl.UNSIGNED_INT, buf);
  gl.readBuffer(gl.COLOR_ATTACHMENT0);
  return buf.filter((_, i) => i % 4 === 0);
}

/** Pixel d'un point en espace vue, origine en haut à gauche. */
function project(p: number[]): [number, number] {
  const f = 1 / Math.tan(FOV / 2);
  const w = -p[2];
  return [((f * p[0]) / w * 0.5 + 0.5) * W, (1 - ((f * p[1]) / w * 0.5 + 0.5)) * H];
}

function line(pts: number[][], radius: number, color = [0.8, 0.5, 0.3]) {
  const n = pts.length;
  return {
    points: Float32Array.from(pts.flat()),
    radii: new Float32Array(n).fill(radius),
    colors: Float32Array.from(Array.from({ length: n }, () => color).flat()),
  };
}

function transportSuite(): void {
  // Une sinusoïde plane : la courbure change de signe, le repère de Frenet s'y retourne.
  const n = 81;
  const pts: number[][] = Array.from({ length: n }, (_, k) => [k * 0.1 - 4, Math.sin(k * 0.1 * 1.5), 0]);
  const P = Float32Array.from(pts.flat());
  const N = transport(P, [0], n);
  // La propriété du transport parallèle : la normale ne tourne jamais plus que la tangente.
  // Sur une courbe plane, elle tourne exactement autant (elle doit rester perpendiculaire) ;
  // c'est le retournement de 180° là où la tangente tourne à peine qui serait la faute.
  const tangent = (k: number) => {
    const a = pts[Math.min(n - 1, k + 1)], b = pts[Math.max(0, k - 1)];
    const t = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    const l = Math.hypot(...t);
    return t.map((x) => x / l);
  };
  const angle = (u: number[], v: number[]) => Math.acos(Math.min(1, u[0] * v[0] + u[1] * v[1] + u[2] * v[2]));
  let worstPerp = 0;
  let excess = -Infinity;
  for (let k = 1; k < n; k++) {
    const t = tangent(k);
    const nk = [N[3 * k], N[3 * k + 1], N[3 * k + 2]];
    worstPerp = Math.max(worstPerp, Math.abs(nk[0] * t[0] + nk[1] * t[1] + nk[2] * t[2]));
    const turnN = angle(nk, [N[3 * k - 3], N[3 * k - 2], N[3 * k - 1]]);
    excess = Math.max(excess, turnN - angle(t, tangent(k - 1)));
  }
  check(
    "transport parallèle : normales ⊥ tangente, et jamais plus de rotation que la tangente",
    worstPerp < 1e-5 && excess < 1e-4,
    `|n·t| max ${worstPerp.toExponential(1)}, excès de rotation max ${((excess * 180) / Math.PI).toFixed(4)}°`,
  );
  // Témoin : la normale de Frenet (vers le centre de courbure) se retourne aux inflexions.
  let flips = 0;
  let prevSign = 0;
  for (let k = 1; k < n - 1; k++) {
    const curv = pts[k + 1][1] - 2 * pts[k][1] + pts[k - 1][1];
    const sign = Math.sign(curv);
    if (prevSign && sign && sign !== prevSign) flips++;
    if (sign) prevSign = sign;
  }
  check("témoin : la normale de Frenet se retourne de 180° à chaque inflexion", flips >= 2, `${flips} retournements`);
}

function tubeSuite(gl: WebGL2RenderingContext, targets: Targets, picker: Picker): void {
  const proj = perspective(FOV, W / H, 0.1, 100);

  // — Un tube droit : il passe par ses points, au bon rayon. —
  {
    const pts = Array.from({ length: 13 }, (_, k) => [-3 + k * 0.5, 0, -6]);
    const t = new TubeRenderer(gl);
    t.upload(line(pts, 0.5));
    clear(gl, targets);
    t.draw(proj, { view: IDENT, height: H });
    const ids = readIds(gl, targets);
    // Demi-largeur à l'écran, au milieu, à la profondeur de l'axe : r / 6 × f.
    const [cx, cy] = project([0, 0, -6]);
    let top = 0;
    let bottom = 0;
    for (let y = Math.round(cy); y >= 0 && ids[(H - 1 - y) * W + Math.round(cx)] !== NO_HIT; y--) top++;
    for (let y = Math.round(cy); y < H && ids[(H - 1 - y) * W + Math.round(cx)] !== NO_HIT; y++) bottom++;
    const f = 1 / Math.tan(FOV / 2);
    // Silhouette d'un cylindre vu en perspective : demi-angle asin(r/d), un peu plus que r/d.
    const expected = Math.tan(Math.asin(0.5 / 6)) * f * 0.5 * H;
    const half = (top + bottom - 1) / 2;
    clear(gl, targets);
    t.draw(proj, { view: IDENT, height: H, cullBack: true });
    const culled = readIds(gl, targets).reduce((s2, v) => s2 + (v !== NO_HIT ? 1 : 0), 0);
    const total = ids.reduce((s2, v) => s2 + (v !== NO_HIT ? 1 : 0), 0);
    check("tube : faces avant vers l'extérieur (l'élimination de l'arrière ne retire rien de la silhouette)", culled === total, `${culled} pixels sur ${total}`);
    check("tube : centré sur ses points, demi-largeur au pixel près", Math.abs(top - bottom) <= 1 && Math.abs(half - expected) <= 1.5, `${half.toFixed(1)} px, attendu ${expected.toFixed(1)}`);
    const id = picker.pick(...project([1.0, 0, -6]));
    check("tube : le picking rend le point le plus proche (point 9 à x = 1)", id === 9, `id ${id}`);
    t.setFade(Array.from({ length: 13 }, (_, k) => (k < 6 ? 0 : 1)));
    clear(gl, targets);
    t.draw(proj, { view: IDENT, height: H });
    const after = readIds(gl, targets);
    const leaked = after.reduce((s, v) => s + (v !== NO_HIT && v <= 6 ? 1 : 0), 0);
    check("tube : un intervalle éteint n'a plus aucun pixel", leaked === 0, `${leaked} pixels des points 1–6`);
    t.dispose();
  }

  // — Pas de fissure entre segments de tessellations différentes. —
  //
  // Points alternativement serrés et espacés : deux segments voisins tombent dans des
  // gabarits différents. Une grande bille d'identifiant 999 est placée derrière : une
  // fissure laisse voir 999 sur l'axe projeté du tube.
  {
    const pts: number[][] = [];
    let x = -3.2;
    for (let k = 0; k < 24; k++) {
      pts.push([x, 0.6 * Math.sin(x * 1.3), -6 + 0.4 * Math.cos(x)]);
      x += k % 2 ? 0.06 : 0.5;
    }
    const back = new ImpostorRenderer(gl);
    back.upload({ centers: Float32Array.from([0, 0, -30]), radii: Float32Array.from([20]), ids: Uint32Array.from([999]), colors: Float32Array.from([0.2, 0.2, 0.2]), count: 1 });
    const t = new TubeRenderer(gl);
    t.upload(line(pts, 0.35));
    const run = (sidesPerLevel: boolean) => {
      clear(gl, targets);
      back.draw(proj);
      t.draw(proj, { view: IDENT, height: H, sidesPerLevel, cullBack: true });
      const ids = readIds(gl, targets);
      // On sonde, autour de chaque point de jonction, un disque de 80 % du rayon projeté :
      // une fissure entre un anneau à 4 côtés et un anneau à 6 s'ouvre près de la
      // silhouette, pas sur l'axe — une première sonde de 7 × 7 pixels sur l'axe n'en
      // voyait que 2. Pas 90 % : un tube à 8 côtés n'a que cos(π/8) = 0,92 r entre deux
      // arêtes, et la sonde mordait alors sur le polygone lui-même (3 « trous » qui n'en
      // étaient pas).
      let probes = 0;
      let leaks = 0;
      const f = 1 / Math.tan(FOV / 2);
      for (const p of pts.slice(1, -1)) {
        const [px, py] = project(p);
        const rpx = Math.floor(((0.35 / -p[2]) * f * 0.5 * H) * 0.8);
        for (let dy = -rpx; dy <= rpx; dy++) {
          for (let dx = -rpx; dx <= rpx; dx++) {
            if (dx * dx + dy * dy > rpx * rpx) continue;
            const ix = Math.round(px + dx);
            const iy = Math.round(py + dy);
            if (ix < 0 || iy < 0 || ix >= W || iy >= H) continue;
            probes++;
            if (ids[(H - 1 - iy) * W + ix] === 999) leaks++;
          }
        }
      }
      return { probes, leaks, levels: t.lastLevels.filter((v) => v > 0).length };
    };
    const ok = run(false);
    check(
      "tube : aucune fissure aux jonctions de gabarits différents",
      ok.leaks === 0 && ok.levels >= 2,
      `${ok.probes} pixels sondés autour de ${pts.length - 2} jonctions, ${ok.levels} gabarits en jeu, ${ok.leaks} trous`,
    );
    const bad = run(true);
    check("témoin : un nombre de côtés par segment fend la surface", bad.leaks > 0, `${bad.leaks} pixels de fond visibles`);
    t.dispose();
    back.dispose();
  }

  // — La tessellation suit la taille à l'écran. —
  {
    const pts = Array.from({ length: 200 }, (_, k) => [Math.cos(k * 0.2) * 1.5, (k - 100) * 0.02, Math.sin(k * 0.2) * 1.5]);
    const t = new TubeRenderer(gl);
    t.upload(line(pts, 0.08));
    const at = (d: number) => {
      const view = new Float32Array(IDENT);
      view[14] = -d;
      clear(gl, targets);
      t.draw(perspective(FOV, W / H, 0.05, 500), { view, height: H });
      return { tris: t.triangles(), sides: t.lastSides, levels: [...t.lastLevels] };
    };
    const far = at(60);
    const near = at(4);
    check(
      "tube : loin, peu de triangles ; près, beaucoup plus",
      near.tris > 8 * far.tris && far.levels[0] === 199,
      `à 60 : ${far.tris} triangles (${far.sides} côtés) ; à 4 : ${near.tris} (${near.sides} côtés)`,
    );
    t.dispose();
  }
}

function ellipsoidSuite(gl: WebGL2RenderingContext, targets: Targets, picker: Picker): void {
  const proj = perspective(FOV, W / H, 0.1, 100);
  const one = (axis: number[], id = 1) => {
    const e = new EllipsoidRenderer(gl);
    e.upload({ centers: Float32Array.from([0, 0, -6]), axes: Float32Array.from(axis), radii: Float32Array.from([1.0, 0.5]), ids: Uint32Array.from([id]), colors: Float32Array.from([0.9, 0.7, 0.3]), count: 1 });
    return e;
  };
  const extent = () => {
    const ids = readIds(gl, targets);
    let minX = W, maxX = -1, minY = H, maxY = -1;
    for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) if (ids[y * W + x] !== NO_HIT) {
      minX = Math.min(minX, x); maxX = Math.max(maxX, x); minY = Math.min(minY, y); maxY = Math.max(maxY, y);
    }
    return { w: maxX - minX + 1, h: maxY - minY + 1 };
  };

  const face = one([0, 0, 1]);
  clear(gl, targets);
  face.draw(proj, { view: IDENT });
  const ef = extent();
  const dFace = picker.depthAt(W / 2, H / 2);
  check("nucléosome de face : un disque (largeur = hauteur)", Math.abs(ef.w - ef.h) <= 1, `${ef.w} × ${ef.h} px`);
  face.dispose();

  const edge = one([0, 1, 0]);
  clear(gl, targets);
  edge.draw(proj, { view: IDENT });
  const ee = extent();
  const dEdge = picker.depthAt(W / 2, H / 2);
  check("nucléosome de profil : une ellipse de rapport 1:2", Math.abs(ee.h / ee.w - 0.5) < 0.06, `${ee.w} × ${ee.h} px, rapport ${(ee.h / ee.w).toFixed(3)}`);

  // Profondeur du point d'impact : de face, le pôle à z = −6 + c ; de profil, l'équateur à −6 + a.
  const ndc = (z: number) => {
    const n = 0.1, f = 100;
    return ((-(f + n) / (f - n)) * z - (2 * f * n) / (f - n)) / -z * 0.5 + 0.5;
  };
  check(
    "nucléosome : la profondeur est celle du point d'impact (pôle de face, équateur de profil)",
    Math.abs(dFace - ndc(-5.5)) < 2e-5 && Math.abs(dEdge - ndc(-5.0)) < 2e-5 && dEdge < dFace,
    `face ${dFace.toFixed(6)} (attendu ${ndc(-5.5).toFixed(6)}), profil ${dEdge.toFixed(6)} (attendu ${ndc(-5.0).toFixed(6)})`,
  );
  edge.dispose();

  // Deux nucléosomes qui s'interpénètrent : même image dans les deux ordres de dessin.
  const pair = (order: number[]) => {
    const e = new EllipsoidRenderer(gl);
    const items = [
      { c: [-0.4, 0, -6], a: [0.3, 0.2, 0.93], id: 11 },
      { c: [0.4, 0.1, -6.1], a: [-0.5, 0.6, 0.62], id: 22 },
    ];
    const o = order.map((k) => items[k]);
    const norm = (v: number[]) => { const l = Math.hypot(...v); return v.map((x) => x / l); };
    e.upload({
      centers: Float32Array.from(o.flatMap((x) => x.c)),
      axes: Float32Array.from(o.flatMap((x) => norm(x.a))),
      radii: Float32Array.from(o.flatMap(() => [0.8, 0.4])),
      ids: Uint32Array.from(o.map((x) => x.id)),
      colors: Float32Array.from(o.flatMap(() => [0.8, 0.8, 0.8])),
      count: 2,
    });
    clear(gl, targets);
    e.draw(proj, { view: IDENT });
    const ids = readIds(gl, targets);
    e.dispose();
    return ids;
  };
  const ab = pair([0, 1]);
  const ba = pair([1, 0]);
  let differ = 0;
  let both = new Set<number>();
  for (let i = 0; i < ab.length; i++) {
    if (ab[i] !== ba[i]) differ++;
    if (ab[i] !== NO_HIT) both.add(ab[i]);
  }
  check("nucléosomes qui s'interpénètrent : l'ordre de dessin ne change rien", differ === 0 && both.size === 2, `${differ} pixels différents`);
}

function run(): void {
  transportSuite();
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
  tubeSuite(gl, targets, picker);
  ellipsoidSuite(gl, targets, picker);
}

try {
  run();
} catch (e) {
  results.push({ name: "exécution de la suite", ok: false, detail: String(e) });
}
const pre = document.createElement("pre");
pre.id = "results";
pre.textContent = JSON.stringify({ passed: results.filter((r) => r.ok).length, failed: results.filter((r) => !r.ok).length, results });
document.body.appendChild(pre);
