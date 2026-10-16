/** Suite de correction de la semaine 10, dans un vrai navigateur.
 *
 *  Caméra et chemin de plongée, fondu tramé entre niveaux, extinction par groupe,
 *  occlusion ambiante. Chaque propriété qui compte pour une plongée « sans saut » est
 *  testée avec un **témoin** construit pour la violer : un test qui passerait aussi sur
 *  le témoin ne prouverait rien.
 *
 *  Comme la suite de la semaine 4 : correction seulement. SwiftShader ne dit rien des fps.
 */

import {
  createTargets,
  dive,
  eye,
  ImpostorRenderer,
  invert,
  lookAt,
  NO_HIT,
  type Orbit,
  perspective,
  Picker,
  RHO,
  Ssao,
  type Targets,
  type Vec3,
  viewMatrix,
  zoomPath,
} from "../src/index.js";

const W = 256;
const H = 256;
const FOV = (45 * Math.PI) / 180;

interface Result {
  name: string;
  ok: boolean;
  detail: string;
}
const results: Result[] = [];
function check(name: string, ok: boolean, detail = ""): void {
  results.push({ name, ok, detail });
}

type Item = { c: Vec3; r: number; id: number; col?: Vec3; col2?: Vec3; group?: number };

function beadsFrom(items: Item[]) {
  const n = items.length;
  const centers = new Float32Array(n * 3);
  const radii = new Float32Array(n);
  const ids = new Uint32Array(n);
  const colors = new Float32Array(n * 3);
  const colors2 = new Float32Array(n * 3);
  const groups = new Uint8Array(n);
  items.forEach((it, i) => {
    centers.set(it.c, i * 3);
    radii[i] = it.r;
    ids[i] = it.id;
    colors.set(it.col ?? [0.8, 0.8, 0.8], i * 3);
    colors2.set(it.col2 ?? [0.8, 0.8, 0.8], i * 3);
    groups[i] = it.group ?? 0;
  });
  return { centers, radii, ids, colors, colors2, groups, count: n };
}

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

function readDepth(gl: WebGL2RenderingContext, t: Targets): Float32Array {
  const buf = new Float32Array(W * H * 4);
  gl.bindFramebuffer(gl.READ_FRAMEBUFFER, t.fbo);
  gl.readBuffer(gl.COLOR_ATTACHMENT2);
  gl.readPixels(0, 0, W, H, gl.RGBA, gl.FLOAT, buf);
  gl.readBuffer(gl.COLOR_ATTACHMENT0);
  return buf.filter((_, i) => i % 4 === 0);
}

function readAo(gl: WebGL2RenderingContext, ssao: Ssao): Uint8Array {
  const buf = new Uint8Array(W * H * 4);
  gl.bindFramebuffer(gl.READ_FRAMEBUFFER, ssao.fbo);
  gl.readPixels(0, 0, W, H, gl.RGBA, gl.UNSIGNED_BYTE, buf);
  return buf.filter((_, i) => i % 4 === 0);
}

/** Point en espace vue → pixel, origine en haut à gauche. */
function project(p: Vec3): [number, number] {
  const f = 1 / Math.tan(FOV / 2);
  const w = -p[2];
  return [((f * p[0]) / w * 0.5 + 0.5) * W, (1 - ((f * p[1]) / w * 0.5 + 0.5)) * H];
}

/** Un amas reproductible de billes qui se chevauchent : le cas difficile d'un fondu. */
function cluster(n: number, idBase: number, seed = 7): Item[] {
  let s = seed;
  const rnd = () => ((s = (s * 1664525 + 1013904223) >>> 0) / 4294967296);
  return Array.from({ length: n }, (_, i) => ({
    c: [(rnd() - 0.5) * 4, (rnd() - 0.5) * 4, -8 - rnd() * 3] as Vec3,
    r: 0.35 + rnd() * 0.4,
    id: idBase + i,
    group: i % 3,
  }));
}

function cameraSuite(): void {
  // lookAt : l'œil va à l'origine, la cible sur l'axe -z à la bonne distance.
  const o: Orbit = { target: [3, -2, 5], yaw: 0.7, pitch: 0.4, distance: 12 };
  const v = viewMatrix(o);
  const apply = (p: Vec3) => [0, 1, 2].map((k) => v[k] * p[0] + v[4 + k] * p[1] + v[8 + k] * p[2] + v[12 + k]);
  const e = apply(eye(o));
  const t = apply(o.target);
  check(
    "vue : l'œil à l'origine, la cible sur -z à la distance d'orbite",
    Math.hypot(...e) < 1e-5 && Math.hypot(t[0], t[1]) < 1e-5 && Math.abs(t[2] + 12) < 1e-4,
    `œil ${Math.hypot(...e).toExponential(1)}, cible z ${t[2].toFixed(5)}`,
  );

  // Une plongée type : 20 µm de champ sur le noyau → 0,8 µm sur une région à 3 µm du centre.
  const p = zoomPath([0, 0, 0], 20_000, [1800, -1200, 2000], 800);
  const end = p.at(p.length);
  check(
    "plongée : les deux extrémités sont exactes",
    Math.hypot(end.center[0] - 1800, end.center[1] + 1200, end.center[2] - 2000) < 1e-6 && Math.abs(end.width - 800) < 1e-6,
    `écart ${Math.hypot(end.center[0] - 1800, end.center[1] + 1200, end.center[2] - 2000).toExponential(1)} nm`,
  );

  // La propriété qui compte : vitesse perçue constante, dans la métrique où un
  // déplacement compte en largeurs d'écran et un zoom en facteurs.
  const speed = (path: typeof p) => {
    const n = 400;
    const h = path.length / n;
    const v: number[] = [];
    for (let k = 1; k < n; k++) {
      const a = path.at((k - 0.5) * h);
      const b = path.at((k + 0.5) * h);
      const du = Math.hypot(b.center[0] - a.center[0], b.center[1] - a.center[1], b.center[2] - a.center[2]);
      const dw = b.width - a.width;
      const w = (a.width + b.width) / 2;
      v.push(Math.sqrt(RHO * RHO * du * du + (dw * dw) / (RHO * RHO)) / w / h);
    }
    return [Math.min(...v), Math.max(...v)];
  };
  const [lo, hi] = speed(p);
  check("plongée : vitesse perçue constante (à 0,5 % près)", hi / lo < 1.005 && Math.abs(lo - 1) < 0.01, `${lo.toFixed(4)} – ${hi.toFixed(4)}`);

  // Témoin : un zoom linéaire en distance sur le même trajet n'a pas cette propriété.
  const linear = {
    length: 1,
    at: (s: number) => ({
      center: [1800 * s, -1200 * s, 2000 * s] as Vec3,
      width: 20_000 + (800 - 20_000) * s,
    }),
  };
  const [llo, lhi] = speed(linear as typeof p);
  check("témoin : un zoom linéaire varie d'un facteur > 10 en vitesse perçue", lhi / llo > 10, `${llo.toFixed(3)} – ${lhi.toFixed(3)}`);

  const pure = zoomPath([1, 2, 3], 1000, [1, 2, 3], 10);
  const mid = pure.at(pure.length / 2);
  check("zoom pur : à mi-chemin, la moyenne géométrique des largeurs", Math.abs(mid.width - 100) < 1e-6, `${mid.width.toFixed(6)}`);

  const d = dive({ target: [0, 0, 0], yaw: 3.0, pitch: 0, distance: 10 }, { target: [1, 0, 0], yaw: -3.0, pitch: 0.3, distance: 1 }, FOV);
  const half = d.at(d.length / 2);
  check("plongée : le lacet prend le court chemin (≈ 0,28 rad, pas 6 rad)", Math.abs(Math.abs(half.yaw) - Math.PI) < 0.2, `lacet à mi-chemin ${half.yaw.toFixed(3)}`);
}

function renderSuite(gl: WebGL2RenderingContext): void {
  const targets = createTargets(gl, W, H);
  const picker = new Picker(gl, targets);
  const proj = perspective(FOV, W / H, 0.1, 100);

  // --- Repère monde : une orbite autour d'une bille la remet au centre, à sa profondeur.
  {
    const r = new ImpostorRenderer(gl);
    r.upload(beadsFrom([{ c: [4, -1, 2], r: 1, id: 9 }]));
    const ref = new ImpostorRenderer(gl);
    ref.upload(beadsFrom([{ c: [0, 0, -6], r: 1, id: 9 }]));
    clear(gl, targets);
    ref.draw(proj);
    const dRef = picker.depthAt(W / 2, H / 2);
    let worst = 0;
    let allCentred = true;
    for (const [yaw, pitch] of [[0, 0], [1.1, 0.5], [-2.5, -1.2], [3.0, 1.5]]) {
      clear(gl, targets);
      r.draw(proj, { view: viewMatrix({ target: [4, -1, 2], yaw, pitch, distance: 6 }) });
      allCentred &&= picker.pick(W / 2, H / 2) === 9;
      worst = Math.max(worst, Math.abs(picker.depthAt(W / 2, H / 2) - dRef));
    }
    check("repère monde : sous toute orbite, la bille est au centre, à la même profondeur", allCentred && worst < 1e-5, `écart ${worst.toExponential(1)}`);
    r.dispose();
    ref.dispose();
  }

  // --- Fondu tramé : deux niveaux identiques à 1 − t et t se partagent l'écran.
  {
    const A = new ImpostorRenderer(gl);
    const B = new ImpostorRenderer(gl);
    A.upload(beadsFrom(cluster(40, 1)));
    B.upload(beadsFrom(cluster(40, 1001)));
    clear(gl, targets);
    A.draw(proj);
    const idsRef = readIds(gl, targets);
    const depthRef = readDepth(gl, targets);
    const covered = idsRef.reduce((s, v) => s + (v !== NO_HIT ? 1 : 0), 0);

    const blend = (t: number, invert: boolean) => {
      clear(gl, targets);
      A.draw(proj, { fade: 1 - t });
      B.draw(proj, { fade: t, invert });
      const ids = readIds(gl, targets);
      const depth = readDepth(gl, targets);
      let holes = 0;
      let fromB = 0;
      let dz = 0;
      for (let i = 0; i < ids.length; i++) {
        if (idsRef[i] === NO_HIT) continue;
        if (ids[i] === NO_HIT) holes++;
        else {
          if (ids[i] > 1000) fromB++;
          dz = Math.max(dz, Math.abs(depth[i] - depthRef[i]));
        }
      }
      return { holes, share: fromB / covered, dz };
    };

    const runs = [0, 0.25, 0.5, 0.75, 1].map((t) => ({ t, ...blend(t, true) }));
    check(
      "fondu tramé : aucun trou à aucun t (0 ; 0,25 ; 0,5 ; 0,75 ; 1)",
      runs.every((r) => r.holes === 0),
      `${covered} pixels couverts, trous ${runs.map((r) => r.holes).join("/")}`,
    );
    check(
      "fondu tramé : la profondeur visible est celle d'un seul niveau, partout",
      runs.every((r) => r.dz < 1e-7),
      `écart max ${Math.max(...runs.map((r) => r.dz)).toExponential(1)}`,
    );
    check(
      "fondu tramé : la part du second niveau suit t à 3 points près",
      runs.every((r) => Math.abs(r.share - r.t) < 0.03),
      runs.map((r) => `${r.t}→${r.share.toFixed(3)}`).join(" "),
    );
    // Témoin : deux trames identiques au lieu de complémentaires laissent des trous.
    const control = blend(0.5, false);
    check(
      "témoin : deux trames non complémentaires trouent l'image à t = 0,5",
      control.holes > 0.3 * covered,
      `${control.holes} trous sur ${covered} (${((100 * control.holes) / covered).toFixed(0)} %)`,
    );

    // Extinction par groupe : un groupe à 0 disparaît, et le picking ne le voit plus.
    clear(gl, targets);
    A.draw(proj, { groupFade: [1, 0, 1] });
    const idsOff = readIds(gl, targets);
    const group1 = new Set(cluster(40, 1).filter((it) => it.group === 1).map((it) => it.id));
    const leaked = idsOff.reduce((s, v) => s + (group1.has(v) ? 1 : 0), 0);
    const before = idsRef.reduce((s, v) => s + (group1.has(v) ? 1 : 0), 0);
    check("extinction par groupe : aucun pixel du groupe éteint", leaked === 0 && before > 0, `${before} pixels avant, ${leaked} après`);
    A.dispose();
    B.dispose();
  }

  // --- Palette mélangée par uniforme : 0 → première, 1 → seconde.
  {
    const r = new ImpostorRenderer(gl);
    r.upload(beadsFrom([{ c: [0, 0, -6], r: 1, id: 1, col: [1, 0, 0], col2: [0, 0, 1] }]));
    clear(gl, targets);
    r.draw(proj, { mix: 0 });
    const a = picker.colorAt(W / 2, H / 2);
    clear(gl, targets);
    r.draw(proj, { mix: 1 });
    const b = picker.colorAt(W / 2, H / 2);
    check("palette : mix 0 rend la première, mix 1 la seconde", a[0] > 100 && a[2] < 5 && b[2] > 100 && b[0] < 5, `rgb(${a[0]},${a[1]},${a[2]}) → rgb(${b[0]},${b[1]},${b[2]})`);
    r.dispose();
  }

  // --- Occlusion ambiante.
  {
    const ssao = new Ssao(gl, targets);
    const scene = (k: number, items: Item[]) => items.map((it) => ({ ...it, c: it.c.map((x) => x * k) as Vec3, r: it.r * k }));
    const render = (items: Item[], projM: Float32Array, radius: number) => {
      const r = new ImpostorRenderer(gl);
      r.upload(beadsFrom(items));
      clear(gl, targets);
      r.draw(projM);
      ssao.occlusion(projM, { radius });
      r.dispose();
    };

    // Une bille seule, de face : rien ne l'occulte.
    render([{ c: [0, 0, -6], r: 1, id: 1 }], proj, 1);
    const alone = ssao.aoAt(W / 2, H / 2);
    const bg = ssao.aoAt(3, 3);
    check("occlusion : une bille seule n'est pas occultée de face", alone > 0.97, `ao ${alone.toFixed(3)}`);
    check("occlusion : le fond vaut 1", bg === 1, `ao ${bg.toFixed(3)}`);

    // Deux billes tangentes : le pli assombrit, le flanc extérieur non.
    const pair: Item[] = [
      { c: [-1, 0, -7], r: 1, id: 1 },
      { c: [1, 0, -7], r: 1, id: 2 },
    ];
    render(pair, proj, 1);
    // Trois sondes sur la bille de gauche, du flanc extérieur vers le contact. Chacune
    // est la moyenne d'une fenêtre 5×5 : un pixel seul est quantifié par pas de 1/16.
    const region = (x: number, y: number) => {
      let s = 0;
      for (let dy = -2; dy <= 2; dy++) for (let dx = -2; dx <= 2; dx++) s += ssao.aoAt(x + dx, y + dy);
      return s / 25;
    };
    const onA = (x: number): Vec3 => [x, 0, -7 + Math.sqrt(1 - (x + 1) ** 2)];
    const [cx, cy] = project(onA(-0.15));
    const probes = [-1.6, -0.6, -0.15].map((x) => region(...project(onA(x))));
    check(
      "occlusion : elle croît vers le contact entre deux billes (flanc ≥ face > pli)",
      probes[0] >= probes[1] && probes[1] > probes[2] && probes[2] < probes[0] - 1 / 32,
      probes.map((v) => v.toFixed(3)).join(" ≥ "),
    );

    // Invariance d'échelle : tout ×10, rayon ×10 → même image d'occlusion.
    const items = cluster(40, 1);
    render(items, proj, 0.6);
    const ao1 = readAo(gl, ssao);
    render(scene(10, items), perspective(FOV, W / H, 1, 1000), 6);
    const ao10 = readAo(gl, ssao);
    // Un échantillon vaut 1/16, soit 16/255 : l'arrondi flottant peut faire basculer un
    // échantillon tangent à la surface d'un côté ou de l'autre du test de profondeur.
    // On compte donc les pixels qui changent, et de combien d'échantillons.
    const compare = (a: Uint8Array, b: Uint8Array) => {
      let one = 0;
      let more = 0;
      let lit = 0;
      for (let i = 0; i < a.length; i++) {
        if (a[i] === 255 && b[i] === 255) continue;
        lit++;
        const d = Math.abs(a[i] - b[i]);
        if (d > 17) more++;
        else if (d > 1) one++;
      }
      return { one, more, lit };
    };
    const inv = compare(ao1, ao10);
    check(
      "occlusion : invariante d'échelle (scène ×10, rayon ×10) — au plus un échantillon de bascule, sur < 1 % des pixels",
      inv.more === 0 && inv.one < 0.01 * inv.lit,
      `${inv.one} pixels à un échantillon près, ${inv.more} au-delà, sur ${inv.lit} occultés`,
    );
    // Témoin : le même agrandissement sans changer le rayon change l'image.
    render(scene(10, items), perspective(FOV, W / H, 1, 1000), 0.6);
    const fixed = compare(ao1, readAo(gl, ssao));
    check(
      "témoin : à rayon fixe, l'agrandissement change l'occlusion de plus d'un échantillon",
      fixed.more > 0.2 * fixed.lit,
      `${fixed.more} pixels sur ${fixed.lit}`,
    );

    // Composition : l'occlusion assombrit la couleur, et seulement elle.
    render(pair, proj, 1);
    const fb = gl.createFramebuffer()!;
    const tex = gl.createTexture()!;
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texStorage2D(gl.TEXTURE_2D, 1, gl.RGBA8, W, H);
    gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex, 0);
    const at = (x: number, y: number) => {
      const b = new Uint8Array(4);
      gl.bindFramebuffer(gl.READ_FRAMEBUFFER, fb);
      gl.readPixels(Math.round(x), Math.round(H - 1 - y), 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, b);
      return b[0];
    };
    ssao.composite(proj, { radius: 1, strength: 0 }, fb);
    const plain = at(cx, cy);
    const plainBg = at(3, 3);
    ssao.composite(proj, { radius: 1, strength: 1 }, fb);
    const shaded = at(cx, cy);
    const shadedBg = at(3, 3);
    check("composition : le pli s'assombrit, le fond ne change pas", shaded < plain - 10 && shadedBg === plainBg, `pli ${plain} → ${shaded}, fond ${plainBg} → ${shadedBg}`);
    check("matrice inverse : P · P⁻¹ = I", (() => {
      const m = invert(proj);
      let e = 0;
      for (let i = 0; i < 4; i++)
        for (let j = 0; j < 4; j++) {
          let s = 0;
          for (let k = 0; k < 4; k++) s += proj[k * 4 + i] * m[j * 4 + k];
          e = Math.max(e, Math.abs(s - (i === j ? 1 : 0)));
        }
      return e < 1e-4;
    })());
    ssao.dispose();
  }

  void lookAt;
}

function run(): void {
  cameraSuite();
  const canvas = document.createElement("canvas");
  canvas.width = W;
  canvas.height = H;
  document.body.appendChild(canvas);
  const gl = canvas.getContext("webgl2", { antialias: false });
  if (!gl) {
    check("WebGL2 disponible", false);
    return;
  }
  renderSuite(gl);
}

try {
  run();
} catch (e) {
  results.push({ name: "exécution de la suite", ok: false, detail: String(e) });
}

const payload = {
  passed: results.filter((r) => r.ok).length,
  failed: results.filter((r) => !r.ok).length,
  results,
};
const pre = document.createElement("pre");
pre.id = "results";
pre.textContent = JSON.stringify(payload);
document.body.appendChild(pre);
