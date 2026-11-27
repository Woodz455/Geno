/** Suite de correction des imposteurs, exécutée dans un vrai navigateur.
 *
 *  Ce qu'elle prouve : la technique est correcte. Ce qu'elle ne prouve pas : la
 *  performance. En headless, Chromium rend via SwiftShader, un rasteriseur
 *  **logiciel** — ses fps ne disent rien d'un vrai GPU. Les mesures de
 *  performance sont dans `apps/spike`, à exécuter sur du matériel réel.
 */

import { createTargets, ImpostorRenderer, NO_HIT, perspective, Picker } from "../src/index.js";

const W = 256;
const H = 256;
const FOV = (45 * Math.PI) / 180;
const NEAR = 0.1;
const FAR = 100;

interface Result {
  name: string;
  ok: boolean;
  detail: string;
}

const results: Result[] = [];

function check(name: string, ok: boolean, detail = ""): void {
  results.push({ name, ok, detail });
}

/** Point en espace vue → pixel écran, même convention que le picking. */
function project(p: [number, number, number]): [number, number] {
  const f = 1 / Math.tan(FOV / 2);
  const w = -p[2];
  const ndcX = (f * p[0]) / w;
  const ndcY = (f * p[1]) / w;
  return [(ndcX * 0.5 + 0.5) * W, (1 - (ndcY * 0.5 + 0.5)) * H];
}

function beadsFrom(
  items: { c: [number, number, number]; r: number; id: number; col: [number, number, number] }[],
) {
  const n = items.length;
  const centers = new Float32Array(n * 3);
  const radii = new Float32Array(n);
  const ids = new Uint32Array(n);
  const colors = new Float32Array(n * 3);
  items.forEach((it, i) => {
    centers.set(it.c, i * 3);
    radii[i] = it.r;
    ids[i] = it.id;
    colors.set(it.col, i * 3);
  });
  return { centers, radii, ids, colors, count: n };
}

function run(): void {
  const canvas = document.createElement("canvas");
  canvas.width = W;
  canvas.height = H;
  document.body.appendChild(canvas);

  const gl = canvas.getContext("webgl2", { antialias: false });
  if (!gl) {
    check("WebGL2 disponible", false, "getContext a rendu null");
    return;
  }
  check("WebGL2 disponible", true, gl.getParameter(gl.VERSION) as string);

  let targets;
  try {
    targets = createTargets(gl, W, H);
    check("cible multiple couleur + identifiants + profondeur", true, "framebuffer complet");
  } catch (e) {
    check("cible multiple couleur + identifiants + profondeur", false, String(e));
    return;
  }

  let renderer: ImpostorRenderer;
  try {
    renderer = new ImpostorRenderer(gl);
    check("shaders compilés et liés", true, "GLSL ES 3.0");
  } catch (e) {
    check("shaders compilés et liés", false, String(e));
    return;
  }

  const proj = perspective(FOV, W / H, NEAR, FAR);
  const picker = new Picker(gl, targets);

  const drawScene = (beads: ReturnType<typeof beadsFrom>) => {
    const r = new ImpostorRenderer(gl);
    r.upload(beads);
    gl.bindFramebuffer(gl.FRAMEBUFFER, targets.fbo);
    gl.viewport(0, 0, W, H);
    gl.clearBufferfv(gl.COLOR, 0, [0, 0, 0, 1]);
    gl.clearBufferuiv(gl.COLOR, 1, [NO_HIT, 0, 0, 0]);
    gl.clearBufferfv(gl.COLOR, 2, [1, 0, 0, 1]);
    gl.clearBufferfi(gl.DEPTH_STENCIL, 0, 1.0, 0);
    r.draw(proj);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    return r;
  };

  // --- 1. Une bille se rend comme un disque, pas comme un carré -------------
  {
    const beads = beadsFrom([{ c: [0, 0, -6], r: 1, id: 7, col: [1, 0.2, 0.2] }]);
    const r = drawScene(beads);
    const [cx, cy] = project([0, 0, -6]);
    const centre = picker.pick(cx, cy);
    // Le coin du quad est à ~51 px du centre sur chaque axe ; en diagonale il
    // est donc hors du disque. Un billboard non découpé le remplirait.
    const coin = picker.pick(cx + 44, cy + 44);
    const dehors = picker.pick(cx + 90, cy);
    check("le centre de la bille est touché", centre === 7, `id=${centre}`);
    check("le coin du quad est jeté (disque, pas carré)", coin === NO_HIT, `id=${coin}`);
    check("hors du disque, aucun identifiant", dehors === NO_HIT, `id=${dehors}`);
    r.dispose();
  }

  // --- 2. La profondeur bombe : c'est tout l'intérêt des imposteurs ---------
  {
    const beads = beadsFrom([{ c: [0, 0, -6], r: 1, id: 1, col: [1, 1, 1] }]);
    const r = drawScene(beads);
    const [cx, cy] = project([0, 0, -6]);
    const dCentre = picker.depthAt(cx, cy);
    const dBord = picker.depthAt(cx + 45, cy);
    // Au centre on touche le pôle (z = -5) ; près du bord, le limbe (z ≈ -6).
    check(
      "la profondeur bombe au centre",
      dCentre < dBord - 1e-3,
      `centre ${dCentre.toFixed(5)} < bord ${dBord.toFixed(5)}`,
    );
    r.dispose();
  }

  // --- 3. L'interpénétration est résolue par fragment, pas par quad ---------
  //
  // Deux billes de même profondeur qui se chevauchent. Si la profondeur était
  // celle du quad, les deux quads seraient coplanaires et l'ordre de dessin
  // déciderait du vainqueur. Avec la profondeur du point d'impact, c'est la
  // surface la plus proche qui gagne — donc le même résultat dans les deux
  // ordres de dessin.
  {
    const A = { c: [-0.5, 0, -6] as [number, number, number], r: 1, id: 11, col: [1, 0, 0] as [number, number, number] };
    const B = { c: [0.5, 0, -6] as [number, number, number], r: 1, id: 22, col: [0, 0, 1] as [number, number, number] };
    const [gxA, gy] = project([-0.35, 0, -6]);
    const [gxB] = project([0.35, 0, -6]);

    const r1 = drawScene(beadsFrom([A, B]));
    const gauche1 = picker.pick(gxA, gy);
    const droite1 = picker.pick(gxB, gy);
    r1.dispose();

    const r2 = drawScene(beadsFrom([B, A]));
    const gauche2 = picker.pick(gxA, gy);
    const droite2 = picker.pick(gxB, gy);
    r2.dispose();

    check(
      "la surface la plus proche gagne à gauche du plan d'intersection",
      gauche1 === 11 && gauche2 === 11,
      `ordre AB ${gauche1}, ordre BA ${gauche2}`,
    );
    check(
      "la surface la plus proche gagne à droite",
      droite1 === 22 && droite2 === 22,
      `ordre AB ${droite1}, ordre BA ${droite2}`,
    );
    check(
      "l'ordre de dessin ne change pas le résultat",
      gauche1 === gauche2 && droite1 === droite2,
      "",
    );
  }

  // --- 4. Le picking rend le bon identifiant, pour chaque bille -------------
  {
    const items = [
      { c: [-2, 1.2, -8] as [number, number, number], r: 0.5, id: 101, col: [1, 0, 0] as [number, number, number] },
      { c: [0, 1.2, -8] as [number, number, number], r: 0.5, id: 202, col: [0, 1, 0] as [number, number, number] },
      { c: [2, 1.2, -8] as [number, number, number], r: 0.5, id: 303, col: [0, 0, 1] as [number, number, number] },
      { c: [-1, -1.2, -8] as [number, number, number], r: 0.5, id: 404, col: [1, 1, 0] as [number, number, number] },
      { c: [1, -1.2, -8] as [number, number, number], r: 0.5, id: 505, col: [0, 1, 1] as [number, number, number] },
    ];
    const r = drawScene(beadsFrom(items));
    const wrong = items.filter((it) => {
      const [x, y] = project(it.c);
      return picker.pick(x, y) !== it.id;
    });
    check("chaque bille rend son propre identifiant", wrong.length === 0, `${items.length - wrong.length}/${items.length}`);
    check("le fond rend NO_HIT", picker.pick(4, 4) === NO_HIT, "");
    check("hors cadre rend NO_HIT", picker.pick(-10, -10) === NO_HIT && picker.pick(9999, 9999) === NO_HIT, "");
    r.dispose();
  }

  // --- 5. À l'échelle, et sur toute la plage d'identifiants ----------------
  //
  // On ne sonde pas les billes du nuage de fond : à 250 000 billes sur 256 px,
  // chacune couvre moins d'un pixel et n'est physiquement pas pointable. Ce
  // n'est pas un défaut du picking, c'est la résolution de l'écran. On sonde
  // donc des billes franches placées devant, et le nuage sert à ce qu'il doit
  // servir — prouver que le coût et l'exactitude ne dépendent pas du nombre.
  {
    const BACKGROUND = 250_000;
    // Plage complète d'uint32 : le dernier tient dans 32 bits et pas dans 24.
    const probes = [
      { c: [-2.2, 1.3, -8] as [number, number, number], r: 0.35, id: 1 },
      { c: [0, 1.3, -8] as [number, number, number], r: 0.35, id: 65_536 },
      { c: [2.2, 1.3, -8] as [number, number, number], r: 0.35, id: 16_777_217 },
      { c: [-1.1, -1.3, -8] as [number, number, number], r: 0.35, id: 4_294_967_295 },
      { c: [1.1, -1.3, -8] as [number, number, number], r: 0.35, id: 999_999 },
    ];

    const n = BACKGROUND + probes.length;
    const centers = new Float32Array(n * 3);
    const radii = new Float32Array(n);
    const ids = new Uint32Array(n);
    const colors = new Float32Array(n * 3);

    for (let i = 0; i < BACKGROUND; i++) {
      centers[i * 3] = (Math.random() - 0.5) * 14;
      centers[i * 3 + 1] = (Math.random() - 0.5) * 14;
      centers[i * 3 + 2] = -26 - Math.random() * 6; // franchement derrière les sondes
      radii[i] = 0.03;
      ids[i] = 1_000_000 + i;
      colors[i * 3] = 0.2;
      colors[i * 3 + 1] = 0.25;
      colors[i * 3 + 2] = 0.3;
    }
    probes.forEach((p, k) => {
      const i = BACKGROUND + k;
      centers.set(p.c, i * 3);
      radii[i] = p.r;
      ids[i] = p.id;
      colors.set([1, 0.6, 0.2], i * 3);
    });

    const r = drawScene({ centers, radii, ids, colors, count: n });
    const bad = probes.filter((p) => {
      const [x, y] = project(p.c);
      return picker.pick(x, y) !== p.id;
    });
    check(
      `picking exact parmi ${n.toLocaleString("fr")} instances`,
      bad.length === 0,
      bad.length ? `manqué : ${bad.map((p) => p.id).join(", ")}` : `${probes.length}/${probes.length} sondes`,
    );

    const big = probes[3];
    const [bx, by] = project(big.c);
    check(
      "un identifiant occupant les 32 bits survit au transport",
      picker.pick(bx, by) === 4_294_967_295,
      `attendu 4294967295, obtenu ${picker.pick(bx, by)}`,
    );
    check(
      "les billes du fond n'occultent pas les sondes",
      picker.pick(...project(probes[0].c)) === 1,
      "",
    );
    r.dispose();
  }

  renderer.dispose();
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
