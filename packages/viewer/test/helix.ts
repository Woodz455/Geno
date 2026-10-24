/** Suite de correction de la semaine 12 : le chemin de l'ADN et la double hélice.
 *
 *  Comme les suites précédentes : chaque propriété qui compte a son **témoin**, construit
 *  pour la violer. La géométrie est vérifiée sans GPU ; le rendu, dans un vrai navigateur.
 */

import {
  B_DNA,
  BASE_COLORS,
  baseOf,
  complement,
  createTargets,
  DnaPath,
  HelixRenderer,
  MIN_SPACING_NM,
  NO_HIT,
  PER_MONOMER,
  perspective,
  phosphates,
  Picker,
  RISE_NM,
  SUPER_PITCH_NM,
  SUPER_RADIUS_NM,
  type Targets,
  WRAP_BP,
  WRAP_TURNS,
} from "../src/index.js";

const W = 320;
const H = 320;
const FOV = (45 * Math.PI) / 180;
const IDENT = new Float32Array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]);

interface Result {
  name: string;
  ok: boolean;
  detail: string;
}
const results: Result[] = [];
const check = (name: string, ok: boolean, detail = "") => results.push({ name, ok, detail });

type V = [number, number, number];
const sub = (a: ArrayLike<number>, b: ArrayLike<number>): V => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a: ArrayLike<number>, b: ArrayLike<number>) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a: ArrayLike<number>, b: ArrayLike<number>): V => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const len = (a: ArrayLike<number>) => Math.hypot(a[0], a[1], a[2]);
const deg = (r: number) => (r * 180) / Math.PI;
const TWIST = 360 / B_DNA.bpPerTurn;

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

/** Une fibre synthétique : marche de pas σ = 46,3 nm, en µm comme le vrai niveau, à volume
 *  exclu comme la chaîne fine (deux monomères non liés restent à 0,85 σ au moins). */
function fibre(n: number, seed = 5): Float32Array {
  const r = mulberry32(seed);
  const pts = new Float32Array(n * 3);
  let p: V = [0, 0, 0];
  let d: V = [1, 0, 0];
  for (let k = 0; k < n; k++) {
    pts.set(p, 3 * k);
    for (let attempt = 0; attempt < 200; attempt++) {
      // Persistance modérée, comme la chaîne fine (angle de liaison ~60°).
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

/** Angle signé de `a` à `b` autour de `t`. */
const signedAngle = (a: V, b: V, t: V) => Math.atan2(dot(cross(a, b), t), dot(a, b));

function pathSuite(): void {
  const N = 30;
  const pts = fibre(N);
  const path = new DnaPath(pts, 2000, 23.1);
  const nm = 0.001;
  const total = path.bpCount;
  const bps = Array.from({ length: total }, (_, b) => path.basePair(b));
  const inCore = new Uint8Array(total);
  const coreOf = new Int32Array(total).fill(-1);
  const cores = path.cores(Array.from({ length: N }, (_, k) => k));
  for (const c of cores) for (let j = 0; j < WRAP_BP; j++) {
    inCore[c.firstBp + j] = 1;
    coreOf[c.firstBp + j] = c.index;
  }

  // 1. Le pas : 0,34 nm sur les liens, et ce qu'impose le cœur (1,67 tour, rayon 4,18 nm) dessus.
  const wrapStep = (WRAP_TURNS * Math.hypot(2 * Math.PI * SUPER_RADIUS_NM, SUPER_PITCH_NM)) / (WRAP_BP - 1);
  let devLink = 0, devWrap = 0, maxJump = 0, maxBend = 0, bent = 0;
  for (let b = 0; b + 1 < total; b++) {
    const s = len(sub(bps[b + 1].p, bps[b].p)) / nm;
    maxJump = Math.max(maxJump, s);
    const same = inCore[b] && inCore[b + 1] && coreOf[b] === coreOf[b + 1];
    if (same) devWrap = Math.max(devWrap, Math.abs(s - wrapStep) / wrapStep);
    else devLink = Math.max(devLink, Math.abs(s - RISE_NM) / RISE_NM);
    const bend = deg(Math.acos(Math.max(-1, Math.min(1, dot(bps[b].t, bps[b + 1].t)))));
    maxBend = Math.max(maxBend, bend);
    if (bend > 20) bent++;
  }
  check(
    "pas d'une base à la suivante : 0,34 nm sur les liens, à 1 % près",
    devLink < 0.01,
    `écart max ${(devLink * 100).toFixed(2)} % sur ${total.toLocaleString("fr")} pas`,
  );
  check(
    `pas sur l'ADN enroulé : ${wrapStep.toFixed(4)} nm, imposé par 1,67 tour de 4,18 nm`,
    devWrap < 0.01,
    `écart max ${(devWrap * 100).toFixed(3)} % — une hélice B de 147 pb à 0,34 nm serait plus longue que l'enroulement`,
  );
  check(
    "aucun saut ni coude : chemin continu à la base près, jonctions et entrées de cœur comprises",
    maxJump < 0.35 && maxBend < 20,
    `plus grand pas ${maxJump.toFixed(3)} nm, plus grand coude ${maxBend.toFixed(1)}° d'une base à la suivante, ${bent} pas au-delà de 20°`,
  );
  check("aucun lien plus court que sa corde", path.infeasible === 0, `${path.infeasible} lien(s) impossible(s)`);

  // 2. Le twist : 360/10,5 par pas, dans le sens droit ; sur un lien, l'écart réparti est borné.
  const twist: number[] = [];
  let wrapDev = 0, linkOver = 0, right = 0;
  for (let b = 0; b + 1 < total; b++) {
    const t = bps[b + 1].t;
    // La dyade de b transportée sur la tangente suivante, puis l'angle autour de celle-ci.
    const d0 = sub(bps[b].d, [t[0] * dot(bps[b].d, t), t[1] * dot(bps[b].d, t), t[2] * dot(bps[b].d, t)]);
    const a = deg(signedAngle(d0, bps[b + 1].d as V, t as V));
    twist.push(a);
    if (a > 0) right++;
    const same = inCore[b] && inCore[b + 1] && coreOf[b] === coreOf[b + 1];
    if (same) wrapDev = Math.max(wrapDev, Math.abs(a - TWIST));
    else if (Math.abs(a - TWIST) > 180 / 8 + 1) linkOver++;
  }
  const mean = twist.reduce((s, x) => s + x, 0) / twist.length;
  check(
    "twist sur l'ADN enroulé : 34,29° par pas, au transport près de la superhélice",
    wrapDev < 0.6,
    `écart max ${wrapDev.toFixed(3)}° (la torsion de la superhélice en ajoute ~0,4)`,
  );
  check(
    "twist sur les liens : l'écart réparti reste sous sa borne, et la moyenne est 34,29°",
    linkOver === 0 && Math.abs(mean - TWIST) < 0.5,
    `moyenne ${mean.toFixed(3)}° sur ${twist.length.toLocaleString("fr")} pas, ${linkOver} au-delà de la borne`,
  );
  check("hélice droite partout", right === twist.length, `${right}/${twist.length} pas tournent à droite`);

  // 3. La superhélice est gauche : l'ADN tourne autour de l'axe du cœur dans un sens, et recule
  //    le long de l'axe dans l'autre.
  let left = 0;
  for (const c of cores) {
    const p0 = bps[c.firstBp + 30].p, p1 = bps[c.firstBp + 31].p;
    const r0 = sub(p0, c.center), r1 = sub(p1, c.center);
    const turn = dot(cross(r0, r1), c.axis);
    const climb = dot(sub(p1, p0), c.axis);
    if (turn * climb < 0) left++;
  }
  check("superhélice gauche sur chaque cœur", left === cores.length, `${left}/${cores.length}`);

  // 4. Les cœurs : 10 par monomère, à 10 nm au moins dans un monomère ; entre monomères voisins,
  //    rien n'est vérifié à la construction — on compte.
  let inside = 0, across = 0;
  for (let a = 0; a < cores.length; a++) {
    for (let b = a + 1; b < cores.length; b++) {
      if (len(sub(cores[a].center, cores[b].center)) / nm < MIN_SPACING_NM) {
        if (cores[a].monomer === cores[b].monomer) inside++;
        else across++;
      }
    }
  }
  check(
    `${PER_MONOMER} cœurs par monomère, aucun à moins de 10 nm d'un autre du même monomère`,
    cores.length === N * PER_MONOMER && inside === 0,
    `${cores.length} cœurs, ${inside} paire(s) trop proches dans un monomère ; ${across} entre monomères différents`,
  );
  check("et pas plus de 2 % de paires trop proches entre monomères", across <= 0.02 * cores.length, `${across} sur ${cores.length} cœurs`);

  // 5. Paresseux et déterministe : l'ordre de calcul ne change rien.
  const again = new DnaPath(pts, 2000, 23.1);
  const probe = [59_999, 31_337, 0, 2_000, 1_999, 12_345];
  let same = true;
  for (const b of probe) {
    const x = again.basePair(b), y = bps[b];
    same &&= len(sub(x.p, y.p)) < 1e-9 && len(sub(x.d, y.d)) < 1e-6;
  }
  check("déterministe : autre instance, autre ordre de calcul, mêmes bases", same);
  const other = new DnaPath(pts, 2000, 23.1, 0.001, 12);
  check("témoin : une autre graine place d'autres nucléosomes", len(sub(other.basePair(12_345).p, bps[12_345].p)) / nm > 1);
}

/** Une hélice droite le long de y, en nm, au centre du champ. */
function straight(n: number, z = -12): Float32Array {
  const f = new Float32Array(n * 9);
  const a = (2 * Math.PI) / B_DNA.bpPerTurn;
  for (let b = 0; b < n; b++) {
    const phi = a * b;
    f.set([0, (b - n / 2) * B_DNA.riseNm, z, 0, 1, 0, Math.cos(phi), 0, -Math.sin(phi)], 9 * b);
  }
  return f;
}

function grooveSuite(): void {
  const n = 80;
  const f = straight(n, 0);
  const { s1, s2 } = phosphates(f, 1);
  const P = (s: Float32Array, i: number): V => [s[3 * i], s[3 * i + 1], s[3 * i + 2]];
  // Rayon des phosphates et pas : exacts par construction.
  check("phosphates à 8,91 Å de l'axe", Math.abs(Math.hypot(s1[0], s1[2]) - 0.891) < 1e-5);
  const pitch = (P(s1, 50)[1] - P(s1, 29)[1]) / 2;
  check("pas de l'hélice : 10,5 pb, 3,57 nm", Math.abs(pitch - 3.57) < 1e-4, `${pitch.toFixed(4)} nm`);
  // Largeurs de sillon : du phosphate 40 du brin 1, le P le plus proche du brin 2 au-dessus et
  // au-dessous, moins un diamètre de phosphate.
  const p = P(s1, 40);
  let up = Infinity, down = Infinity;
  for (let i = 0; i < n; i++) {
    const q = P(s2, i);
    const d = len(sub(q, p));
    if (q[1] > p[1]) up = Math.min(up, d);
    else down = Math.min(down, d);
  }
  const [minor, major] = [up, down].map((x) => x - B_DNA.phosphateDiameterNm).sort((a, b) => a - b);
  check(
    "petit sillon ~5,7 Å, grand ~11,7 Å (valeurs publiées pour la B)",
    minor > 0.45 && minor < 0.65 && major > 1.05 && major < 1.3,
    `mesurés ${(minor * 10).toFixed(1)} Å et ${(major * 10).toFixed(1)} Å`,
  );
  // Témoin : sans le décalage axial ±2,85 Å, les deux sillons se rapprochent.
  const saved = B_DNA.phosphateRiseNm;
  B_DNA.phosphateRiseNm = 0;
  const flat = phosphates(f, 1);
  B_DNA.phosphateRiseNm = saved;
  let u2 = Infinity, d2 = Infinity;
  const p2 = P(flat.s1, 40);
  for (let i = 0; i < n; i++) {
    const q = P(flat.s2, i);
    const d = len(sub(q, p2));
    if (q[1] > p2[1]) u2 = Math.min(u2, d);
    else d2 = Math.min(d2, d);
  }
  const ratio = Math.max(u2, d2) / Math.min(u2, d2);
  check("témoin : sans décalage axial, les sillons s'égalisent", ratio < (major + 0.58) / (minor + 0.58), `rapport ${ratio.toFixed(2)} contre ${((major + 0.58) / (minor + 0.58)).toFixed(2)}`);
  // Un barreau ne traverse jamais le squelette de l'autre brin, hors de ses propres extrémités.
  const segDist = (a: V, b: V, c: V, d: V) => {
    let best = Infinity;
    for (let i = 0; i <= 20; i++) {
      const x: V = [a[0] + (b[0] - a[0]) * (i / 20), a[1] + (b[1] - a[1]) * (i / 20), a[2] + (b[2] - a[2]) * (i / 20)];
      for (let j = 0; j <= 20; j++) {
        const y: V = [c[0] + (d[0] - c[0]) * (j / 20), c[1] + (d[1] - c[1]) * (j / 20), c[2] + (d[2] - c[2]) * (j / 20)];
        best = Math.min(best, len(sub(x, y)));
      }
    }
    return best;
  };
  let crossing = 0;
  for (let i = 5; i < n - 5; i++) {
    const a = P(s1, i), b = P(s2, i);
    // Le barreau privé de ses 0,3 nm d'extrémité, contre les squelettes non voisins.
    const ab = sub(b, a);
    const l = len(ab);
    const a1: V = [a[0] + (ab[0] / l) * 0.3, a[1] + (ab[1] / l) * 0.3, a[2] + (ab[2] / l) * 0.3];
    const b1: V = [b[0] - (ab[0] / l) * 0.3, b[1] - (ab[1] / l) * 0.3, b[2] - (ab[2] / l) * 0.3];
    for (const s of [s1, s2]) {
      for (let j = 0; j + 1 < n; j++) {
        if (Math.abs(j - i) < 2) continue;
        if (segDist(a1, b1, P(s, j), P(s, j + 1)) < 0.21 + 0.16) crossing++;
      }
    }
  }
  check("aucun barreau ne traverse un squelette", crossing === 0, `${crossing} contact(s)`);
}

function clear(gl: WebGL2RenderingContext, t: Targets): void {
  gl.bindFramebuffer(gl.FRAMEBUFFER, t.fbo);
  gl.viewport(0, 0, W, H);
  gl.clearBufferfv(gl.COLOR, 0, [0, 0, 0, 1]);
  gl.clearBufferuiv(gl.COLOR, 1, [NO_HIT, 0, 0, 0]);
  gl.clearBufferfv(gl.COLOR, 2, [1, 0, 0, 1]);
  gl.clearBufferfi(gl.DEPTH_STENCIL, 0, 1.0, 0);
}

function read(gl: WebGL2RenderingContext, t: Targets): { ids: Uint32Array; rgb: Uint8Array } {
  const buf = new Uint32Array(W * H * 4);
  gl.bindFramebuffer(gl.READ_FRAMEBUFFER, t.fbo);
  gl.readBuffer(gl.COLOR_ATTACHMENT1);
  gl.readPixels(0, 0, W, H, gl.RGBA_INTEGER, gl.UNSIGNED_INT, buf);
  gl.readBuffer(gl.COLOR_ATTACHMENT0);
  const rgb = new Uint8Array(W * H * 4);
  gl.readPixels(0, 0, W, H, gl.RGBA, gl.UNSIGNED_BYTE, rgb);
  return { ids: buf.filter((_, i) => i % 4 === 0), rgb };
}

function hue(r: number, g: number, b: number): { h: number; s: number } {
  const mx = Math.max(r, g, b), mn = Math.min(r, g, b);
  const d = mx - mn;
  if (d === 0) return { h: 0, s: 0 };
  let h = mx === r ? ((g - b) / d) % 6 : mx === g ? (b - r) / d + 2 : (r - g) / d + 4;
  h *= 60;
  if (h < 0) h += 360;
  return { h, s: d / mx };
}

const HUES = Object.fromEntries(
  (["A", "C", "G", "T"] as const).map((b) => [b, hue(...(BASE_COLORS[b].map((x) => x * 255) as [number, number, number])).h]),
);
function classify(r: number, g: number, b: number): string | null {
  const { h, s } = hue(r, g, b);
  if (s < 0.35) return null;           // squelette gris, ou ombre
  let best = "N", gap = Infinity;
  for (const [k, v] of Object.entries(HUES)) {
    const dh = Math.min(Math.abs(h - v), 360 - Math.abs(h - v));
    if (dh < gap) {
      gap = dh;
      best = k;
    }
  }
  return best;
}

function renderSuite(gl: WebGL2RenderingContext, targets: Targets, picker: Picker): void {
  const n = 40;
  const r = mulberry32(3);
  const seq = Array.from({ length: n }, () => "ACGT"[Math.floor(r() * 4)]).join("");
  const helix = new HelixRenderer(gl, 1);
  const FIRST = 1_000_000;          // rang dans une région fictive : l'identifiant doit le porter
  helix.upload({ ranges: [[FIRST, n]], frames: straight(n), bases: seq });
  const proj = perspective(FOV, W / H, 1, 40);
  clear(gl, targets);
  helix.draw(proj, { view: IDENT, height: H, idBase: 0 });
  const { ids, rgb } = read(gl, targets);

  const seen = new Set<number>();
  let outside = 0;
  for (const id of ids) {
    if (id === NO_HIT) continue;
    const { rank } = baseOf(id);
    if (rank < FIRST || rank >= FIRST + n) outside++;
    else seen.add(id);
  }
  check("chaque pixel nomme une base de la fenêtre, et seulement elle", outside === 0 && seen.size >= 50, `${seen.size} bases distinctes visibles, ${outside} pixel(s) hors fenêtre`);

  // Fidélité : la couleur de chaque pixel de barreau, classée par teinte, contre la base que
  // son identifiant désigne dans la séquence — brin + tel quel, brin − complémenté.
  const expected = (id: number, shift = 0) => {
    const { rank, strand } = baseOf(id);
    const q = rank - FIRST + shift;
    if (q < 0 || q >= n) return null;
    return strand === 0 ? seq[q] : complement(seq[q]);
  };
  let agree = 0, total = 0, agreeShift = 0, totalShift = 0;
  for (let i = 0; i < ids.length; i++) {
    if (ids[i] === NO_HIT) continue;
    const got = classify(rgb[4 * i], rgb[4 * i + 1], rgb[4 * i + 2]);
    if (!got) continue;
    const e = expected(ids[i]);
    if (e) {
      total++;
      if (e === got) agree++;
    }
    const s = expected(ids[i], 1);
    if (s) {
      totalShift++;
      if (s === got) agreeShift++;
    }
  }
  check("chaque barreau a la couleur de sa base", total > 500 && agree / total > 0.99, `${((agree / total) * 100).toFixed(2)} % de ${total} pixels`);
  check("témoin : la séquence décalée d'une base ne colle plus", agreeShift / totalShift < 0.5, `${((agreeShift / totalShift) * 100).toFixed(1)} %`);

  // Picking : chaque pixel cliqué rend la base qu'un lancer de rayon **analytique** (cylindres
  // idéaux des barreaux et des squelettes, sans GPU) trouve la plus proche. La première
  // version visait le milieu de chaque demi-barreau et comptait faux quand un voisin, à
  // 0,34 nm, le cachait légitimement : c'était l'oracle qui se trompait, pas le picking.
  const { s1, s2 } = phosphates(straight(n), 1);
  const caps: { a: number[]; b: number[]; r: number; id: number; rung: boolean }[] = [];
  for (let q = 0; q < n; q++) {
    const a = [s1[3 * q], s1[3 * q + 1], s1[3 * q + 2]];
    const b = [s2[3 * q], s2[3 * q + 1], s2[3 * q + 2]];
    const m = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2];
    caps.push({ a, b: m, r: 0.16, id: 2 * (FIRST + q) + 1, rung: true }, { a: b, b: m, r: 0.16, id: 2 * (FIRST + q) + 2, rung: true });
    if (q + 1 < n) {
      // Un segment de squelette appartient pour moitié à chacun de ses deux points : c'est la
      // règle du tube (« un pixel appartient au point le plus proche »).
      for (const [s, st] of [[s1, 0], [s2, 1]] as const) {
        const p0 = [s[3 * q], s[3 * q + 1], s[3 * q + 2]];
        const p1 = [s[3 * q + 3], s[3 * q + 4], s[3 * q + 5]];
        const mid = [(p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2, (p0[2] + p1[2]) / 2];
        caps.push({ a: p0, b: mid, r: 0.21, id: 2 * (FIRST + q) + st + 1, rung: false }, { a: mid, b: p1, r: 0.21, id: 2 * (FIRST + q + 1) + st + 1, rung: false });
      }
    }
  }
  // Rayon depuis l'origine vers le pixel ; profondeur d'entrée approchée dans chaque cylindre.
  const hits = (x: number, y: number) => {
    const f = 1 / Math.tan(FOV / 2);
    const dir = [((x + 0.5) / W * 2 - 1) / f, (1 - (y + 0.5) / H * 2) / f, -1];
    const dl = Math.hypot(dir[0], dir[1], dir[2]);
    const d = dir.map((v) => v / dl);
    const out: { t: number; id: number; depth: number; rung: boolean }[] = [];
    for (const c of caps) {
      const u = sub(c.b, c.a);
      const w0 = [-c.a[0], -c.a[1], -c.a[2]];
      const A = dot(d, d), B = dot(d, u), C = dot(u, u), D = dot(d, w0), E = dot(u, w0);
      const den = A * C - B * B;
      let sc = den > 1e-12 ? (B * E - C * D) / den : 0;
      let tc = den > 1e-12 ? (A * E - B * D) / den : 0;
      tc = Math.max(0, Math.min(1, tc));
      sc = Math.max(0, (B * tc - D) / A);
      const pr = [d[0] * sc, d[1] * sc, d[2] * sc];
      const ps = [c.a[0] + u[0] * tc, c.a[1] + u[1] * tc, c.a[2] + u[2] * tc];
      const dist = len(sub(pr, ps));
      if (dist <= c.r) out.push({ t: sc - Math.sqrt(c.r * c.r - dist * dist), id: c.id, depth: dist / c.r, rung: c.rung });
    }
    return out.sort((p, q) => p.t - q.t);
  };
  let exact = 0, tried = 0, spine = 0, spineOk = 0;
  const misses: string[] = [];
  for (let y = 0; y < H; y += 5) {
    for (let x = 0; x < W; x += 5) {
      const h = hits(x, y);
      if (!h.length) continue;
      // Pixels sans ambiguïté seulement : bien à l'intérieur du cylindre touché (le tube rendu
      // est un polygone), et rien d'autre à moins de 0,05 nm derrière.
      if (h[0].depth > 0.8 || (h.length > 1 && h[1].t - h[0].t < 0.05 && h[1].id !== h[0].id)) continue;
      const got = picker.pick(x, y);
      // Un segment de squelette entre deux phosphates n'appartient à une base que par
      // convention ; le tube rendu est une spline, l'oracle un segment droit. On le compte à
      // part : ce qui doit être exact, c'est le barreau d'une base.
      if (!h[0].rung) {
        spine++;
        if (got === h[0].id) spineOk++;
        continue;
      }
      tried++;
      if (got === h[0].id) exact++;
      else if (misses.length < 6) misses.push(`(${x},${y}) ${baseOf(h[0].id).rank - FIRST}/${baseOf(h[0].id).strand}→${got === NO_HIT ? "rien" : `${baseOf(got).rank - FIRST}/${baseOf(got).strand}`}`);
    }
  }
  check(
    "picking exact à la base : un clic sur un barreau rend sa base (lancer de rayon analytique)",
    tried >= 200 && exact / tried >= 0.99,
    `${exact}/${tried} pixels de barreau${misses.length ? ` ; manqués ${misses.join(", ")}` : ""} · squelette, par convention : ${spineOk}/${spine}`,
  );
}

function run(): void {
  pathSuite();
  grooveSuite();
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
  renderSuite(gl, targets, picker);
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
