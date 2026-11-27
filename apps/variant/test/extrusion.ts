/** Le port TypeScript de l'extrusion 1D contre Python, en statistique (les générateurs
 *  aléatoires diffèrent). Exécuté sous Node par `tools/test-variant.mjs`. */

import { readFileSync } from "node:fs";
import { mean, rows, simulate } from "../extrusion.js";

const [, , expectedPath] = process.argv;
const E = JSON.parse(readFileSync(expectedPath, "utf8")) as {
  n: number; params: any; seeds: number; law: { alpha: number; d0: number };
  sites: { pos: number[]; strand: number[]; occupancy: number[] };
  sources: number[]; landmarks: Record<string, number[]>; foot: number[]; rows: number[][]; occupancy_lef: number;
  foot_b: number[]; rows_b: number[][]; batches: Record<string, number[]>;
};

let passed = 0;
let failed = 0;
const check = (name: string, ok: boolean, detail = "") => {
  if (ok) passed++;
  else failed++;
  console.log(`  ${ok ? "\x1b[32mPASS\x1b[0m" : "\x1b[31mFAIL\x1b[0m"}  ${name}${detail ? `  \x1b[2m${detail}\x1b[0m` : ""}`);
};
const pearson = (a: ArrayLike<number>, b: ArrayLike<number>) => {
  const n = a.length;
  let ma = 0, mb = 0;
  for (let i = 0; i < n; i++) { ma += a[i]; mb += b[i]; }
  ma /= n; mb /= n;
  let sab = 0, saa = 0, sbb = 0;
  for (let i = 0; i < n; i++) { const x = a[i] - ma, y = b[i] - mb; sab += x * y; saa += x * x; sbb += y * y; }
  return sab / Math.sqrt(saa * sbb);
};

const sites = { pos: Int32Array.from(E.sites.pos), strand: Int8Array.from(E.sites.strand), occupancy: Float32Array.from(E.sites.occupancy) };
const n = E.n;
const t0 = performance.now();
const runs = Array.from({ length: E.seeds }, (_, s) => simulate(n, sites, E.params, 70_000 + s));
const tSim = performance.now() - t0;

const foot = new Float64Array(n);
let loaded = 0, slots = 0;
for (const r of runs) {
  for (let f = 0; f < r.snapshots; f++) for (let k = 0; k < r.nLef; k++) {
    const i = r.legs[(f * r.nLef + k) * 2], j = r.legs[(f * r.nLef + k) * 2 + 1];
    slots++;
    if (i >= 0) { foot[i] += 1; foot[j] += 1; loaded++; }
  }
}
for (let i = 0; i < n; i++) foot[i] /= runs.length * runs[0].snapshots;
// Le juge est le bruit : deux lots Python de graines disjointes ne s'accordent pas mieux que
// ça, et TypeScript contre Python ne doit pas faire nettement pire (marge : 1,5 × l'écart).
const rf = pearson(foot, E.foot);
const rfPy = pearson(E.foot_b, E.foot);
check("densité de pieds : TS contre Python, pas pire que Python contre Python", 1 - rf <= 1.5 * (1 - rfPy),
      `Pearson TS–Py ${rf.toFixed(3)} · Py–Py ${rfPy.toFixed(3)} sur ${n} monomères, ${E.seeds} exécutions par lot`);
const occ = loaded / slots;
check("occupation des cohésines", Math.abs(occ - E.occupancy_lef) <= 0.02, `TS ${occ.toFixed(3)} · Python ${E.occupancy_lef.toFixed(3)}`);

const t1 = performance.now();
const R = rows(n, runs, E.sources, E.law);
const tRows = performance.now() - t1;
const flat = Float64Array.from(E.rows.flat());
const lr = pearson(Array.from(R, (x) => Math.log(x)), Array.from(flat, (x) => Math.log(x)));
check("contacts prédits (lignes depuis l'enhancer et PAX3) : TS ≈ Python", lr >= 0.98, `Pearson des log ${lr.toFixed(4)}`);

const at = new Map(E.sources.map((s, k) => [s, k]));
const enh = E.landmarks.enhancer.filter((i) => at.has(i));
const pick = (M: ArrayLike<number>, a: number[], b: number[]) => {
  const sub = new Float32Array(a.length * n);
  a.forEach((i, k) => { for (let j = 0; j < n; j++) sub[k * n + j] = M[at.get(i)! * n + j]; });
  return mean(sub, n, a.length, b);
};
// Contacts enhancer → promoteur : 16 lots de 40 exécutions de chaque côté, test de Welch. Un
// écart systématique de 3,6 % (mulberry32) passait inaperçu sur un lot et ressortait ici.
const K = E.batches.PAX3.length;
const ts: Record<string, number[]> = { EPHA4: [], PAX3: [], IHH: [] };
for (let b = 0; b < K; b++) {
  const rr = Array.from({ length: E.seeds }, (_, s) => simulate(n, sites, E.params, 500_000 + 1000 * b + s));
  const Rb = rows(n, rr, E.sources, E.law);
  for (const g of Object.keys(ts)) ts[g].push(pick(Rb, enh, E.landmarks[g]));
}
const stat = (x: number[]) => { const m = x.reduce((a, b) => a + b, 0) / x.length; return [m, x.reduce((a, b) => a + (b - m) ** 2, 0) / (x.length - 1)]; };
for (const g of Object.keys(ts)) {
  const [ma, va] = stat(ts[g]), [mb, vb] = stat(E.batches[g]);
  const t = (ma - mb) / Math.sqrt(va / K + vb / K);
  // |t| < 2,75 : p > 0,01 pour une trentaine de degrés de liberté.
  check(`enhancer → ${g} : TS et Python indiscernables (Welch, ${K} lots chacun)`, Math.abs(t) < 2.75,
        `TS ${ma.toExponential(3)} · Python ${mb.toExponential(3)} · écart ${((ma / mb - 1) * 100).toFixed(1)} % · t = ${t.toFixed(2)}`);
}
console.log(`\n  ${E.seeds} extrusions en ${tSim.toFixed(0)} ms, ${E.sources.length} lignes de contact en ${tRows.toFixed(0)} ms (Node)`);
console.log(`\n${passed} passés, ${failed} échoués`);
process.exit(failed ? 1 : 0);
