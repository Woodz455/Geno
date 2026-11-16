/** L'extrusion de boucles en 1D, et les contacts qu'elle prédit — dans le navigateur.
 *
 *  Port de `pipeline/geno_pipeline/fine/extrusion.py` (Fudenberg 2016, cinq règles) et de
 *  `sv/predict.py` (distance effective par les boucles tenues, loi `min(1, (d/d0)^−α + c)` calée
 *  sur le 3D). Mêmes règles, même ordre des opérations par pas ; le générateur aléatoire diffère
 *  (numpy n'est pas là), donc l'accord avec Python se mesure en statistique, pas au tirage près
 *  (`tools/test-variant.mjs`).
 */

export interface Params {
  separation: number;      // pb par cohésine
  processivity: number;    // pb
  velocity: number;
  release: number;
  burnin: number;
  snapshots: number;
  stride: number;
  bpPerBead: number;
}

export interface SiteSet {
  pos: Int32Array;          // index de monomère, sur l'allèle simulé
  strand: Int8Array;        // +1 forward, −1 reverse
  occupancy: Float32Array;
}

/** sfc32, graine en quatrième mot, 15 tirages de mise en route. mulberry32 (un seul mot d'état)
 *  donnait des contacts biaisés de −3,6 % par rapport à Python sur 20 lots × 40 exécutions
 *  (Welch p = 0,004) ; sfc32, aucun écart significatif (VALIDATION § S15). */
export function rng(seed: number): () => number {
  let a = 0x9e3779b9, b = 0x243f6a88, c = 0xb7e15162, d = seed >>> 0;
  const next = () => {
    const t = (((a + b) >>> 0) + d) >>> 0;
    d = (d + 1) >>> 0;
    a = b ^ (b >>> 9);
    b = (c + (c << 3)) >>> 0;
    c = (c << 21) | (c >>> 11);
    c = (c + t) >>> 0;
    return t >>> 0;
  };
  for (let i = 0; i < 15; i++) next();
  return () => next() / 4294967296;
}

export function barriers(n: number, s: SiteSet): [Float32Array, Float32Array] {
  const left = new Float32Array(n);
  const right = new Float32Array(n);
  for (let k = 0; k < s.pos.length; k++) {
    const t = s.strand[k] > 0 ? left : right;
    t[s.pos[k]] = Math.max(t[s.pos[k]], s.occupancy[k]);
  }
  return [left, right];
}

/** Les pieds des cohésines, instantané par instantané : (snapshots × nLef × 2), −1 = décroché. */
export function simulate(n: number, sites: SiteSet, p: Params, seed: number): { legs: Int32Array; nLef: number; snapshots: number } {
  const R = rng(seed);
  const [stopL, stopR] = barriers(n, sites);
  const nLef = Math.max(1, Math.round((n * p.bpPerBead) / p.separation));
  const lifetime = Math.max(1, p.processivity / (2 * p.velocity * p.bpPerBead));
  const pUnload = 1 / lifetime;
  const left = new Int32Array(nLef).fill(-1);
  const right = new Int32Array(nLef).fill(-1);
  const heldL = new Uint8Array(nLef);
  const heldR = new Uint8Array(nLef);
  const busy = new Uint8Array(n);

  const unload = (k: number) => {
    if (left[k] >= 0) {
      busy[left[k]] = 0;
      busy[right[k]] = 0;
    }
    left[k] = right[k] = -1;
    heldL[k] = heldR[k] = 0;
  };
  const load = (k: number) => {
    for (let t = 0; t < 8; t++) {
      const i = Math.floor(R() * (n - 1));
      if (!busy[i] && !busy[i + 1]) {
        left[k] = i;
        right[k] = i + 1;
        busy[i] = busy[i + 1] = 1;
        heldL[k] = R() < stopL[i] ? 1 : 0;
        heldR[k] = R() < stopR[i + 1] ? 1 : 0;
        return;
      }
    }
  };
  const advance = () => {
    for (let k = 0; k < nLef; k++) if (left[k] >= 0 && R() < pUnload) unload(k);
    for (let k = 0; k < nLef; k++) if (left[k] < 0) load(k);
    for (let k = 0; k < nLef; k++) {
      if (heldL[k] && R() < p.release) heldL[k] = 0;
      if (heldR[k] && R() < p.release) heldR[k] = 0;
    }
    for (let k = 0; k < nLef; k++) {
      if (left[k] < 0) continue;
      const ml = R() < p.velocity;
      const mr = R() < p.velocity;
      if (ml && !heldL[k]) {
        const d = left[k] - 1;
        if (d >= 0 && !busy[d]) {
          busy[left[k]] = 0;
          busy[d] = 1;
          left[k] = d;
          if (R() < stopL[d]) heldL[k] = 1;
        }
      }
      if (mr && !heldR[k]) {
        const d = right[k] + 1;
        if (d < n && !busy[d]) {
          busy[right[k]] = 0;
          busy[d] = 1;
          right[k] = d;
          if (R() < stopR[d]) heldR[k] = 1;
        }
      }
    }
  };
  for (let t = 0; t < p.burnin; t++) advance();
  const legs = new Int32Array(p.snapshots * nLef * 2);
  for (let f = 0; f < p.snapshots; f++) {
    for (let t = 0; t < Math.max(1, p.stride); t++) advance();
    for (let k = 0; k < nLef; k++) {
      legs[(f * nLef + k) * 2] = left[k];
      legs[(f * nLef + k) * 2 + 1] = right[k];
    }
  }
  return { legs, nLef, snapshots: p.snapshots };
}

/** Plus court chemin (parcours en largeur, arêtes de longueur 1) depuis `src`, chaîne + boucles
 *  tenues de l'instantané `f`. `dist`, `queue`, `jumpTo` sont des tampons réutilisés. */
function bfs(n: number, legs: Int32Array, nLef: number, f: number, src: number, dist: Int32Array, queue: Int32Array,
             jumpTo: Int32Array): void {
  // Les boucles de l'instantané : un monomère porte au plus un pied (les jambes ne se
  // chevauchent pas), donc un saut au plus.
  jumpTo.fill(-1);
  for (let k = 0; k < nLef; k++) {
    const i = legs[(f * nLef + k) * 2], j = legs[(f * nLef + k) * 2 + 1];
    if (i >= 0 && j >= 0 && i !== j) {
      jumpTo[i] = j;
      jumpTo[j] = i;
    }
  }
  dist.fill(-1);
  let head = 0, tail = 0;
  dist[src] = 0;
  queue[tail++] = src;
  while (head < tail) {
    const u = queue[head++];
    const du = dist[u] + 1;
    if (u > 0 && dist[u - 1] < 0) { dist[u - 1] = du; queue[tail++] = u - 1; }
    if (u < n - 1 && dist[u + 1] < 0) { dist[u + 1] = du; queue[tail++] = u + 1; }
    const v = jumpTo[u];
    if (v >= 0 && dist[v] < 0) { dist[v] = du; queue[tail++] = v; }
  }
}

export interface Law {
  alpha: number;
  d0: number;
  /** Plancher : le mélange dû au confinement (`sv/predict.py`). */
  c?: number;
}

/** Lignes de contact prédit : (sources × n), moyennées sur les instantanés de toutes les
 *  exécutions données. */
export function rows(n: number, runs: { legs: Int32Array; nLef: number; snapshots: number }[], sources: number[], law: Law): Float32Array {
  const out = new Float32Array(sources.length * n);
  const dist = new Int32Array(n);
  const queue = new Int32Array(n);
  const jump = new Int32Array(n);
  // P(d) en table : les distances sont entières et bornées par n.
  const table = new Float32Array(n + 1);
  for (let d = 0; d <= n; d++) table[d] = d === 0 ? 1 : Math.min(1, Math.pow(d / law.d0, -law.alpha) + (law.c ?? 0));
  let frames = 0;
  for (const r of runs) {
    for (let f = 0; f < r.snapshots; f++) {
      frames++;
      sources.forEach((s, k) => {
        bfs(n, r.legs, r.nLef, f, s, dist, queue, jump);
        const o = k * n;
        for (let j = 0; j < n; j++) out[o + j] += table[dist[j]];
      });
    }
  }
  for (let i = 0; i < out.length; i++) out[i] /= frames;
  return out;
}

/** Contact moyen entre deux groupes de monomères, depuis des lignes calculées pour `a`. */
export function mean(rowsOfA: Float32Array, n: number, nA: number, b: number[]): number {
  let s = 0;
  for (let k = 0; k < nA; k++) for (const j of b) s += rowsOfA[k * n + j];
  return s / (nA * b.length);
}
