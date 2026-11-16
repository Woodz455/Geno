/** Début, fin et dispersion de **toutes** les billes d'un niveau, lus une fois (semaine 13).
 *
 *  Pour qu'un clic ne coûte plus de requête : `G3D.identify` relit les colonnes de la bille
 *  à chaque appel. À part de `G3D`, comme la séquence : la page du premier rendu embarque
 *  cette classe-là tout entière, et la mesure S9 a vu chaque kilo-octet de plus.
 */

import { decodeColumn } from "./format.js";
import type { G3D } from "./reader.js";

/** Début et fin (0-based, demi-ouvert) et dispersion de chaque bille, ordre des identifiants. */
export interface Intervals {
  start: Uint32Array;
  end: Uint32Array;
  variabilityNm: Float64Array;
}

export class BeadIndex {
  private readonly spans = new Map<string, Intervals>();

  constructor(private readonly g: G3D) {}

  /** Une lecture couvrant les colonnes de tous les chunks, la première fois ; la mémoire ensuite. */
  async load(name: string): Promise<Intervals> {
    const hit = this.spans.get(name);
    if (hit) return hit;
    const lv = this.g.level(name);
    const specs = lv.chunks.flatMap((c) => [c.columns.start, c.columns.end, c.columns.variability]);
    const lo = Math.min(...specs.map((s) => s.offset));
    const hi = Math.max(...specs.map((s) => s.offset + s.length));
    const bytes = await this.g.readSpan(this.g.dataStart + lo, hi - lo);
    const cols = await Promise.all(specs.map((s) => decodeColumn(bytes.subarray(s.offset - lo, s.offset - lo + s.length), s, this.g.verifying)));
    const start = new Uint32Array(lv.n_beads);
    const end = new Uint32Array(lv.n_beads);
    const variabilityNm = new Float64Array(lv.n_beads);
    lv.chunks.forEach((c, k) => {
      start.set(cols[3 * k] as Uint32Array, c.first);
      end.set(cols[3 * k + 1] as Uint32Array, c.first);
      const v = cols[3 * k + 2];
      for (let i = 0; i < c.n; i++) variabilityNm[c.first + i] = v[i] * lv.variability.step;
    });
    const out = { start, end, variabilityNm };
    this.spans.set(name, out);
    return out;
  }

  /** Comme `load`, sans attendre : `null` si le niveau n'a pas encore été lu. */
  now(name: string): Intervals | null {
    return this.spans.get(name) ?? null;
  }
}
