/** La séquence de référence d'un `.g3d` (semaine 12), lue par fenêtre.
 *
 *  À part de `G3D` : la page du premier rendu embarque cette classe-là tout entière, et la
 *  séquence n'y sert à rien. Ne lit **que** les blocs qui couvrent la fenêtre demandée, en
 *  une lecture pour des blocs contigus, et vérifie chacun contre sa table d'empreintes.
 */

import { CorruptError, hex, type SequenceSpec, sha256 } from "./format.js";
import type { G3D } from "./reader.js";

export class SequenceReader {
  /** Blocs décodés, et ceux en cours de lecture (clé : intervalle, bloc). */
  private readonly bases = new Map<string, string>();
  private readonly pending = new Map<string, Promise<void>>();
  private readonly hashTables = new Map<number, Promise<Uint8Array>>();
  /** Octets de séquence lus depuis l'ouverture : de quoi prouver qu'on ne lit que la fenêtre. */
  bytes = 0;

  constructor(private readonly g: G3D) {}

  /** Les intervalles couverts. */
  get intervals(): SequenceSpec[] {
    return this.g.header.sequences ?? [];
  }

  /** L'intervalle de séquence qui contient `[start, end)`, ou -1. */
  private sequenceIndex(chrom: string, start: number, end: number): number {
    return (this.g.header.sequences ?? []).findIndex((q) => q.chrom === chrom && q.start <= start && end <= q.end);
  }

  /** Les blocs à lire pour `[start, end)`, et ceux qui manquent encore. */
  private blocks(si: number, start: number, end: number): { spec: SequenceSpec; ks: number[] } {
    const spec = this.g.header.sequences![si];
    const k0 = Math.floor((start - spec.start) / spec.chunk_bp);
    const k1 = Math.floor((Math.max(start, end - 1) - spec.start) / spec.chunk_bp);
    const ks: number[] = [];
    for (let k = k0; k <= k1; k++) ks.push(k);
    return { spec, ks };
  }

  /** Vrai si la séquence de `[start, end)` est déjà là : une image peut la dessiner
   *  sans attendre. Synchrone, pour que le rendu sache s'il montre un chargement. */
  has(chrom: string, start: number, end: number): boolean {
    const si = this.sequenceIndex(chrom, start, end);
    if (si < 0) return false;
    return this.blocks(si, start, end).ks.every((k) => this.bases.has(`${si}:${k}`));
  }

  /** Lit les blocs de `[start, end)` qui manquent — **seulement eux**, en une lecture
   *  pour les blocs contigus — et vérifie chacun contre sa table d'empreintes. */
  async prefetch(chrom: string, start: number, end: number): Promise<void> {
    const si = this.sequenceIndex(chrom, start, end);
    if (si < 0) throw new Error(`pas de séquence pour ${chrom}:${start}-${end}`);
    const { spec, ks } = this.blocks(si, start, end);
    const wait: Promise<void>[] = [];
    const todo: number[] = [];
    for (const k of ks) {
      const key = `${si}:${k}`;
      if (this.bases.has(key)) continue;
      const p = this.pending.get(key);
      if (p) wait.push(p);
      else todo.push(k);
    }
    if (todo.length) {
      const read = (async () => {
        let table = this.hashTables.get(si);
        if (!table) {
          table = (async () => {
            const t = await this.g.readSpan(this.g.dataStart + spec.hashes.offset, spec.hashes.length);
            this.bytes += t.length;
            if (this.g.verifying && hex(await sha256(t)).slice(0, 32) !== spec.hashes.sha256) {
              throw new CorruptError("séquence : table d'empreintes divergente");
            }
            return t;
          })();
          this.hashTables.set(si, table);
        }
        const lo = todo[0], hi = todo[todo.length - 1];
        const from = spec.offset + lo * spec.chunk_bytes;
        const to = Math.min(spec.offset + (hi + 1) * spec.chunk_bytes, spec.offset + spec.length);
        const [t, bytes] = await Promise.all([table, this.g.readSpan(this.g.dataStart + from, to - from)]);
        this.bytes += bytes.length;
        for (let k = lo; k <= hi; k++) {
          const a = (k - lo) * spec.chunk_bytes;
          const block = bytes.subarray(a, Math.min(a + spec.chunk_bytes, bytes.length));
          if (this.g.verifying && hex(await sha256(block)).slice(0, 32) !== hex(t.subarray(16 * k, 16 * k + 16))) {
            throw new CorruptError(`séquence : bloc ${k} divergent`);
          }
          this.bases.set(`${si}:${k}`, unpack(block, spec, k));
        }
      })();
      for (const k of todo) this.pending.set(`${si}:${k}`, read);
      read.finally(() => todo.forEach((k) => this.pending.delete(`${si}:${k}`))).catch(() => {});
      wait.push(read);
    }
    await Promise.all(wait);
  }

  /** Les bases de `[start, end)` (0-based, demi-ouvert), en majuscules, N dans les trous. */
  async read(chrom: string, start: number, end: number): Promise<string> {
    await this.prefetch(chrom, start, end);
    return this.now(chrom, start, end)!;      // tous les blocs viennent d'être lus
  }

  /** Comme `sequence`, sans attendre : `null` si un bloc manque. */
  now(chrom: string, start: number, end: number): string | null {
    if (!this.has(chrom, start, end)) return null;
    const si = this.sequenceIndex(chrom, start, end);
    const { spec, ks } = this.blocks(si, start, end);
    const joined = ks.map((k) => this.bases.get(`${si}:${k}`)!).join("");
    const base = spec.start + ks[0] * spec.chunk_bp;
    return joined.slice(start - base, end - base);
  }
}

const LETTERS = "ACGT";

/** Un bloc de séquence décodé, avec ses N remis en place. */
function unpack(block: Uint8Array, spec: SequenceSpec, k: number): string {
  const first = spec.start + k * spec.chunk_bp;
  const n = Math.min(spec.chunk_bp, spec.end - first);
  const out = new Array<string>(n);
  for (let i = 0; i < n; i++) out[i] = LETTERS[(block[i >> 2] >> (6 - 2 * (i & 3))) & 3];
  for (const [a, b] of spec.n_runs) {
    for (let p = Math.max(a, first); p < Math.min(b, first + n); p++) out[p - first] = "N";
  }
  return out.join("");
}

