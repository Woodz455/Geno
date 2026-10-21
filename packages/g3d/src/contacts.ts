/** La carte de contacts du modèle (semaine 14), lue tuile par tuile pour HiGlass.
 *
 *  Miroir de `pipeline/geno_pipeline/export/linked.py`. Deux jeux de tuiles — `noyau`
 *  (24 Mb → 750 kb, tout le génome) et `fin` (512 kb → 2 kb, les régions fines) —, une
 *  tuile par requête Range, vérifiée contre son empreinte, gardée ensuite. Une tuile
 *  `z.x.y` : colonnes = bins x, lignes = bins y ; seules les tuiles x ≤ y existent.
 */

import { CorruptError, hex, inflateRaw, sha256 } from "./format.js";
import type { Source } from "./source.js";

const MAGIC = [0x89, 0x47, 0x43, 0x54, 0x0d, 0x0a, 0x1a, 0x0a];
const PREAMBLE = 64;
export const TILE = 256;

interface BlobSpec {
  offset: number;
  length: number;
  sha256: string;
}

export interface ContactLevel {
  z: number;
  resolution: number;
  pixels: number;
  tiles: Record<string, BlobSpec>;
  total: number;
  max: number;
}

export interface Region {
  name: string;
  chrom: string;
  start: number;
  end: number;
  abs_start: number;
  abs_end: number;
}

export interface Tileset {
  uid: string;
  label: string;
  evidence: string;
  resolutions: number[];
  bins_per_dimension: number;
  min_pos: [number, number];
  max_pos: [number, number];
  value: string;
  sources: { note: string; structures: number; cutoff: number; bins: number; total: number }[];
  regions: Region[] | null;
  levels: ContactLevel[];
}

export interface ContactsHeader {
  format: string;
  version: number;
  assembly: string;
  cell_type: string;
  chromsizes: [string, number][];
  tilesets: Tileset[];
  provenance: Record<string, unknown>;
}

/** float16 → float32, par table : 65 536 entrées calculées une fois. */
let F16: Float32Array | null = null;
function f16table(): Float32Array {
  if (F16) return F16;
  F16 = new Float32Array(65536);
  for (let h = 0; h < 65536; h++) {
    const s = h & 0x8000 ? -1 : 1;
    const e = (h >> 10) & 0x1f;
    const m = h & 0x3ff;
    F16[h] = e === 0 ? s * m * 2 ** -24 : e === 31 ? (m ? NaN : s * Infinity) : s * (1 + m / 1024) * 2 ** (e - 15);
  }
  return F16;
}

export class ContactsReader {
  private readonly tiles = new Map<string, Float32Array | null>();
  private readonly pending = new Map<string, Promise<Float32Array | null>>();
  /** Octets lus et requêtes faites depuis l'ouverture. */
  bytes = 0;
  requests = 0;
  /** Décalage absolu de chaque chromosome, dans l'ordre de `chromsizes`. */
  readonly offsets = new Map<string, number>();
  readonly total: number;

  private constructor(
    private readonly source: Source,
    readonly header: ContactsHeader,
    private readonly base: number,
    private readonly verify: boolean,
  ) {
    let at = 0;
    for (const [c, n] of header.chromsizes) {
      this.offsets.set(c, at);
      at += n;
    }
    this.total = at;
  }

  static async open(source: Source, opts: { verify?: boolean } = {}): Promise<ContactsReader> {
    const pre = await source.read(0, PREAMBLE);
    if (!MAGIC.every((b, i) => pre[i] === b)) throw new Error("ce n'est pas un fichier de contacts Geno");
    const dv = new DataView(pre.buffer, pre.byteOffset, pre.byteLength);
    const version = dv.getUint32(8, true);
    if (version !== 1) throw new Error(`contacts version ${version}, ce lecteur lit la 1`);
    const hlen = dv.getUint32(12, true);
    const hb = await source.read(PREAMBLE, hlen);
    const verify = opts.verify ?? true;
    if (verify && hex(await sha256(hb)) !== hex(pre.subarray(16, 48))) throw new CorruptError("contacts : en-tête divergent");
    const header = JSON.parse(new TextDecoder().decode(await inflateRaw(hb))) as ContactsHeader;
    const r = new ContactsReader(source, header, PREAMBLE + hlen, verify);
    r.bytes = PREAMBLE + hlen;
    r.requests = 2;
    return r;
  }

  tileset(uid: string): Tileset {
    const t = this.header.tilesets.find((s) => s.uid === uid);
    if (!t) throw new Error(`pas de jeu de tuiles ${uid}`);
    return t;
  }

  /** chrom:pos → position absolue (celle de HiGlass) ; `null` hors du génome du modèle. */
  absolute(chrom: string, pos: number): number | null {
    const o = this.offsets.get(chrom);
    return o === undefined ? null : o + pos;
  }

  /** Position absolue → chrom:pos (0-based). */
  relative(abs: number): { chrom: string; pos: number } | null {
    let at = 0;
    for (const [c, n] of this.header.chromsizes) {
      if (abs < at + n) return { chrom: c, pos: Math.max(0, Math.floor(abs - at)) };
      at += n;
    }
    return null;
  }

  /** La tuile (x ≤ y), 256 × 256 en lignes ; `null` si elle est vide (non écrite). */
  async tile(uid: string, z: number, x: number, y: number): Promise<Float32Array | null> {
    const key = `${uid}.${z}.${x}.${y}`;
    if (this.tiles.has(key)) return this.tiles.get(key)!;
    let p = this.pending.get(key);
    if (!p) {
      const spec = this.tileset(uid).levels[z]?.tiles[`${z}.${x}.${y}`];
      p = (spec ? this.decode(spec) : Promise.resolve(null)).then((t) => {
        this.tiles.set(key, t);
        this.pending.delete(key);
        return t;
      });
      this.pending.set(key, p);
    }
    return p;
  }

  private async decode(spec: BlobSpec): Promise<Float32Array> {
    const d = await this.source.read(this.base + spec.offset, spec.length);
    this.bytes += d.length;
    this.requests++;
    if (this.verify && hex(await sha256(d)).slice(0, 32) !== spec.sha256) throw new CorruptError("contacts : tuile divergente");
    const raw = await inflateRaw(d);
    // Octets regroupés par rang : les poids faibles d'abord, puis les forts.
    const n = raw.length >> 1;
    const t = f16table();
    const out = new Float32Array(n);
    for (let i = 0; i < n; i++) out[i] = t[raw[i] | (raw[n + i] << 8)];
    return out;
  }

  /** La valeur d'un pixel, miroir compris (0 si la tuile n'existe pas). */
  async pixel(uid: string, z: number, row: number, col: number): Promise<number> {
    if (col > row) [row, col] = [col, row];
    const t = await this.tile(uid, z, Math.floor(col / TILE), Math.floor(row / TILE));
    return t ? t[(row % TILE) * TILE + (col % TILE)] : 0;
  }
}
