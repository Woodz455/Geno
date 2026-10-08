/** Le format `.g3d` côté lecture : préambule, en-tête, colonnes.
 *
 *  Miroir exact de `pipeline/geno_pipeline/export/g3d.py`. Rien ici ne dépend d'une
 *  bibliothèque : `DecompressionStream("deflate-raw")` et `crypto.subtle` sont dans
 *  tous les navigateurs récents et dans Node 22, et c'est tout ce qu'il faut.
 */

export const MAGIC = [0x89, 0x47, 0x33, 0x44, 0x0d, 0x0a, 0x1a, 0x0a];
export const VERSION = 1;
export const PREAMBLE = 64;

export interface ColumnSpec {
  offset: number;
  length: number;
  raw: number;
  dtype: "u1" | "u2" | "u4";
  planes: number;
  filters: string[];
  sha256: string;
}

export interface ChunkSpec {
  first: number;
  n: number;
  bbox: number[];
  /** Centre x, y, z puis rayon, en nm dans le repère du niveau. */
  sphere: [number, number, number, number];
  columns: Record<string, ColumnSpec>;
}

export interface LevelSpec {
  name: string;
  bp_per_bead: number;
  n_beads: number;
  evidence: "measured" | "simulated" | "deterministic" | "predicted";
  copies: string[];
  radius_nm: number[];
  frame: { units: "nm"; transform: number[]; fit: Record<string, unknown> | null };
  quant: { origin: number[]; step: number[]; max_error_nm: number; rule: string };
  variability: { kind: string; unit: "nm"; step: number; [k: string]: unknown };
  index: { copy_offsets: number[]; columns: Record<string, ColumnSpec> };
  notes: string[];
  chunks: ChunkSpec[];
}

export interface Header {
  format: "geno-g3d";
  version: number;
  first: string;
  levels: LevelSpec[];
  assembly?: string;
  cell_type?: string;
  karyotype?: string;
  provenance?: { pipeline: string; sources: Record<string, unknown>[] };
  warnings?: string[];
}

export interface Preamble {
  headerLength: number;
  firstFrameEnd: number;
  headerSha256: Uint8Array;
}

export class CorruptError extends Error {}

export function parsePreamble(bytes: Uint8Array): Preamble {
  if (bytes.length < PREAMBLE) throw new Error("préfixe plus court que le préambule");
  for (let i = 0; i < MAGIC.length; i++) {
    if (bytes[i] !== MAGIC[i]) throw new Error("ce n'est pas un .g3d (signature)");
  }
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const version = dv.getUint32(8, true);
  if (version !== VERSION) throw new Error(`.g3d version ${version}, ce lecteur lit la ${VERSION}`);
  return {
    headerLength: dv.getUint32(12, true),
    firstFrameEnd: dv.getUint32(16, true),
    headerSha256: bytes.slice(24, 56),
  };
}

export async function inflateRaw(data: Uint8Array): Promise<Uint8Array> {
  const stream = new Blob([data as BlobPart]).stream().pipeThrough(new DecompressionStream("deflate-raw"));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

async function sha256(data: Uint8Array): Promise<Uint8Array> {
  return new Uint8Array(await crypto.subtle.digest("SHA-256", data as BufferSource));
}

const hex = (b: Uint8Array) => Array.from(b, (x) => x.toString(16).padStart(2, "0")).join("");

export async function parseHeader(bytes: Uint8Array, pre: Preamble, verify = true): Promise<Header> {
  if (verify) {
    const got = await sha256(bytes);
    if (hex(got) !== hex(pre.headerSha256)) throw new CorruptError("en-tête : empreinte divergente");
  }
  return JSON.parse(new TextDecoder().decode(await inflateRaw(bytes))) as Header;
}

export type Column = Uint8Array | Uint16Array | Uint32Array;

function unshuffle(raw: Uint8Array, itemsize: number): Uint8Array {
  if (itemsize === 1) return raw;
  const n = raw.length / itemsize;
  const out = new Uint8Array(raw.length);
  for (let b = 0; b < itemsize; b++) {
    const plane = b * n;
    for (let i = 0; i < n; i++) out[i * itemsize + b] = raw[plane + i];
  }
  return out;
}

/** Somme cumulée par plan. L'affectation dans un tableau typé non signé replie modulo
 *  2^16 ou 2^32 : c'est exactement l'inverse du delta modulaire de l'écriture. */
function undelta(a: Column, planes: number): void {
  const n = a.length / planes;
  for (let p = 0; p < planes; p++) {
    const base = p * n;
    for (let i = 1; i < n; i++) a[base + i] = a[base + i] + a[base + i - 1];
  }
}

export async function decodeColumn(data: Uint8Array, spec: ColumnSpec, verify = true): Promise<Column> {
  if (verify) {
    const got = hex(await sha256(data)).slice(0, 32);
    if (got !== spec.sha256) throw new CorruptError("colonne : empreinte divergente");
  }
  const size = { u1: 1, u2: 2, u4: 4 }[spec.dtype];
  let raw = await inflateRaw(data);
  if (spec.filters.includes("shuffle")) raw = unshuffle(raw, size);
  // Copie dans un tampon neuf : l'alignement d'un Uint32Array sur un sous-tableau n'est pas garanti.
  const buf = raw.slice().buffer;
  const col: Column =
    size === 1 ? new Uint8Array(buf) : size === 2 ? new Uint16Array(buf) : new Uint32Array(buf);
  if (spec.filters.includes("delta")) undelta(col, spec.planes);
  return col;
}
