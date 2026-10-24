/** L'index de recherche (semaine 16) : un gène, une position, un rsID → un intervalle hg38.
 *
 *  Miroir de `pipeline/geno_pipeline/search/` — même grammaire (`query.py`), même fichier
 *  (`build.py`), même réponse (`resolve.py`). Deux tables triées par clé, en blocs : les gènes
 *  (nom en majuscules et ENSG) et les rsID de ClinVar et du catalogue GWAS — **pas** dbSNP. Une
 *  recherche lit l'en-tête une fois, puis un bloc par requête Range, vérifié et gardé.
 */

import { CorruptError, hex, inflateRaw, sha256 } from "./format.js";
import type { Source } from "./source.js";

const MAGIC = [0x89, 0x47, 0x53, 0x52, 0x0d, 0x0a, 0x1a, 0x0a];
const PREAMBLE = 64;

interface Block {
  first: string | number;
  last?: number;
  n: number;
  offset: number;
  length: number;
  sha256: string;
}

export interface SearchHeader {
  format: string;
  version: number;
  assembly: string;
  chroms: [string, number][];
  model_chroms: string[];
  tables: { genes: { blocks: Block[]; rows: number }; rs: { blocks: Block[]; rows: number; distinct: number } };
  scope: string;
  provenance: Record<string, unknown>;
}

export interface Query {
  kind: "rs" | "position" | "interval" | "gene-id" | "gene" | "malformed";
  key: string;
  chrom: string;
  start: number;
  end: number;
  error: string;
}

export interface Hit {
  label: string;
  chrom: string;
  /** 0-based, demi-ouvert. */
  start: number;
  end: number;
  in_model: boolean;
  id?: string;
  strand?: string;
  type?: string;
  sources?: string[];
}

export interface Resolution {
  kind: Query["kind"];
  hits: Hit[];
  /** "" si au moins une cible ; sinon malformed, unknown-gene, unknown-rs, unknown-chrom, out-of-range. */
  error: string;
}

const RS = /^rs([0-9]{1,15})$/i;
const POS = /^(?:chr)?([0-9]{1,2}|x|y|m|mt)\s*[:\s]\s*([0-9,_   ]+?)(?:\s*[-–]\s*([0-9,_   ]+))?$/i;
const ENSG = /^(ENSG[0-9]{11})(?:\.[0-9]+)?$/i;
const NAME = /^[A-Za-z0-9][A-Za-z0-9._@-]{0,39}$/;
const DIGITS = /[,_   ]/g;

const chromName = (c: string): string => {
  const u = c.toUpperCase();
  return u === "M" || u === "MT" ? "chrM" : `chr${u}`;
};

/** La requête telle qu'on la tape — 1-based, bornes incluses — en intervalle 0-based demi-ouvert. */
export function parseQuery(text: string): Query {
  const q = (kind: Query["kind"], o: Partial<Query> = {}): Query => ({ kind, key: "", chrom: "", start: 0, end: 0, error: "", ...o });
  const t = text.trim();
  if (!t) return q("malformed", { error: "requête vide" });
  let m = RS.exec(t);
  if (m) return q("rs", { key: String(Number(m[1])) });
  m = POS.exec(t);
  if (m) {
    const chrom = chromName(m[1]);
    const a = Number(m[2].replace(DIGITS, ""));
    if (m[3] === undefined) {
      if (a < 1) return q("malformed", { error: "les coordonnées 1-based commencent à 1" });
      return q("position", { chrom, start: a - 1, end: a });
    }
    const b = Number(m[3].replace(DIGITS, ""));
    if (a < 1) return q("malformed", { error: "les coordonnées 1-based commencent à 1" });
    if (b < a) return q("malformed", { error: "fin avant début" });
    return q("interval", { chrom, start: a - 1, end: b });
  }
  m = ENSG.exec(t);
  if (m) return q("gene-id", { key: m[1].toUpperCase() });
  if (NAME.test(t) && !/^[0-9]+$/.test(t)) return q("gene", { key: t.toUpperCase() });
  return q("malformed", { error: "ni un gène, ni une position chrom:pos, ni un rsID" });
}

const fmt = (n: number) => n.toLocaleString("en-US");

export class SearchIndex {
  private readonly blocks = new Map<string, Promise<unknown>>();
  readonly sizes: Map<string, number>;
  readonly model: Set<string>;
  bytes = 0;
  requests = 0;

  private constructor(
    private readonly source: Source,
    readonly header: SearchHeader,
    private readonly base: number,
    private readonly verify: boolean,
  ) {
    this.sizes = new Map(header.chroms);
    this.model = new Set(header.model_chroms);
  }

  static async open(source: Source, opts: { verify?: boolean } = {}): Promise<SearchIndex> {
    const pre = await source.read(0, PREAMBLE);
    if (!MAGIC.every((b, i) => pre[i] === b)) throw new Error("ce n'est pas un index de recherche Geno");
    const dv = new DataView(pre.buffer, pre.byteOffset, pre.byteLength);
    const version = dv.getUint32(8, true);
    if (version !== 1) throw new Error(`index version ${version}, ce lecteur lit la 1`);
    const hlen = dv.getUint32(12, true);
    const hb = await source.read(PREAMBLE, hlen);
    const verify = opts.verify ?? true;
    if (verify && hex(await sha256(hb)) !== hex(pre.subarray(16, 48))) throw new CorruptError("recherche : en-tête divergent");
    const header = JSON.parse(new TextDecoder().decode(await inflateRaw(hb))) as SearchHeader;
    const r = new SearchIndex(source, header, PREAMBLE + hlen, verify);
    r.bytes = PREAMBLE + hlen;
    r.requests = 2;
    return r;
  }

  private async raw(spec: Block): Promise<Uint8Array> {
    const d = await this.source.read(this.base + spec.offset, spec.length);
    this.bytes += d.length;
    this.requests++;
    if (this.verify && hex(await sha256(d)).slice(0, 32) !== spec.sha256) throw new CorruptError("recherche : bloc divergent");
    return inflateRaw(d);
  }

  private cached<T>(key: string, make: () => Promise<T>): Promise<T> {
    let p = this.blocks.get(key) as Promise<T> | undefined;
    if (!p) {
      p = make();
      this.blocks.set(key, p);
      p.catch(() => this.blocks.delete(key));
    }
    return p;
  }

  /** Les lignes d'un nom (en majuscules) ou d'un ENSG : [clé, nom, id, chrom, début, fin, brin, biotype]. */
  async genes(key: string): Promise<(string | number)[][]> {
    const bl = this.header.tables.genes.blocks;
    // Le dernier bloc qui commence avant la clé peut finir par elle ; les suivants qui
    // commencent par elle la continuent (même règle que Python).
    let lo = 0, hi = bl.length;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if ((bl[mid].first as string) < key) lo = mid + 1;
      else hi = mid;
    }
    const out: (string | number)[][] = [];
    for (let k = Math.max(0, lo - 1); k < bl.length && (bl[k].first as string) <= key; k++) {
      const rows = await this.cached(`g${k}`, async () => JSON.parse(new TextDecoder().decode(await this.raw(bl[k]))) as (string | number)[][]);
      for (const r of rows) if (r[0] === key) out.push(r);
    }
    return out;
  }

  /** Les positions d'un rsID : [chrom, pos 0-based, sources (bit 1 ClinVar, bit 2 GWAS)]. */
  async rs(n: number): Promise<[string, number, number][]> {
    const bl = this.header.tables.rs.blocks;
    let lo = 0, hi = bl.length;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if ((bl[mid].first as number) <= n) lo = mid + 1;
      else hi = mid;
    }
    const k = lo - 1;
    if (k < 0 || n > (bl[k].last as number)) return [];
    const b = bl[k];
    const cols = await this.cached(`r${k}`, async () => {
      const raw = await this.raw(b);
      const dv = new DataView(raw.buffer, raw.byteOffset, raw.byteLength);
      const N = b.n;
      const rs = new Float64Array(N), pos = new Uint32Array(N);
      for (let i = 0; i < N; i++) {
        rs[i] = (b.first as number) + dv.getUint32(4 * i, true);
        pos[i] = dv.getUint32(5 * N + 4 * i, true);
      }
      return { rs, chrom: raw.slice(4 * N, 5 * N), pos, src: raw.slice(9 * N, 10 * N) };
    });
    const out: [string, number, number][] = [];
    for (let i = 0; i < b.n; i++) if (cols.rs[i] === n) out.push([this.header.chroms[cols.chrom[i]][0], cols.pos[i], cols.src[i]]);
    return out;
  }

  /** La même réponse que `resolve()` en Python. */
  async resolve(text: string): Promise<Resolution> {
    const q = parseQuery(text);
    const hits: Hit[] = [];
    let error = "";
    if (q.kind === "malformed") return { kind: q.kind, hits, error: "malformed" };
    if (q.kind === "gene" || q.kind === "gene-id") {
      for (const r of await this.genes(q.key)) {
        const [, name, id, chrom, s, e, strand, type] = r as [string, string, string, string, number, number, string, string];
        hits.push({ label: name, id, chrom, start: s, end: e, strand, type, in_model: false });
      }
      error = hits.length ? "" : "unknown-gene";
    } else if (q.kind === "rs") {
      for (const [chrom, pos, src] of await this.rs(Number(q.key))) {
        const sources = [["clinvar", 1], ["gwas", 2]].filter(([, b]) => src & (b as number)).map(([s]) => s as string);
        hits.push({ label: `rs${q.key}`, chrom, start: pos, end: pos + 1, sources, in_model: false });
      }
      error = hits.length ? "" : "unknown-rs";
    } else {
      const size = this.sizes.get(q.chrom);
      if (size === undefined) return { kind: q.kind, hits, error: "unknown-chrom" };
      if (q.end > size) return { kind: q.kind, hits, error: "out-of-range" };
      const label = q.kind === "position" ? `${q.chrom}:${fmt(q.start + 1)}` : `${q.chrom}:${fmt(q.start + 1)}-${fmt(q.end)}`;
      hits.push({ label, chrom: q.chrom, start: q.start, end: q.end, in_model: false });
    }
    for (const h of hits) h.in_model = this.model.has(h.chrom);
    return { kind: q.kind, hits, error };
  }
}
