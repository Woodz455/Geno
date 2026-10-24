/** Le fichier de fiches (semaine 13) : ce que l'on sait de chaque bille, lu par bloc.
 *
 *  À part du `.g3d`, pour que le premier rendu n'en paie rien (`export/fiche` côté
 *  pipeline). Une fiche = l'enregistrement de son intervalle génomique (partagé par les
 *  homologues) + les champs du modèle propres à sa copie. Les enregistrements sont rangés
 *  par chromosome en blocs de 64 : un bloc par requête Range, vérifié contre son empreinte,
 *  gardé ensuite.
 *
 *  `card()` est **synchrone** : un clic dont le bloc est déjà là compose sa fiche sans
 *  attendre le réseau. C'est ce que le préchargement achète, et ce que la mesure compte.
 */

import { CorruptError, hex, inflateRaw, sha256 } from "./format.js";
import type { Source } from "./source.js";

const MAGIC = [0x89, 0x46, 0x43, 0x48, 0x0d, 0x0a, 0x1a, 0x0a];
const PREAMBLE = 64;

export interface FicheField {
  /** Évidence : mesure, autre, simule, reference, annotation. */
  e: string;
  /** Source : identifiant du manifeste (ou `g3d:…`, `ensemble:…` pour le modèle). */
  s: string;
  v?: any;
  absent?: string;
  [extra: string]: any;
}
export type Card = Record<string, FicheField>;

interface BlobSpec {
  offset: number;
  length: number;
  sha256: string;
}

export interface FicheLevel {
  name: string;
  n_beads: number;
  records: number;
  fine: boolean;
  record_of: BlobSpec & { dtype: string };
  model: BlobSpec;
  blocks: (BlobSpec & { first: number; n: number; chrom: string })[];
}

export interface FicheHeader {
  format: string;
  version: number;
  assembly: string;
  cell_type: string;
  g3d: { path: string; sha256: string };
  pipeline: string;
  evidence: Record<string, string>;
  sections: Record<string, string[]>;
  fine_only: Record<string, string[]>;
  labels: Record<string, string>;
  sources: Record<string, { url: string; sha256: string; license: string; notes: string }>;
  tissues: string[];
  clinvar_classes: Record<string, string>;
  repli_terciles: number[];
  nuclear_radius_nm: number;
  /** Par chaîne : chromosome cible, sa taille, brin, score (le plus haut gagne, comme `liftOver`). */
  chains: { hs1: { targets: [string, number, string, number][]; regions: { chrom: string; start: number; end: number; blocks: [number, number, number, number][] }[] } };
  levels: FicheLevel[];
  /** Ce que les fiches partagent : gènes (un bloc par chromosome), maladies, traits GWAS. */
  tables: {
    genes: (BlobSpec & { chrom: string; first: number; n: number })[];
    diseases: (BlobSpec & { first: number; n: number })[];
    traits: (BlobSpec & { first: number; n: number })[];
  };
}

export interface Gene {
  id: string;
  name: string;
  type: string;
  strand: string;
  start: number;
  end: number;
  mane: boolean;
  gtex: [number, number, number] | null;
  gnomad: [number, number] | null;
}

type TableName = "genes" | "diseases" | "traits";
const GENE_KEYS = ["id", "name", "type", "strand", "start", "end", "mane"] as const;

export interface Lifted {
  chrom: string;
  pos: number;
  strand: "+" | "-";
}

export class FicheReader {
  private readonly blocks = new Map<string, any[]>();
  private readonly pending = new Map<string, Promise<any>>();
  private readonly recordOf = new Map<string, Uint32Array>();
  private readonly models = new Map<string, Card[]>();
  private readonly rows = new Map<string, any[]>();       // `table:bloc` → lignes décodées
  /** Octets lus et requêtes faites depuis l'ouverture. */
  bytes = 0;
  requests = 0;
  private queue: Promise<void> = Promise.resolve();

  private constructor(
    private readonly source: Source,
    readonly header: FicheHeader,
    private readonly base: number,
    private readonly verify: boolean,
  ) {}

  /** Préambule et en-tête : deux lectures, la seconde de la taille exacte de l'en-tête. */
  static async open(source: Source, opts: { verify?: boolean } = {}): Promise<FicheReader> {
    const pre = await source.read(0, PREAMBLE);
    if (!MAGIC.every((b, i) => pre[i] === b)) throw new Error("ce n'est pas un fichier de fiches");
    const dv = new DataView(pre.buffer, pre.byteOffset, pre.byteLength);
    const version = dv.getUint32(8, true);
    if (version !== 1) throw new Error(`fiche version ${version}, ce lecteur lit la 1`);
    const hlen = dv.getUint32(12, true);
    const hsha = pre.subarray(16, 48);
    const hb = await source.read(PREAMBLE, hlen);
    const verify = opts.verify ?? true;
    if (verify && hex(await sha256(hb)) !== hex(hsha)) throw new CorruptError("fiche : en-tête divergent");
    const header = JSON.parse(new TextDecoder().decode(await inflateRaw(hb))) as FicheHeader;
    const r = new FicheReader(source, header, PREAMBLE + hlen, verify);
    r.bytes = PREAMBLE + hlen;
    r.requests = 2;
    return r;
  }

  level(name: string): FicheLevel | undefined {
    return this.header.levels.find((l) => l.name === name);
  }

  private async blob(spec: BlobSpec): Promise<Uint8Array> {
    const d = await this.source.read(this.base + spec.offset, spec.length);
    this.bytes += d.length;
    this.requests++;
    if (this.verify && hex(await sha256(d)).slice(0, 32) !== spec.sha256) throw new CorruptError("fiche : bloc divergent");
    return inflateRaw(d);
  }

  private async json(spec: BlobSpec): Promise<any> {
    return JSON.parse(new TextDecoder().decode(await this.blob(spec)));
  }

  /** Index et colonne du modèle d'un niveau : petits, lus une fois. */
  async openLevel(name: string): Promise<void> {
    const lv = this.level(name);
    if (!lv || this.models.has(name)) return;
    const [ro, model] = await Promise.all([this.blob(lv.record_of), this.json(lv.model)]);
    this.recordOf.set(name, new Uint32Array(ro.buffer.slice(ro.byteOffset, ro.byteOffset + ro.byteLength)));
    this.models.set(name, model as Card[]);
  }

  private blockOf(name: string, beadId: number): number {
    const lv = this.level(name)!;
    const r = this.recordOf.get(name)![beadId - 1];
    return Math.floor(r / (lv.blocks[0]?.n || 64));
  }

  /** Le bloc `b` du niveau, lu au plus une fois même demandé deux fois en même temps. */
  private async block(name: string, b: number): Promise<any[]> {
    const key = `${name}:${b}`;
    const hit = this.blocks.get(key);
    if (hit) return hit;
    let p = this.pending.get(key);
    if (!p) {
      p = this.json(this.level(name)!.blocks[b]).then((recs: any[]) => {
        this.blocks.set(key, recs);
        this.pending.delete(key);
        return recs;
      });
      this.pending.set(key, p);
    }
    return p;
  }

  /** Le bloc d'une table qui contient la ligne `i`. */
  private tableBlock(t: TableName, i: number): number {
    const specs = this.header.tables[t];
    let lo = 0;
    let hi = specs.length - 1;
    while (lo < hi) {
      const m = (lo + hi + 1) >> 1;
      if (specs[m].first <= i) lo = m;
      else hi = m - 1;
    }
    return lo;
  }

  private row(t: TableName, i: number): any {
    const b = this.tableBlock(t, i);
    return this.rows.get(`${t}:${b}`)![i - this.header.tables[t][b].first];
  }

  private async tableRows(t: TableName, b: number): Promise<any[]> {
    const key = `${t}:${b}`;
    const hit = this.rows.get(key);
    if (hit) return hit;
    let p = this.pending.get(key);
    if (!p) {
      const spec = this.header.tables[t][b];
      p = this.json(spec).then((part: any) => {
        const rows = t === "genes"
          ? Array.from({ length: spec.n }, (_, i) => ({
              id: part.id[i], name: part.name[i], type: part.type[i], strand: part.strand[i],
              start: part.start[i], end: part.end[i], mane: part.mane[i] === "1", gtex: part.gtex[i], gnomad: part.gnomad[i],
            }))
          : part;
        this.rows.set(key, rows);
        this.pending.delete(key);
        return rows;
      });
      this.pending.set(key, p);
    }
    return p;
  }

  /** Les blocs de tables qu'un enregistrement référence. */
  private needs(rec: any): [TableName, number][] {
    const out = new Set<string>();
    const gi: number[] = [...rec.genes.v.genes, ...(rec.genes.v.nearest ?? []).map((n: any) => n.g)];
    for (const i of gi) out.add(`genes:${this.tableBlock("genes", i)}`);
    const cv = rec.clinvar.v;
    for (const x of [...(cv.variants ?? []), ...(cv.plp ?? [])]) for (const d of x.diseases) out.add(`diseases:${this.tableBlock("diseases", d)}`);
    for (const h of rec.gwas.v.top) out.add(`traits:${this.tableBlock("traits", h.t)}`);
    return [...out].map((k) => {
      const [t, b] = k.split(":");
      return [t as TableName, Number(b)];
    });
  }

  private record(name: string, beadId: number): any | null {
    if (!this.recordOf.has(name)) return null;
    const lv = this.level(name)!;
    const r = this.recordOf.get(name)![beadId - 1];
    const per = lv.blocks[0].n;
    return this.blocks.get(`${name}:${Math.floor(r / per)}`)?.[r % per] ?? null;
  }

  /** Vrai si la fiche est composable sans attendre : son bloc, et les blocs de tables
   *  qu'elle référence. */
  has(name: string, beadId: number): boolean {
    const rec = this.record(name, beadId);
    return !!rec && this.needs(rec).every(([t, b]) => this.rows.has(`${t}:${b}`));
  }

  /** La fiche, sans attendre : `null` si un bloc manque. */
  card(name: string, beadId: number): Card | null {
    if (!this.has(name, beadId)) return null;
    return { ...this.expand(this.record(name, beadId)), ...this.models.get(name)![beadId - 1] };
  }

  /** La fiche, en lisant ce qui manque. */
  async cardAsync(name: string, beadId: number): Promise<Card> {
    await this.openLevel(name);
    await this.block(name, this.blockOf(name, beadId));
    await Promise.all(this.needs(this.record(name, beadId)).map(([t, b]) => this.tableRows(t, b)));
    return this.card(name, beadId)!;
  }

  /** Les indices remplacés par ce qu'ils désignent — `fiche.build.expand` côté pipeline. */
  private expand(rec: any): Card {
    const gene = (i: number) => {
      const g = this.row("genes", i) as Gene;
      return Object.fromEntries(GENE_KEYS.map((k) => [k, g[k]]));
    };
    const r: any = { ...rec };
    const gv = rec.genes.v;
    r.genes = { ...rec.genes, v: { genes: gv.genes.map(gene) } };
    if (gv.nearest) r.genes.v.nearest = gv.nearest.map((n: any) => ({ ...gene(n.g), distance: n.distance, side: n.side }));
    if (rec.expression.v) {
      r.expression = { ...rec.expression, v: Object.fromEntries(rec.expression.v.map((i: number) => {
        const g = this.row("genes", i) as Gene;
        return [g.id, { lcl: g.gtex![0], tissue: g.gtex![1], max: g.gtex![2] }];
      })) };
    }
    if (rec.contrainte.v) {
      r.contrainte = { ...rec.contrainte, v: Object.fromEntries(rec.contrainte.v.map((i: number) => {
        const g = this.row("genes", i) as Gene;
        return [g.id, { pli: g.gnomad![0], loeuf: g.gnomad![1] }];
      })) };
    }
    const cv = { ...rec.clinvar.v };
    const variant = (x: any) => ({ ...x, diseases: x.diseases.map((d: number) => this.row("diseases", d)) });
    if (cv.variants) {
      cv.variants = cv.variants.map(variant);
      cv.plp = cv.variants.filter((x: any) => x.sig === "plp").sort((a: any, b: any) => b.stars - a.stars || a.pos - b.pos);
    } else {
      cv.plp = cv.plp.map(variant);
    }
    r.clinvar = { ...rec.clinvar, v: cv };
    r.gwas = { ...rec.gwas, v: { ...rec.gwas.v, top: rec.gwas.v.top.map((h: any) => {
      const [trait, mapped, uri] = this.row("traits", h.t);
      const { t: _t, ...rest } = h;
      return { ...rest, trait, mapped, uri };
    }) } };
    return r as Card;
  }

  /** Lit en arrière-plan, **une requête à la fois**, les blocs des chromosomes demandés
   *  d'abord puis tous les autres, pour chaque niveau. Rend une promesse résolue quand
   *  tout est là ; les clics n'attendent pas cette promesse, ils lisent `has`. */
  prefetch(order: { level: string; first?: string[] }[], progress?: (step: string) => void): Promise<void> {
    this.queue = this.queue.then(async () => {
      // Les tables d'abord : sans elles, aucune fiche n'est composable.
      const firstChroms = order[0]?.first ?? [];
      const g = this.header.tables.genes.map((b, k) => ({ k, r: firstChroms.includes(b.chrom) ? firstChroms.indexOf(b.chrom) : firstChroms.length }));
      for (const { k } of g.sort((a, b) => a.r - b.r || a.k - b.k)) await this.tableRows("genes", k);
      for (const t of ["diseases", "traits"] as const) {
        for (let k = 0; k < this.header.tables[t].length; k++) await this.tableRows(t, k);
      }
      progress?.("tables");
      for (const { level, first = [] } of order) {
        await this.openLevel(level);
        const lv = this.level(level);
        if (!lv) continue;
        const rank = (c: string) => (first.includes(c) ? first.indexOf(c) : first.length);
        const ks = lv.blocks.map((b, k) => ({ k, r: rank(b.chrom) })).sort((a, b) => a.r - b.r || a.k - b.k);
        for (const { k } of ks) await this.block(level, k);
        progress?.(level);
      }
    });
    return this.queue;
  }

  /** Taille des tables partagées, octets. */
  tablesBytes(): number {
    const t = this.header.tables;
    return [...t.genes, ...t.diseases, ...t.traits].reduce((a, b) => a + b.length, 0);
  }

  /** Taille totale des blocs d'un niveau, octets. */
  levelBytes(name: string): number {
    const lv = this.level(name);
    return lv ? lv.blocks.reduce((a, b) => a + b.length, 0) + lv.record_of.length + lv.model.length : 0;
  }

  /** Une position hg38 (0-based) en T2T-CHM13, par les blocs de chaîne embarqués pour les
   *  régions fines. `null` hors de ces régions ou hors de tout bloc — pas d'interpolation. */
  liftHs1(chrom: string, pos: number): Lifted | null {
    const c = this.header.chains.hs1;
    const reg = c.regions.find((r) => r.chrom === chrom && r.start <= pos && pos < r.end);
    if (!reg) return null;
    let best: { k: number; q: number } | null = null;
    for (const [t0, t1, q0, k] of reg.blocks) {
      if (pos >= t0 && pos < t1 && (best === null || c.targets[k][3] > c.targets[best.k][3])) best = { k, q: q0 + (pos - t0) };
    }
    if (!best) return null;
    const [qName, qSize, qStrand] = c.targets[best.k];
    return { chrom: qName, pos: qStrand === "-" ? qSize - 1 - best.q : best.q, strand: qStrand as "+" | "-" };
  }
}
