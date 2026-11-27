/** Le locus EPHA4 côté navigateur : construire un allèle, prédire ses contacts.
 *
 *  Miroir de `pipeline/geno_pipeline/sv/locus.py` (`Locus.sites`, `landmarks`) : un monomère
 *  dont le milieu tombe dans la délétion disparaît, les sites CTCF qu'il portait aussi, sauf
 *  ceux qu'un témoin « épargne », posés au point de recollement.
 */

import { type Law, type Params, rows, simulate, type SiteSet } from "./extrusion.js";

export interface Gene { name: string; start: number; end: number; strand: string; tss: number }
export interface Peak { start: number; end: number; motif: "+" | "-" | null; signal: number; occupancy: number | null }
export interface Deletion { name?: string; kind?: string; start: number; end: number; keep?: [number, number][] }
export interface AlleleDoc {
  name: string; kind: string | null; variant: Deletion | null; frame: number; frames: number; n: number;
  sigma_nm: number; coords: string; mids_first: number;
}

export interface VariantDoc {
  format: string;
  region: { chrom: string; start: number; end: number; bp_per_bead: number; n: number };
  genes: Gene[];
  ctcf: Peak[];
  enhancer: [number, number];
  boundary: [number, number];
  tad: [number, number];
  promoter_half: number;
  variants: Deletion[];
  params: Params;
  law: Law & { r_log: number; pairs: number };
  alleles: Record<string, AlleleDoc>;
  clinical: { clinvar_plp: number[]; gwas: number[]; genes: { name: string; start: number; end: number; strand: string; pLI: number | null; LOEUF: number | null }[] };
  evidence: Record<string, string>;
  source: string;
}

export interface Allele {
  n: number;
  /** Milieu génomique (hg38) de chaque monomère de l'allèle. */
  mids: Float64Array;
  /** Index du monomère dans le type sauvage. */
  wild: Int32Array;
  sites: SiteSet;
  deletion: Deletion | null;
}

export function allele(doc: VariantDoc, del: Deletion | null): Allele {
  const { start, bp_per_bead: bp, n } = doc.region;
  const keepMask = new Uint8Array(n).fill(1);
  const mid = (k: number) => start + (k + 0.5) * bp;
  if (del) for (let k = 0; k < n; k++) if (mid(k) >= del.start && mid(k) < del.end) keepMask[k] = 0;
  const index = new Int32Array(n).fill(-1);
  const wild: number[] = [];
  for (let k = 0; k < n; k++) if (keepMask[k]) { index[k] = wild.length; wild.push(k); }
  const pos: number[] = [], strand: number[] = [], occ: number[] = [];
  for (const p of doc.ctcf) {
    if (p.motif === null || p.occupancy === null) continue;
    const m = Math.floor((p.start + p.end) / 2);
    const k = Math.floor((m - start) / bp);
    if (k < 0 || k >= n) continue;
    let j = index[k];
    if (j < 0) {
      if (!del?.keep?.some(([a, b]) => a <= m && m < b)) continue;
      j = -1;
      for (let q = k - 1; q >= 0; q--) if (index[q] >= 0) { j = index[q]; break; }
      if (j < 0) j = 0;
    }
    pos.push(j);
    strand.push(p.motif === "+" ? 1 : -1);
    occ.push(p.occupancy);
  }
  const order = pos.map((_, i) => i).sort((a, b) => pos[a] - pos[b] || a - b);
  return {
    n: wild.length,
    mids: Float64Array.from(wild, mid),
    wild: Int32Array.from(wild),
    sites: { pos: Int32Array.from(order, (i) => pos[i]), strand: Int8Array.from(order, (i) => strand[i]), occupancy: Float32Array.from(order, (i) => occ[i]) },
    deletion: del,
  };
}

export const between = (a: Allele, s: number, e: number): number[] => {
  const out: number[] = [];
  for (let i = 0; i < a.n; i++) if (a.mids[i] >= s && a.mids[i] < e) out.push(i);
  return out;
};

export function landmarks(doc: VariantDoc, a: Allele): Record<string, number[]> {
  const out: Record<string, number[]> = { enhancer: between(a, doc.enhancer[0], doc.enhancer[1]) };
  for (const g of doc.genes) out[g.name] = between(a, g.tss - doc.promoter_half, g.tss + doc.promoter_half);
  return out;
}

export interface Prediction {
  allele: Allele;
  lm: Record<string, number[]>;
  /** 4C virtuels (n) : depuis le promoteur de PAX3 et depuis la grappe d'enhancers. */
  from: { PAX3: Float32Array; enhancer: Float32Array };
  /** Contact moyen enhancer → promoteur, par gène. */
  contact: Record<string, number>;
  ms: { extrusion: number; contacts: number };
}

/** L'allèle, `runs` exécutions d'extrusion (graines `seed0 + i`), et ses 4C virtuels. */
export function predict(doc: VariantDoc, del: Deletion | null, runs = 16, seed0 = 1): Prediction {
  const a = allele(doc, del);
  const lm = landmarks(doc, a);
  const t0 = performance.now();
  const sims = Array.from({ length: runs }, (_, i) => simulate(a.n, a.sites, doc.params, seed0 + i));
  const t1 = performance.now();
  const enh = lm.enhancer;
  const pax3 = lm.PAX3;
  const sources = [...enh, ...pax3];
  const R = rows(a.n, sims, sources, doc.law);
  const t2 = performance.now();
  const avg = (ks: number[], off: number) => {
    const out = new Float32Array(a.n);
    for (let q = 0; q < ks.length; q++) for (let j = 0; j < a.n; j++) out[j] += R[(off + q) * a.n + j] / ks.length;
    return out;
  };
  const fromEnh = avg(enh, 0);
  const fromPax3 = avg(pax3, enh.length);
  const contact: Record<string, number> = {};
  for (const g of doc.genes) {
    const ks = lm[g.name];
    contact[g.name] = ks.length ? ks.reduce((s, j) => s + fromEnh[j], 0) / ks.length : NaN;
  }
  return { allele: a, lm, from: { PAX3: fromPax3, enhancer: fromEnh }, contact, ms: { extrusion: t1 - t0, contacts: t2 - t1 } };
}

export function decodeCoords(b64: string): Float32Array {
  const bin = atob(b64);
  const u = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i);
  return new Float32Array(u.buffer);
}
