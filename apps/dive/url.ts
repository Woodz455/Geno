/** Semaine 16 : la vue dans l'URL, et les signets.
 *
 *  Le fragment porte tout ce qui décide de l'image — l'orbite, la sélection, la surimpression,
 *  la fiche ouverte et son repère — et rien d'autre : `?target=` (la région fine chargée) reste
 *  dans la requête, et la page le relit au chargement comme avant.
 *
 *  Les nombres s'écrivent par `String(x)`, la plus courte écriture qui se relit **au bit près**
 *  en double précision : l'orbite décodée est celle qui a été encodée, pas une approximation.
 *
 *    #v=1&o=x,y,z,lacet,tangage,distance&s=chr17:43044291-43170245&q=BRCA1&l=clinvar
 *        &c=noyau.1234.43044291.-1&w=bille+du+noyau&p=bille+de+750+kb&m=x,y,z
 */

import type { Orbit } from "@geno/viewer";
import type { Region } from "./dive.js";

export const URL_VERSION = 1;

export interface ViewState {
  orbit: Orbit;
  /** Ce qui est surligné (haploïde, 0-based demi-ouvert), ou null. */
  sel: Region[] | null;
  /** La requête qui l'a produit, rendue au champ de recherche. */
  query: string;
  /** La surimpression clinique : "", "clinvar", "gwas", "loeuf". */
  layer: string;
  /** La fiche ouverte : niveau, bille, base (−1 si aucune), et ce que le clic en disait — l'objet
   *  (« base », « nucléosome »…), le brin du barreau, la précision et le détail de la fiche. */
  card: CardState | null;
  /** Le repère du point cliqué, µm. */
  marker: number[] | null;
}

export interface CardState {
  level: string;
  bead: number;
  pos: number;
  what: string;
  strand: 0 | 1 | -1;
  precision: string;
  detail: string;
}

const LAYERS = new Set(["", "clinvar", "gwas", "loeuf"]);
const nums = (s: string, n: number): number[] | null => {
  const v = s.split(",").map(Number);
  return v.length === n && v.every(Number.isFinite) ? v : null;
};

export function encode(s: ViewState): string {
  const o = s.orbit;
  const p = new URLSearchParams();
  p.set("v", String(URL_VERSION));
  p.set("o", [...o.target, o.yaw, o.pitch, o.distance].map(String).join(","));
  if (s.sel?.length) p.set("s", s.sel.map((r) => `${r.chrom}:${r.start}-${r.end}`).join(";"));
  if (s.query) p.set("q", s.query);
  if (s.layer) p.set("l", s.layer);
  if (s.card) {
    p.set("c", `${s.card.level}.${s.card.bead}.${s.card.pos}.${s.card.strand}`);
    p.set("w", s.card.what);
    p.set("p", s.card.precision);
    if (s.card.detail) p.set("d", s.card.detail);
  }
  if (s.marker) p.set("m", s.marker.map(String).join(","));
  return p.toString();
}

/** Le fragment → l'état, ou null s'il ne décrit pas une vue (absent, autre version, illisible). */
export function decode(hash: string): ViewState | null {
  const p = new URLSearchParams(hash.replace(/^#/, ""));
  if (p.get("v") !== String(URL_VERSION)) return null;
  const o = nums(p.get("o") ?? "", 6);
  if (!o) return null;
  const sel: Region[] = [];
  for (const part of (p.get("s") ?? "").split(";").filter(Boolean)) {
    const m = /^(chr[0-9XYM]{1,2}):(\d+)-(\d+)$/.exec(part);
    if (!m) return null;
    sel.push({ chrom: m[1], start: Number(m[2]), end: Number(m[3]) });
  }
  const layer = p.get("l") ?? "";
  if (!LAYERS.has(layer)) return null;
  let card: ViewState["card"] = null;
  const c = p.get("c");
  if (c) {
    const m = /^([\w-]+)\.(\d+)\.(-?\d+)\.(-1|0|1)$/.exec(c);
    if (!m) return null;
    card = { level: m[1], bead: Number(m[2]), pos: Number(m[3]), strand: Number(m[4]) as 0 | 1 | -1,
             what: p.get("w") ?? "", precision: p.get("p") ?? "", detail: p.get("d") ?? "" };
  }
  const marker = p.has("m") ? nums(p.get("m")!, 3) : null;
  if (p.has("m") && !marker) return null;
  return {
    orbit: { target: [o[0], o[1], o[2]], yaw: o[3], pitch: o[4], distance: o[5] },
    sel: sel.length ? sel : null, query: p.get("q") ?? "", layer, card, marker,
  };
}

// — Signets : un nom et un fragment, dans le stockage de ce navigateur seulement —

export interface Bookmark {
  name: string;
  hash: string;
  at: number;
}

const KEY = "geno.signets.v1";

export function loadBookmarks(): Bookmark[] {
  try {
    const v = JSON.parse(localStorage.getItem(KEY) ?? "[]") as Bookmark[];
    return Array.isArray(v) ? v.filter((b) => typeof b?.name === "string" && typeof b?.hash === "string") : [];
  } catch {
    return [];
  }
}

export function saveBookmarks(list: Bookmark[]): boolean {
  try {
    localStorage.setItem(KEY, JSON.stringify(list));
    return true;
  } catch {
    return false;
  }
}
