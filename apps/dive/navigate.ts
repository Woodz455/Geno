/** Semaine 16 : où va la caméra pour une cible de recherche, et ce qu'elle y montre.
 *
 *  Partout, la **bille du noyau** (750 kb) qui porte la cible, sur la copie `a` (l'X actif, `Xa`,
 *  pour le chromosome X du modèle 46,XX) ; l'autre homologue reste surligné ailleurs. Dans la région
 *  fine chargée (TP53 par défaut), plus bas : le monomère de 2 kb, et la paire de bases si elle
 *  est dans la fenêtre où l'hélice se dessine. La profondeur atteinte est rendue, et dite.
 */

import { distanceFor, type Orbit } from "@geno/viewer";
import { type DiveScene, FOV, sideOrbit } from "./dive.js";

export type Depth = "bille" | "fibre" | "base";

export interface Target {
  chrom: string;
  /** 0-based, demi-ouvert : ce qu'on surligne. */
  start: number;
  end: number;
  /** La base qui choisit la bille ou le monomère dont on ouvre la fiche (TSS d'un gène). */
  pos: number;
}

export interface Destination {
  orbit: Orbit;
  depth: Depth;
  /** La fiche à ouvrir : niveau, bille (≥ 1), copie, son intervalle. */
  level: string;
  bead: number;
  copy: string;
  start: number;
  end: number;
  /** Le point visé, µm. */
  centre: [number, number, number];
}

/** Le champ d'une bille : celui de l'arrivée de la semaine 10. */
const BEAD_FIELD = 2.6;
const FIBRE_FIELD = 0.6;
/** Hors de la fenêtre isolée, la fibre s'efface sous 1,4 µm de champ (`SCHEDULE.loopUm`) : on
 *  s'arrête là où elle est encore entière. */
const FIBRE_FIELD_OUTSIDE = 1.5;
const BASE_FIELD = 0.004;

export async function destination(scene: DiveScene, t: Target, from: Orbit): Promise<Destination | null> {
  const R = scene.region;
  const F = scene.fiber;
  if (R && F && t.chrom === R.chrom && t.pos >= R.start && t.pos < R.end) {
    const bp = scene.g3d.level(R.level).bp_per_bead;
    const k = Math.floor((t.pos - R.start) / bp);
    const [b] = (await scene.g3d.locate(R.level, t.chrom, t.pos, t.pos + 1)).find((h) => h.copy === R.copy)?.ids ?? [];
    const bead = b ?? k + 1;
    const info = await scene.g3d.identify(R.level, bead);
    const S = scene.sequence;
    const inWindow = F.windowNear([0, 0, 0]).includes(k);
    if (S && inWindow && t.pos >= S.regionStart) {
      const orbit = sideOrbit(F, t.pos - S.regionStart, BASE_FIELD);
      return { orbit, depth: "base", level: R.level, bead, copy: R.copy, start: info.start, end: info.end,
               centre: orbit.target as [number, number, number] };
    }
    const c: [number, number, number] = [F.points[3 * k], F.points[3 * k + 1], F.points[3 * k + 2]];
    return { orbit: { target: c, yaw: from.yaw, pitch: from.pitch, distance: distanceFor(inWindow ? FIBRE_FIELD : FIBRE_FIELD_OUTSIDE, FOV) },
             depth: "fibre", level: R.level, bead, copy: R.copy, start: info.start, end: info.end, centre: c };
  }
  const lay = scene.fine;
  if (!lay) return null;
  const hits = await scene.g3d.locate(lay.name, t.chrom, t.start, Math.max(t.end, t.start + 1));
  const copy = [`${t.chrom}:a`, `${t.chrom}:Xa`].find((c) => hits.some((h) => h.copy === c)) ?? hits[0]?.copy;
  const span = hits.find((h) => h.copy === copy);
  const at = (await scene.g3d.locate(lay.name, t.chrom, t.pos, t.pos + 1)).find((h) => h.copy === copy);
  if (!span?.ids.length || !at?.ids.length) return null;
  const P = lay.frame.positions;
  const pts = span.ids.map((i) => [P[3 * (i - 1)] / 1000, P[3 * (i - 1) + 1] / 1000, P[3 * (i - 1) + 2] / 1000]);
  const c = [0, 1, 2].map((k) => pts.reduce((s, p) => s + p[k], 0) / pts.length) as [number, number, number];
  const spread = Math.max(...pts.map((p) => Math.hypot(p[0] - c[0], p[1] - c[1], p[2] - c[2])));
  // Un gène de plusieurs billes : le champ les contient toutes.
  const field = Math.min(24, Math.max(BEAD_FIELD, 2.4 * spread + 1));
  const bead = at.ids[0];
  const info = await scene.g3d.identify(lay.name, bead);
  return { orbit: { target: c, yaw: from.yaw, pitch: from.pitch, distance: distanceFor(field, FOV) },
           depth: "bille", level: lay.name, bead, copy, start: info.start, end: info.end, centre: c };
}

/** La profondeur, dite à l'écran. */
export function depthText(d: Destination, scene: DiveScene): string {
  if (d.depth === "base") return "profondeur atteinte : la paire de bases (région fine, séquence hg38)";
  if (d.depth === "fibre") return "profondeur atteinte : le monomère de 2 kb (région fine) — l'hélice ne se dessine que dans la fenêtre de la boucle";
  const fine = scene.region ? ` ; seule la région ${scene.region.chrom}:${(scene.region.start / 1e6).toFixed(2)}–${(scene.region.end / 1e6).toFixed(2)} Mb descend plus bas` : "";
  return `profondeur atteinte : la bille de 750 kb du noyau, copie ${d.copy}${fine}`;
}
