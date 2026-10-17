/** « Une boucle CTCF se lit comme telle, ancres comprises, à l'échelle où elle existe. »
 *
 *  Traduit en propriétés mesurables sur l'image de la vue de boucle, par le tampon
 *  d'identifiants du picking — ce qui est mesuré est ce qui est **visible**, pas ce qui est
 *  dessiné, puisque le tampon ne garde que le fragment le plus proche :
 *
 *  1. les deux **ancres** — la flèche CTCF `+` au pied gauche, la `−` au pied droit — ont
 *     chacune au moins `MIN_PX` pixels non occultés ;
 *  2. l'**anneau** de cohésine à la base, et l'**arc** de la paire, aussi ;
 *  3. au moins 80 % des monomères du **corps** de la boucle ont un pixel visible ;
 *  4. les deux flèches suivent leur brin **vers l'intérieur** de la boucle, en 3D et à l'écran
 *     — et celles de tous les arcs dessinés aussi.
 *
 *  Les seuils sont fixés ici, avant la mesure. Le témoin est la même vue sans isoler la
 *  boucle du reste de la région (`noIsolate`) : il dit ce que l'isolement achète.
 */

import type { DiveScene } from "./dive.js";
import { IDS } from "./fiber.js";

export const MIN_PX = 50;
export const MIN_BODY = 0.8;

export interface Legibility {
  ok: boolean;
  loop: { i: number; j: number; kb: number } | null;
  anchorPx: [number, number];
  ringPx: number;
  arcPx: number;
  arcWeight: number | null;
  body: number;
  inward: [boolean, boolean];
  inwardDot: [number, number];
  /** Le même produit scalaire en 3D, avant projection : sépare une erreur de flèche d'une
   *  erreur de mesure. */
  inward3d: [number, number];
  /** Ancres de tous les arcs dessinés, et celles dont la flèche suit son brin vers l'intérieur. */
  arcAnchors: number;
  arcAnchorsInward: number;
}

export function measureLoop(scene: DiveScene): Legibility {
  const F = scene.fiber;
  const focus = F?.focus;
  const none: Legibility = { ok: false, loop: null, anchorPx: [0, 0], ringPx: 0, arcPx: 0, arcWeight: null, body: 0, inward: [false, false], inwardDot: [0, 0], inward3d: [0, 0], arcAnchors: 0, arcAnchorsInward: 0 };
  if (!F || !focus || !scene.loopView) return none;
  scene.arrival = 1;
  scene.render(scene.loopView);
  const ids = scene.ids();

  const count = new Map<number, number>();
  for (const v of ids) count.set(v, (count.get(v) ?? 0) + 1);
  const px = (id: number) => count.get(id) ?? 0;

  const siteAt = (k: number, strand: 1 | -1) =>
    F.sites.findIndex((s) => s.strand === strand && Math.abs(s.monomer - k) <= 1);
  const left = siteAt(focus.i, 1);
  const right = siteAt(focus.j, -1);
  const anchorPx: [number, number] = [left >= 0 ? px(IDS.glyph + left + 1) : 0, right >= 0 ? px(IDS.glyph + right + 1) : 0];

  let ringPx = 0;
  F.held.forEach((l, q) => {
    if (Math.abs(l.i - focus.i) <= 2 && Math.abs(l.j - focus.j) <= 2) ringPx += px(IDS.ring + q + 1);
  });
  const arcIndex = F.arcs.findIndex((a) => Math.abs(a.i - focus.i) <= 2 && Math.abs(a.j - focus.j) <= 2);
  const arcPx = arcIndex >= 0 ? px(IDS.arc + arcIndex + 1) : 0;

  let seen = 0;
  for (let k = focus.i; k <= focus.j; k++) if (px(IDS.fiber + k + 1) > 0) seen++;
  const body = seen / (focus.j - focus.i + 1);

  // Orientation : chaque flèche doit suivre **son brin vers l'intérieur de la boucle** —
  // un site `+` au pied gauche vers les coordonnées croissantes, un `−` au pied droit vers
  // les décroissantes ; deux flèches qui se font face le long de la chaîne, c'est la règle
  // de convergence. On compare la flèche à la corde d'un monomère dans la boucle, en 3D et à
  // l'écran, pour les deux ancres de **chaque** arc dessiné. Une erreur de signe sur le brin
  // — la faute que la semaine 8 gardait en vue dans le 1D — mettrait tout en négatif.
  //
  // Une première version comparait la flèche à la corde vers le monomère situé trois plus
  // loin : sur la boucle visée, la chaîne fait demi-tour en 6 kb à sa base, et la mesure
  // donnait −0,35 et −0,65 en 3D. Elle mesurait la pelote, pas la flèche.
  const P = (k: number) => [F.points[3 * k], F.points[3 * k + 1], F.points[3 * k + 2]];
  const along = (q: number, inward: 1 | -1) => {
    const st = F.sites[q];
    const p = P(st.monomer);
    const t = P(Math.max(0, Math.min(F.n - 1, st.monomer + inward)));
    const d = [t[0] - p[0], t[1] - p[1], t[2] - p[2]];
    const l = Math.hypot(d[0], d[1], d[2]) || 1;
    const d3 = (st.dir[0] * d[0] + st.dir[1] * d[1] + st.dir[2] * d[2]) / l;
    const a = scene.project(p);
    const b = scene.project([p[0] + st.dir[0] * F.radiusUm, p[1] + st.dir[1] * F.radiusUm, p[2] + st.dir[2] * F.radiusUm]);
    const c = scene.project(t);
    const u = [b[0] - a[0], b[1] - a[1]];
    const v = [c[0] - a[0], c[1] - a[1]];
    const d2 = (u[0] * v[0] + u[1] * v[1]) / ((Math.hypot(u[0], u[1]) * Math.hypot(v[0], v[1])) || 1);
    return { d3, d2 };
  };
  const inwardDot = [left >= 0 ? along(left, 1).d2 : -1, right >= 0 ? along(right, -1).d2 : -1] as [number, number];
  const inward3d = [left >= 0 ? along(left, 1).d3 : -1, right >= 0 ? along(right, -1).d3 : -1] as [number, number];
  const inward: [boolean, boolean] = [inward3d[0] > 0 && inwardDot[0] > 0, inward3d[1] > 0 && inwardDot[1] > 0];
  let arcAnchors = 0;
  let arcAnchorsInward = 0;
  for (const a of F.arcs) {
    const l = siteAt(a.i, 1);
    const r = siteAt(a.j, -1);
    for (const [q, dir] of [[l, 1], [r, -1]] as [number, 1 | -1][]) {
      if (q < 0) continue;
      arcAnchors++;
      if (along(q, dir).d3 > 0) arcAnchorsInward++;
    }
  }

  return {
    ok: anchorPx[0] >= MIN_PX && anchorPx[1] >= MIN_PX && ringPx >= MIN_PX && arcPx >= MIN_PX && body >= MIN_BODY && inward[0] && inward[1],
    loop: { i: focus.i, j: focus.j, kb: (focus.j - focus.i) * 2 },
    anchorPx,
    ringPx,
    arcPx,
    arcWeight: arcIndex >= 0 ? F.arcs[arcIndex].weight : null,
    body,
    inward,
    inwardDot,
    inward3d,
    arcAnchors,
    arcAnchorsInward,
  };
}
