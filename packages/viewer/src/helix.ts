/** La double hélice B, procédurale, autour d'un chemin d'ADN (semaine 12).
 *
 *  Une base pair est donnée par son axe : position `p`, tangente `t`, et la direction de sa
 *  dyade `d` (⟂ t), qui tourne de 360°/10,5 d'un pas au suivant. Les deux phosphates de la
 *  paire sont posés par les coordonnées cylindriques du **modèle de fibre B d'Arnott** :
 *  rayon 8,91 Å, angle ±94,9° autour de la dyade, et ±2,85 Å le long de l'axe. C'est ce
 *  décalage, pas un réglage, qui ouvre un grand et un petit sillon : sur l'hélice idéale à
 *  10,5 pb/tour, le test mesure 5,4 Å et 11,8 Å (P–P le plus court à travers le sillon,
 *  moins 5,8 Å), pour 5,7 et 11,7 Å publiés.
 *
 *  Dessin, en tubes (`TubeRenderer`) : deux squelettes qui passent par les phosphates, et
 *  un barreau par paire, du phosphate d'un brin à celui de l'autre, coupé en son milieu —
 *  chaque moitié a la couleur de sa base. L'identifiant d'une moitié (et du squelette au
 *  même pas) est `2 × rang + brin` : le picking nomme une base, pas une paire.
 *
 *  La géométrie est générée **pour la seule fenêtre demandée** : quelques milliers de
 *  paires, jamais la région entière.
 */

import { TubeRenderer } from "./tube.js";

export const B_DNA = {
  bpPerTurn: 10.5,
  riseNm: 0.34,
  phosphateRadiusNm: 0.891,
  phosphateAngleDeg: 94.9,
  phosphateRiseNm: 0.285,
  /** Diamètre de van der Waals d'un phosphate, retranché pour une largeur de sillon. */
  phosphateDiameterNm: 0.58,
};

const BACKBONE_NM = 0.21;
const RUNG_NM = 0.16;

/** Palette des bases : Okabe–Ito, lisible par les trois daltonismes courants. */
export const BASE_COLORS: Record<string, [number, number, number]> = {
  A: [0.0, 0.62, 0.45],      // vert bleuté
  C: [0.34, 0.71, 0.91],     // bleu ciel
  G: [0.94, 0.89, 0.26],     // jaune
  T: [0.84, 0.37, 0.0],      // vermillon
  N: [0.5, 0.5, 0.5],
};
const BACKBONE = [[0.78, 0.8, 0.84], [0.62, 0.65, 0.72]] as const;
const COMPLEMENT: Record<string, string> = { A: "T", C: "G", G: "C", T: "A", N: "N" };
export const complement = (b: string) => COMPLEMENT[b] ?? "N";

/** Les deux phosphates de chaque paire. `frames` : p, t, d (9 flottants par paire). */
export function phosphates(frames: Float32Array, unitsPerNm: number): { s1: Float32Array; s2: Float32Array } {
  const n = frames.length / 9;
  const s1 = new Float32Array(n * 3);
  const s2 = new Float32Array(n * 3);
  const a = (B_DNA.phosphateAngleDeg * Math.PI) / 180;
  const r = B_DNA.phosphateRadiusNm * unitsPerNm;
  const z = B_DNA.phosphateRiseNm * unitsPerNm;
  for (let i = 0; i < n; i++) {
    const o = 9 * i;
    const p = [frames[o], frames[o + 1], frames[o + 2]];
    const t = [frames[o + 3], frames[o + 4], frames[o + 5]];
    const d = [frames[o + 6], frames[o + 7], frames[o + 8]];
    const e = [t[1] * d[2] - t[2] * d[1], t[2] * d[0] - t[0] * d[2], t[0] * d[1] - t[1] * d[0]];
    for (let c = 0; c < 3; c++) {
      s1[3 * i + c] = p[c] + r * (Math.cos(a) * d[c] + Math.sin(a) * e[c]) + z * t[c];
      s2[3 * i + c] = p[c] + r * (Math.cos(a) * d[c] - Math.sin(a) * e[c]) - z * t[c];
    }
  }
  return { s1, s2 };
}

export interface HelixWindow {
  /** Rang dans la région de la première paire de chaque plage, et le nombre de paires. */
  ranges: [number, number][];
  /** p, t, d par paire, plages mises bout à bout. */
  frames: Float32Array;
  /** Bases du brin +, une lettre par paire, plages mises bout à bout. */
  bases: string;
  /** Rang de la base à mettre en évidence, ou −1. */
  highlight?: number;
}

export interface HelixDrawOptions {
  view: Float32Array;
  height: number;
  fade?: number;
  invert?: boolean;
  idBase?: number;
}

export class HelixRenderer {
  private readonly backbone: TubeRenderer;
  private readonly rungs: TubeRenderer;
  private empty = true;
  /** Phosphates de la dernière fenêtre, pour l'étiquetage et la mesure. */
  last: { s1: Float32Array; s2: Float32Array; window: HelixWindow } | null = null;

  constructor(gl: WebGL2RenderingContext, private readonly unitsPerNm = 0.001) {
    this.backbone = new TubeRenderer(gl);
    this.rungs = new TubeRenderer(gl);
  }

  upload(w: HelixWindow): void {
    const { s1, s2 } = phosphates(w.frames, this.unitsPerNm);
    this.last = { s1, s2, window: w };
    const n = w.bases.length;
    this.empty = n === 0;
    if (this.empty) return;
    const idOf = (q: number, strand: 0 | 1) => {
      let rank = 0;
      let acc = 0;
      for (const [b0, len] of w.ranges) {
        if (q < acc + len) {
          rank = b0 + (q - acc);
          break;
        }
        acc += len;
      }
      return 2 * rank + strand + 1;
    };
    const ids = new Uint32Array(n * 2);
    for (let q = 0; q < n; q++) {
      ids[2 * q] = idOf(q, 0);
      ids[2 * q + 1] = idOf(q, 1);
    }

    // Squelettes : une polyligne par brin et par plage.
    {
      const pts = new Float32Array(n * 6);
      const radii = new Float32Array(n * 2).fill(BACKBONE_NM * this.unitsPerNm);
      const colors = new Float32Array(n * 6);
      const bid = new Uint32Array(n * 2);
      const breaks: number[] = [];
      let o = 0;
      let acc = 0;
      for (const [, len] of w.ranges) {
        for (const strand of [0, 1] as const) {
          if (o) breaks.push(o);
          const s = strand ? s2 : s1;
          for (let q = acc; q < acc + len; q++, o++) {
            pts.set(s.subarray(3 * q, 3 * q + 3), 3 * o);
            colors.set(BACKBONE[strand], 3 * o);
            bid[o] = ids[2 * q + strand];
          }
        }
        acc += len;
      }
      this.backbone.upload({ points: pts, radii, colors, ids: bid, breaks });
    }

    // Barreaux : deux demi-barreaux par paire, du phosphate au milieu.
    {
      const pts = new Float32Array(n * 12);
      const radii = new Float32Array(n * 4);
      const colors = new Float32Array(n * 12);
      const rid = new Uint32Array(n * 4);
      const breaks: number[] = [];
      for (let q = 0; q < n; q++) {
        const a = s1.subarray(3 * q, 3 * q + 3);
        const b = s2.subarray(3 * q, 3 * q + 3);
        const m = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2];
        const base = w.bases[q];
        const rank = (ids[2 * q] - 1) >> 1;
        const hot = w.highlight === rank;
        const r = (hot ? 1.35 : 1) * RUNG_NM * this.unitsPerNm;
        for (const [h, from, letter] of [[0, a, base], [1, b, complement(base)]] as const) {
          const o = 4 * q + 2 * h;
          if (o) breaks.push(o);
          pts.set(from, 3 * o);
          pts.set(m, 3 * (o + 1));
          radii[o] = radii[o + 1] = r;
          let c = BASE_COLORS[letter] ?? BASE_COLORS.N;
          if (hot) c = [Math.min(1, c[0] * 1.25 + 0.15), Math.min(1, c[1] * 1.25 + 0.15), Math.min(1, c[2] * 1.25 + 0.15)];
          colors.set(c, 3 * o);
          colors.set(c, 3 * (o + 1));
          rid[o] = rid[o + 1] = ids[2 * q + h];
        }
      }
      this.rungs.upload({ points: pts, radii, colors, ids: rid, breaks });
    }
  }

  draw(proj: Float32Array, o: HelixDrawOptions): void {
    if (this.empty) return;
    const opts = { view: o.view, height: o.height, fade: o.fade ?? 1, invert: o.invert ?? false, idBase: o.idBase ?? 0 };
    this.backbone.draw(proj, opts);
    this.rungs.draw(proj, opts);
  }

  triangles(): number {
    return this.empty ? 0 : this.backbone.triangles() + this.rungs.triangles();
  }
}

/** Décode un identifiant de base : rang dans la région et brin (0 : +, 1 : −). */
export function baseOf(id: number): { rank: number; strand: 0 | 1 } {
  const v = id - 1;
  return { rank: v >> 1, strand: (v & 1) as 0 | 1 };
}
