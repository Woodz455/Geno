/** Nucléosomes **illustratifs** le long de la fibre de 2 kb.
 *
 *  Ce n'est pas un modèle. Aucune donnée de positionnement des nucléosomes n'est lisible
 *  ici (pas de MNase-seq, pas de Micro-C), et le modèle fin s'arrête au monomère de 2 kb.
 *  Ce qui est juste, c'est l'échelle et la forme : environ un nucléosome par 200 pb (ordre
 *  de grandeur de la longueur de répétition), un cœur de 11 nm sur 5,5 nm, un lien d'ADN
 *  entre deux cœurs, et une **chaîne désordonnée** confinée dans le volume du monomère —
 *  pas une fibre régulière de 30 nm, que la ChromEMT n'observe pas *in situ* (Ou 2017).
 *
 *  Déterministe : la graine de chaque nucléosome dépend de son monomère et de son rang,
 *  pas de l'ordre de génération ni de la caméra. Revenir au même endroit montre les mêmes
 *  nucléosomes ; zoomer ne les fait pas scintiller.
 */

export const PER_MONOMER = 10;
export const CORE_RADIUS_UM = 0.0055;       // 11 nm de diamètre
export const CORE_HALF_THICKNESS_UM = 0.00275;
export const MIN_SPACING_UM = 0.010;        // deux cœurs ne s'interpénètrent pas

export interface Nucleosomes {
  centers: Float32Array;
  axes: Float32Array;
  radii: Float32Array;
  /** Monomère porteur de chaque nucléosome. */
  monomer: Uint32Array;
  count: number;
}

function mulberry32(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** Les nucléosomes de toute la polyligne : `PER_MONOMER` par monomère, sur la même spline
 *  d'Hermite que le tube, décalés au hasard dans 55 % du rayon du monomère. */
export function nucleosomes(points: Float32Array, monomerRadiusUm: number, seed = 11): Nucleosomes {
  const n = points.length / 3;
  const total = Math.max(0, n - 1) * PER_MONOMER;
  const centers = new Float32Array(total * 3);
  const axes = new Float32Array(total * 3);
  const radii = new Float32Array(total * 2);
  const monomer = new Uint32Array(total);
  const P = (k: number, c: number) => points[3 * Math.max(0, Math.min(n - 1, k)) + c];
  let prev: number[] | null = null;
  let q = 0;
  for (let k = 0; k < n - 1; k++) {
    const rnd = mulberry32(seed ^ Math.imul(k + 1, 0x9e3779b1));
    for (let m = 0; m < PER_MONOMER; m++) {
      const u = (m + 0.5) / PER_MONOMER;
      const u2 = u * u, u3 = u2 * u;
      const h00 = 2 * u3 - 3 * u2 + 1, h10 = u3 - 2 * u2 + u, h01 = -2 * u3 + 3 * u2, h11 = u3 - u2;
      const base = [0, 1, 2].map((c) => {
        const m0 = k === 0 ? P(1, c) - P(0, c) : 0.5 * (P(k + 1, c) - P(k - 1, c));
        const m1 = k + 2 >= n ? P(k + 1, c) - P(k, c) : 0.5 * (P(k + 2, c) - P(k, c));
        return h00 * P(k, c) + h10 * m0 + h01 * P(k + 1, c) + h11 * m1;
      });
      // Un décalage dans la boule de 55 % du rayon, et une distance minimale au précédent :
      // une chaîne de cœurs qui se touchent ne se lit plus comme une chaîne.
      let c = base;
      for (let attempt = 0; attempt < 12; attempt++) {
        let x, y, z, l;
        do {
          x = rnd() * 2 - 1;
          y = rnd() * 2 - 1;
          z = rnd() * 2 - 1;
          l = x * x + y * y + z * z;
        } while (l > 1);
        const s = 0.55 * monomerRadiusUm;
        const cand = [base[0] + x * s, base[1] + y * s, base[2] + z * s];
        c = cand;
        if (!prev || Math.hypot(cand[0] - prev[0], cand[1] - prev[1], cand[2] - prev[2]) >= MIN_SPACING_UM) break;
      }
      let ax, ay, az, al;
      do {
        ax = rnd() * 2 - 1;
        ay = rnd() * 2 - 1;
        az = rnd() * 2 - 1;
        al = Math.hypot(ax, ay, az);
      } while (al > 1 || al < 1e-3);
      centers.set(c, 3 * q);
      axes.set([ax / al, ay / al, az / al], 3 * q);
      radii.set([CORE_RADIUS_UM, CORE_HALF_THICKNESS_UM], 2 * q);
      monomer[q] = k;
      prev = c;
      q++;
    }
  }
  return { centers, axes, radii, monomer, count: q };
}
