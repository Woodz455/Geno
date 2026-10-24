/** Le chemin de l'ADN le long d'une fibre de monomères : des nucléosomes **qui portent une
 *  longueur exacte** (semaine 12).
 *
 *  Pour qu'une base ait une position, l'ADN doit être une courbe continue dont chaque pas
 *  de base pair est connu. Les nucléosomes de la semaine 11 étaient des cœurs posés au
 *  hasard le long de la fibre ; ici on construit l'ADN d'abord :
 *
 *  - un monomère de `bpPerMonomer` (2 kb) porte `PER_MONOMER` (10) nucléosomes de 147 pb,
 *    reliés par des liens ; il va du point de jonction J(k−1), à mi-chemin du monomère
 *    précédent, au point J(k). Un J est la première base de son monomère ;
 *  - l'ADN enroulé suit une superhélice **gauche** de 1,67 tour, rayon 4,18 nm, pas
 *    2,39 nm (Luger 1997). 146 pas sur cette longueur font 0,30 nm par pas, pas 0,34 :
 *    une hélice B idéale de 147 pb serait plus longue que l'enroulement du cristal. On
 *    garde la géométrie du cœur, et le pas en découle — dit, et vérifié par le test ;
 *  - un lien est une courbe d'Hermite entre la sortie d'un cœur et l'entrée du suivant,
 *    tangente à l'un et à l'autre, dont on règle les tangentes par bissection pour que sa
 *    longueur d'arc soit **exactement** `pas × 0,34 nm` ;
 *  - les longueurs de lien varient (la somme par monomère est fixée) : chaque lien reçoit
 *    au moins ce que sa corde demande, et le reste est partagé. Les vrais liens varient de
 *    20 à 90 pb ; aucune donnée de positionnement n'est lue ici, donc les positions sont
 *    **illustratives** — la géométrie et les longueurs, elles, sont tenues.
 *
 *  Le repère de chaque base : la tangente, et une direction « dyade » qui tourne de
 *  360°/10,5 par pas. Sur l'ADN enroulé, elle part de la normale tournée vers l'octamère ;
 *  sur un lien, d'une normale transportée parallèlement, et l'écart d'angle avec le cœur
 *  suivant est réparti le long du lien — le twist d'un lien s'écarte donc de 34,29° d'au
 *  plus π / pas, pour que la phase se raccorde sans mémoire de toute la région.
 *
 *  **Paresseux, monomère par monomère** : la graine d'un monomère ne dépend que de son rang,
 *  et les cœurs d'un monomère ne dépendent que de ses deux jonctions. On ne calcule que ce
 *  qu'on regarde, et revenir au même endroit montre les mêmes nucléosomes.
 */

export const PER_MONOMER = 10;
export const WRAP_BP = 147;
export const WRAP_TURNS = 1.67;
export const SUPER_RADIUS_NM = 4.18;
export const SUPER_PITCH_NM = 2.39;
export const RISE_NM = 0.34;
export const BP_PER_TURN = 10.5;
/** L'octamère d'histones, sous l'ADN : rayon équatorial et demi-hauteur, nm. */
export const OCTAMER_NM: [number, number] = [3.2, 2.9];
/** Deux cœurs ne s'interpénètrent pas (semaine 11). */
export const MIN_SPACING_NM = 10;

type V = [number, number, number];
const sub = (a: V, b: V): V => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const add = (a: V, b: V): V => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const mul = (a: V, k: number): V => [a[0] * k, a[1] * k, a[2] * k];
const dot = (a: V, b: V) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const len = (a: V) => Math.hypot(a[0], a[1], a[2]);
const norm = (a: V): V => mul(a, 1 / (len(a) || 1));
const cross = (a: V, b: V): V => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];

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

/** Une normale à `t`, définie sans historique : l'axe du monde le moins aligné, projeté. */
function anyNormal(t: V): V {
  const a = Math.abs(t[0]), b = Math.abs(t[1]), c = Math.abs(t[2]);
  const w: V = a <= b && a <= c ? [1, 0, 0] : b <= c ? [0, 1, 0] : [0, 0, 1];
  return norm(cross(t, w));
}

/** Rotation de `v` autour de l'axe unitaire `k` (Rodrigues). */
function rotate(v: V, k: V, angle: number): V {
  const c = Math.cos(angle), s = Math.sin(angle);
  return add(add(mul(v, c), mul(cross(k, v), s)), mul(k, dot(k, v) * (1 - c)));
}

/** Position et tangente d'un lien au paramètre u : Hermite + bosse, dérivée analytique. */
function linkerPoint(L: Linker, u: number): V {
  const u2 = u * u;
  const bump = Math.sin(Math.PI * u) ** 2;
  return add(add(add(mul(L.p0, 2 * u2 * u - 3 * u2 + 1), mul(L.t0, L.lambda * (u2 * u - 2 * u2 + u))),
                 add(mul(L.p1, -2 * u2 * u + 3 * u2), mul(L.t1, L.lambda * (u2 * u - u2)))), mul(L.w, L.beta * bump));
}
function linkerTangent(L: Linker, u: number): V {
  const u2 = u * u;
  const dp = add(add(add(mul(L.p0, 6 * u2 - 6 * u), mul(L.t0, L.lambda * (3 * u2 - 4 * u + 1))),
                     add(mul(L.p1, -6 * u2 + 6 * u), mul(L.t1, L.lambda * (3 * u2 - 2 * u)))),
                 mul(L.w, L.beta * Math.PI * Math.sin(2 * Math.PI * u)));
  return len(dp) < 1e-12 ? norm(sub(L.p1, L.p0)) : norm(dp);
}
/** Longueur d'arc au paramètre u, par la table. */
function arcAt(table: Float64Array, u: number): number {
  const x = u * SAMPLES;
  const q = Math.min(SAMPLES - 1, Math.floor(x));
  return table[q] + (table[q + 1] - table[q]) * (x - q);
}

export interface Core {
  /** Rang global du nucléosome (monomère × PER_MONOMER + rang). */
  index: number;
  monomer: number;
  center: V;
  /** Axe de la superhélice, unitaire. */
  axis: V;
  e1: V;
  e2: V;
  /** Première base enroulée, rang dans la région. */
  firstBp: number;
}

interface Linker {
  firstBp: number;        // la base de départ (sortie d'un cœur, ou jonction), exclue
  steps: number;          // pas jusqu'à la base d'arrivée
  p0: V; t0: V; n0: V;
  p1: V; t1: V; n1: V;
  lambda: number;
  /** Bosse ajoutée, sin²(πu) le long de `w` : allonge sans toucher aux tangentes d'extrémité. */
  beta: number;
  w: V;
  /** Longueurs cumulées aux échantillons de u. */
  table: Float64Array;
  /** Normales transportées aux échantillons, et l'écart final à répartir — calculées à la
   *  première demande d'un repère : le tube de 1 nm n'a besoin que des positions. */
  normals: V[] | null;
  delta: number;
  start: { n: V };
  end: { n: V };
}

interface Monomer {
  cores: Core[];
  /** 11 liens : J(k−1) → cœur 0, …, cœur 9 → J(k). */
  steps: number[];
  linkers: Linker[];
}

const SAMPLES = 128;
/** Coude maximal admis d'une base à la suivante sur un lien, degrés. L'ADN enroulé, lui, tourne
 *  de 4,1° par base autour de l'octamère. */
export const BEND_MAX = 18;

export interface BasePair {
  p: V;
  t: V;
  /** Direction de la dyade, tournée de 360°/10,5 par pas, unitaire, ⟂ t. */
  d: V;
}

export class DnaPath {
  readonly n: number;
  readonly nm: number;            // µm par nm dans le repère des points
  private readonly monomers = new Map<number, Monomer>();
  /** Liens dont la corde dépasse la longueur disponible (attendu : 0). */
  infeasible = 0;
  /** Monomères dont les liens n'ont pas pu recevoir 5 % de mou (attendu : 0). */
  tight = 0;
  /** Monomères dont aucune réalisation ne tient le coude sous BEND_MAX (attendu : 0). */
  bent = 0;

  /** `points` : centres des monomères (xyz entrelacés) ; `unitsPerNm` : unité des points
   *  par nanomètre (0,001 pour des µm). */
  constructor(
    readonly points: Float32Array,
    readonly bpPerMonomer: number,
    readonly monomerRadiusNm: number,
    readonly unitsPerNm = 0.001,
    readonly seed = 11,
  ) {
    this.n = points.length / 3;
    this.nm = unitsPerNm;
    const spare = bpPerMonomer - 1 - PER_MONOMER * WRAP_BP + (PER_MONOMER + 1);
    if (spare < (PER_MONOMER + 1) * 8) throw new Error(`${bpPerMonomer} pb par monomère ne laissent pas de place aux liens`);
  }

  get bpCount(): number {
    return this.n * this.bpPerMonomer;
  }

  private P(k: number): V {
    const i = Math.max(0, Math.min(this.n - 1, k));
    return [this.points[3 * i], this.points[3 * i + 1], this.points[3 * i + 2]];
  }

  /** La spline d'Hermite du tube (tangentes de Catmull-Rom), au paramètre u ∈ [−½, n−½]. */
  spline(u: number): { p: V; t: V } {
    const n = this.n;
    if (n === 1) return { p: this.P(0), t: [1, 0, 0] };
    if (u < 0 || u > n - 1) {
      // Prolongement droit aux deux bouts, dans l'axe du dernier segment.
      const [a, b] = u < 0 ? [0, 1] : [n - 2, n - 1];
      const d = sub(this.P(b), this.P(a));
      const p = add(u < 0 ? this.P(0) : this.P(n - 1), mul(d, u < 0 ? u : u - (n - 1)));
      return { p, t: norm(d) };
    }
    const k = Math.min(n - 2, Math.floor(u));
    const s = u - k;
    const m0 = k === 0 ? sub(this.P(1), this.P(0)) : mul(sub(this.P(k + 1), this.P(k - 1)), 0.5);
    const m1 = k + 2 >= n ? sub(this.P(k + 1), this.P(k)) : mul(sub(this.P(k + 2), this.P(k)), 0.5);
    const s2 = s * s, s3 = s2 * s;
    const p = add(add(mul(this.P(k), 2 * s3 - 3 * s2 + 1), mul(m0, s3 - 2 * s2 + s)),
                  add(mul(this.P(k + 1), -2 * s3 + 3 * s2), mul(m1, s3 - s2)));
    const dp = add(add(mul(this.P(k), 6 * s2 - 6 * s), mul(m0, 3 * s2 - 4 * s + 1)),
                   add(mul(this.P(k + 1), -6 * s2 + 6 * s), mul(m1, 3 * s2 - 2 * s)));
    return { p, t: norm(dp) };
  }

  /** La jonction J(k) : première base du monomère k + 1, à mi-chemin des billes k et k + 1. */
  junction(k: number): { p: V; t: V; n: V } {
    const { p, t } = this.spline(k + 0.5);
    return { p, t, n: anyNormal(t) };
  }

  /** Position, tangente et normale (vers l'octamère) de la base enroulée j ∈ [0, 146]. */
  wrapped(c: Core, j: number): { p: V; t: V; n: V } {
    const turns = WRAP_TURNS * 2 * Math.PI;
    const th = turns * (j / (WRAP_BP - 1) - 0.5);
    const R = SUPER_RADIUS_NM * this.nm;
    const h = (SUPER_PITCH_NM * this.nm) / (2 * Math.PI);
    const radial = add(mul(c.e1, Math.cos(th)), mul(c.e2, Math.sin(th)));
    // Gauche : l'axe recule quand l'angle avance, dans un repère (e1, e2, axe) direct.
    const p = add(add(c.center, mul(radial, R)), mul(c.axis, -h * th));
    const t = norm(add(mul(add(mul(c.e1, -Math.sin(th)), mul(c.e2, Math.cos(th))), R), mul(c.axis, -h)));
    return { p, t, n: mul(radial, -1) };
  }

  /** Les cœurs d'un monomère, posés une fois pour toutes. */
  monomer(k: number): Monomer {
    const hit = this.monomers.get(k);
    if (hit) return hit;
    const nm = this.nm;
    const r = mulberry32(this.seed ^ Math.imul(k + 1, 0x9e3779b1));
    const J0 = this.junction(k - 1);
    const J1 = this.junction(k);
    const jitter = 0.6 * this.monomerRadiusNm * nm;
    const margin = 0.55 * MIN_SPACING_NM * nm;
    const ball = (): V => {
      for (;;) {
        const v: V = [r() * 2 - 1, r() * 2 - 1, r() * 2 - 1];
        if (dot(v, v) <= 1) return v;
      }
    };
    const unit = (): V => {
      for (;;) {
        const v = ball();
        const l = len(v);
        if (l > 1e-3) return mul(v, 1 / l);
      }
    };
    const frame = (axis: V): [V, V] => {
      const e1 = norm(cross(axis, unit()));
      return [e1, cross(axis, e1)];
    };
    const budget = this.bpPerMonomer - 1 - PER_MONOMER * WRAP_BP + (PER_MONOMER + 1);
    const rise = RISE_NM * nm;
    // La corde visée : la longueur moyenne d'un lien, moins 15 % de mou. Viser la corde la
    // plus courte, comme la première version, laissait aux liens deux fois leur corde en
    // longueur : ils se repliaient en épingle, à 115° d'une base à la suivante.
    const SLACK = 1.15;
    const aim0 = (budget * rise) / ((PER_MONOMER + 1) * SLACK);
    // Le coût d'un lien : l'écart de sa corde à la cible, et l'angle de ses deux tangentes
    // d'extrémité avec elle — un lien tangent à sa corde n'a presque pas à se courber.
    const linkCost = (p0: V, t0: V, p1: V, t1: V, aim: number) => {
      const c = sub(p1, p0);
      const l = len(c) || 1e-12;
      const u = mul(c, 1 / l);
      // Et une quasi-interdiction : une tangente d'extrémité qui tourne le dos à la corde fait
      // s'annuler la dérivée de l'Hermite près du bout — un point de rebroussement (97° d'une
      // base à la suivante, vu par le test). Avec λ ≤ corde et deux tangentes dans le sens de
      // la corde, la courbe avance toujours le long de celle-ci et ne rebrousse jamais.
      const back = Math.max(0, 0.3 - dot(t0, u)) + Math.max(0, 0.3 - dot(t1, u));
      return (10 * nm * Math.abs(l - aim)) / aim + 6 * nm * (2 - dot(t0, u) - dot(t1, u)) + 200 * nm * back;
    };

    let best: { cores: Core[]; steps: number[]; linkers: Linker[]; scale: number; bend: number; ok: boolean; impossible: number } | null = null;
    for (let round = 0; round < 8; round++) {
      const aim = aim0 * Math.pow(0.92, round);
      const cores: Core[] = [];
      const chords: number[] = [];
      let prev = J0.p;
      let prevT = J0.t;
      for (let m = 0; m < PER_MONOMER; m++) {
        const T = this.spline(k - 0.5 + (m + 0.5) / PER_MONOMER).p;
        const nextT = m < PER_MONOMER - 1 ? this.spline(k - 0.5 + (m + 1.5) / PER_MONOMER).p : J1.p;
        let pick: { c: Core; chord: number; score: number } | null = null;
        let crowded: { c: Core; chord: number; room: number } | null = null;
        for (let attempt = 0; attempt < 160 && !(pick && attempt >= 16); attempt++) {
          // Si la place manque, on élargit le tirage, jusqu'à 0,75 rayon de monomère : un
          // cœur de 5,5 nm y reste entier.
          const reach = Math.min(0.75 * this.monomerRadiusNm * nm, jitter * (1 + 0.1 * Math.floor(attempt / 20)));
          const center = add(T, mul(ball(), reach));
          // Chaque cœur reste du bon côté des deux plans de jonction, à une demi-distance
          // minimale : deux monomères voisins ne peuvent alors pas se toucher, et chacun se
          // construit sans connaître l'autre — la paresse est gardée.
          if (dot(sub(center, J1.p), J1.t) > -margin || dot(sub(center, J0.p), J0.t) < margin) continue;
          const room = Math.min(Infinity, ...cores.map((o) => len(sub(o.center, center))));
          for (let o = 0; o < 8; o++) {
            const axis = unit();
            const [e1, e2] = frame(axis);
            const c: Core = { index: k * PER_MONOMER + m, monomer: k, center, axis, e1, e2, firstBp: 0 };
            const inn = this.wrapped(c, 0);
            const chord = len(sub(inn.p, prev));
            if (room < MIN_SPACING_NM * nm) {
              if (!crowded || room > crowded.room) crowded = { c, chord, room };
              continue;
            }
            let score = linkCost(prev, prevT, inn.p, inn.t, aim);
            const out = this.wrapped(c, WRAP_BP - 1);
            if (m === PER_MONOMER - 1) {
              score += linkCost(out.p, out.t, J1.p, J1.t, aim);
            } else {
              // Anticipation : la sortie doit regarder vers la cible suivante. Sans elle, le
              // cœur suivant héritait d'une tangente tournée à l'envers (−0,88 sur la corde),
              // qu'aucun de ses tirages ne pouvait plus rattraper.
              const ahead = norm(sub(nextT, out.p));
              score += 6 * nm * (1 - dot(out.t, ahead)) + 200 * nm * Math.max(0, 0.3 - dot(out.t, ahead));
            }
            if (!pick || score < pick.score) pick = { c, chord, score };
          }
        }
        // Aucune place libre : le moins serré des essais, compté.
        const chosen = pick ?? crowded ?? (() => {
          const axis = unit();
          const [e1, e2] = frame(axis);
          const c: Core = { index: k * PER_MONOMER + m, monomer: k, center: T, axis, e1, e2, firstBp: 0 };
          return { c, chord: len(sub(this.wrapped(c, 0).p, prev)) };
        })();
        cores.push(chosen.c);
        chords.push(chosen.chord);
        const out = this.wrapped(chosen.c, WRAP_BP - 1);
        prev = out.p;
        prevT = out.t;
      }
      chords.push(len(sub(J1.p, prev)));
      // Les pas, proportionnels aux cordes, remis à l'échelle du budget du monomère.
      const scale = budget / chords.reduce((s, c) => s + (SLACK * c) / rise, 0);
      const raw = chords.map((c) => Math.max(8, (SLACK * c * scale) / rise));
      const steps = raw.map(Math.floor);
      let extra = budget - steps.reduce((a, b) => a + b, 0);
      // L'arrondi : les pas manquants vont aux plus grands restes, les pas en trop sont repris
      // aux liens qui en ont le plus.
      const order = raw.map((v, i) => [v - Math.floor(v), i]).sort((a, b) => b[0] - a[0]).map((x) => x[1]);
      for (let q = 0; extra > 0; q = (q + 1) % order.length, extra--) steps[order[q]]++;
      while (extra < 0) {
        const i = steps.indexOf(Math.max(...steps));
        steps[i]--;
        extra++;
      }
      let bp = k * this.bpPerMonomer;          // J(k−1)
      cores.forEach((c, m) => {
        bp += steps[m];
        c.firstBp = bp;
        bp += WRAP_BP - 1;
      });
      // Les liens, construits ici et jugés : une réalisation dont un lien se coude de plus de
      // BEND_MAX d'une base à la suivante est retirée, avec une autre graine de tour. Régler
      // les paramètres jusqu'à ce qu'un tirage passe ne garantissait rien — chaque réglage
      // déplaçait le défaut ailleurs.
      const linkers: Linker[] = [];
      let bend = 0;
      let impossible = 0;
      for (let i = 0; i <= PER_MONOMER; i++) {
        const start = i === 0 ? J0 : this.wrapped(cores[i - 1], WRAP_BP - 1);
        const end = i === PER_MONOMER ? J1 : this.wrapped(cores[i], 0);
        const near: V[] = [];
        if (i > 0) near.push(cores[i - 1].center);
        if (i < PER_MONOMER) near.push(cores[i].center);
        const firstBp = i === 0 ? k * this.bpPerMonomer : cores[i - 1].firstBp + WRAP_BP - 1;
        const built = this.buildLinker(start, end, steps[i], near, firstBp);
        linkers.push(built.L);
        bend = Math.max(bend, built.bend);
        if (built.infeasible) impossible++;
      }
      const ok = scale * SLACK >= 1.05 && impossible === 0;
      const better = !best || (ok && !best.ok) || (ok === best.ok && bend < best.bend);
      if (better) best = { cores, steps, linkers, scale, bend, ok, impossible };
      if (ok && bend <= BEND_MAX) break;
    }
    const { cores, steps, linkers, ok, bend, impossible } = best!;
    if (!ok) this.tight++;
    this.infeasible += impossible;
    if (bend > BEND_MAX) this.bent++;
    const mono: Monomer = { cores, steps, linkers };
    this.monomers.set(k, mono);
    return mono;
  }

  /** Les cœurs de tous les monomères de `ks`. */
  cores(ks: Iterable<number>): Core[] {
    const out: Core[] = [];
    for (const k of ks) out.push(...this.monomer(k).cores);
    return out;
  }

  private linker(k: number, i: number): Linker {
    return this.monomer(k).linkers[i]!;
  }

  /** Un lien : la courbe, sa table de longueurs, ses normales, et son plus fort coude (degrés
   *  d'une base à la suivante). */
  private buildLinker(
    start: { p: V; t: V; n: V },
    end: { p: V; t: V; n: V },
    steps: number,
    near: V[],
    firstBp: number,
  ): { L: Linker; bend: number; infeasible: boolean } {
    const nm = this.nm;
    const target = steps * RISE_NM * nm;
    const { p: p0, t: t0 } = start;
    const { p: p1, t: t1 } = end;
    // La courbe : une Hermite de tangentes naturelles (× corde), plus une bosse en sin²(πu)
    // perpendiculaire à la corde, tournée à l'opposé des cœurs voisins. Allonger une Hermite
    // par ses seules tangentes fait des points de rebroussement — la première version y
    // retournait la tangente de 174° d'une base à la suivante ; la bosse, elle, laisse les
    // tangentes d'extrémité intactes et ne rebrousse pas.
    const chord = len(sub(p1, p0));
    const mid = mul(add(p0, p1), 0.5);
    const away = near.length ? sub(mid, mul(near.reduce((a, b) => add(a, b), [0, 0, 0] as V), 1 / near.length)) : anyNormal(sub(p1, p0));
    const cdir = norm(sub(p1, p0));
    // La bosse se fait dans le plan ⟂ corde, dans la direction qui contrarie le moins les deux
    // tangentes d'extrémité (une bosse le long d'une tangente coude le lien à son bout), et du
    // côté opposé aux cœurs.
    const perp = (v: V): V => sub(v, mul(cdir, dot(v, cdir)));
    const a0 = perp(t0), a1 = perp(t1);
    const awayP = perp(away);
    const cands = [cross(cdir, add(a0, a1)), cross(cdir, sub(a0, a1)), awayP, cross(cdir, awayP)]
      .filter((v) => len(v) > 1e-9 * nm)
      .map(norm);
    if (!cands.length) cands.push(anyNormal(cdir));
    let w = cands.reduce((b, v) => (Math.abs(dot(v, t0)) + Math.abs(dot(v, t1)) < Math.abs(dot(b, t0)) + Math.abs(dot(b, t1)) - 1e-9 ? v : b));
    if (dot(w, awayP) < 0) w = mul(w, -1);
    const curve = (lambda: number, beta: number, u: number): V => {
      const u2 = u * u, u3 = u2 * u;
      const bump = Math.sin(Math.PI * u) ** 2;
      return add(add(add(mul(p0, 2 * u3 - 3 * u2 + 1), mul(t0, lambda * (u3 - 2 * u2 + u))),
                     add(mul(p1, -2 * u3 + 3 * u2), mul(t1, lambda * (u3 - u2)))), mul(w, beta * bump));
    };
    // La bissection mesure sur 32 échantillons, la table finale sur 128 : l'écart entre les
    // deux longueurs reste sous 0,1 %, et la génération d'un monomère passe de ~11 à ~2 ms.
    const length = (lambda: number, beta: number, samples = 32) => {
      let s = 0;
      let a = p0;
      for (let q = 1; q <= samples; q++) {
        const b = curve(lambda, beta, q / samples);
        s += len(sub(b, a));
        a = b;
      }
      return s;
    };
    const bisect = (f: (x: number) => number, lo: number, hi: number) => {
      for (let it = 0; it < 20; it++) {
        const m = (lo + hi) / 2;
        if (f(m) < target) lo = m;
        else hi = m;
      }
      return (lo + hi) / 2;
    };
    let lambda = chord;
    let beta = 0;
    const infeasible = chord >= target;
    if (infeasible) {
      lambda = 0;
    } else if (length(chord, 0) >= target) {
      lambda = bisect((x) => length(x, 0), 0, chord);
    } else {
      let hi = target;
      while (length(chord, hi) < target) hi *= 2;
      beta = bisect((x) => length(chord, x), 0, hi);
    }
    const table = new Float64Array(SAMPLES + 1);
    let a = p0;
    for (let q = 1; q <= SAMPLES; q++) {
      const b = curve(lambda, beta, q / SAMPLES);
      table[q] = table[q - 1] + len(sub(b, a));
      a = b;
    }
    // Le coude par base, sur la dérivée analytique : l'angle entre deux tangentes voisines,
    // rapporté à 0,34 nm de chemin.
    const L0: Linker = { firstBp, steps, p0, t0, n0: start.n, p1, t1, n1: end.n, lambda, beta, w, table, normals: null, delta: 0, start, end };
    let bend = 0;
    const COARSE = 64;
    let tPrev = t0;
    for (let q = 1; q <= COARSE; q++) {
      const u = q / COARSE;
      const t = q === COARSE ? t1 : linkerTangent(L0, u);
      const ds = Math.max(1e-30, arcAt(table, u) - arcAt(table, (q - 1) / COARSE));
      bend = Math.max(bend, (Math.acos(Math.max(-1, Math.min(1, dot(tPrev, t)))) * RISE_NM * nm) / ds);
      tPrev = t;
    }
    return { L: L0, bend: (bend * 180) / Math.PI, infeasible };
  }

  /** Les normales transportées d'un lien, à la première demande d'un repère. */
  private frames(L: Linker): void {
    if (L.normals) return;
    const normals: V[] = [L.start.n];
    let n = L.start.n;
    let tPrev = L.t0;
    for (let q = 1; q <= SAMPLES; q++) {
      const t = q === SAMPLES ? L.t1 : linkerTangent(L, q / SAMPLES);
      const axis = cross(tPrev, t);
      const s = len(axis);
      if (s > 1e-12) n = rotate(n, mul(axis, 1 / s), Math.atan2(s, dot(tPrev, t)));
      n = norm(sub(n, mul(t, dot(n, t))));
      normals.push(n);
      tPrev = t;
    }
    // L'écart, autour de t1, entre la normale transportée et celle qu'attend l'arrivée.
    const nT = normals[SAMPLES];
    L.delta = Math.atan2(dot(cross(nT, L.end.n), L.t1), dot(nT, L.end.n));
    L.normals = normals;
  }


  /** La base `b` (rang dans la région) : position, tangente, dyade. */
  basePair(b: number): BasePair {
    const k = Math.max(0, Math.min(this.n - 1, Math.floor(b / this.bpPerMonomer)));
    const mono = this.monomer(k);
    const twist = (2 * Math.PI * b) / BP_PER_TURN;
    const finish = (p: V, t: V, nref: V): BasePair => {
      const d = rotate(nref, t, twist);
      return { p, t, d: norm(sub(d, mul(t, dot(d, t)))) };
    };
    for (let m = 0; m < PER_MONOMER; m++) {
      const c = mono.cores[m];
      if (b >= c.firstBp && b < c.firstBp + WRAP_BP) {
        const w = this.wrapped(c, b - c.firstBp);
        return finish(w.p, w.t, w.n);
      }
    }
    const { L, s, u, q, f } = this.onLinker(k, b);
    const p = linkerPoint(L, u);
    const t = linkerTangent(L, u);
    this.frames(L);
    const n0 = L.normals![q - 1], n1 = L.normals![q];
    let nref = norm(add(mul(n0, 1 - f), mul(n1, f)));
    nref = norm(sub(nref, mul(t, dot(nref, t))));
    nref = rotate(nref, t, L.delta * (s / L.steps));
    return finish(p, t, nref);
  }

  /** Le lien qui porte la base `b`, et son paramètre. */
  private onLinker(k: number, b: number): { L: Linker; s: number; u: number; q: number; f: number } {
    const mono = this.monomer(k);
    let i = 0;
    while (i < PER_MONOMER && b >= mono.cores[i].firstBp) i++;
    const L = this.linker(k, i);
    const s = b - L.firstBp;                       // 0 : la base de départ (jonction)
    const want = (s / L.steps) * L.table[SAMPLES];
    let q = 1;
    while (q < SAMPLES && L.table[q] < want) q++;
    const f = (want - L.table[q - 1]) / Math.max(1e-30, L.table[q] - L.table[q - 1]);
    return { L, s, u: (q - 1 + f) / SAMPLES, q, f };
  }

  /** La position seule de la base `b` : sans repère, donc sans transport de normales. */
  basePoint(b: number): V {
    const k = Math.max(0, Math.min(this.n - 1, Math.floor(b / this.bpPerMonomer)));
    const mono = this.monomer(k);
    for (const c of mono.cores) {
      if (b >= c.firstBp && b < c.firstBp + WRAP_BP) return this.wrapped(c, b - c.firstBp).p;
    }
    const { L, u } = this.onLinker(k, b);
    return linkerPoint(L, u);
  }

  /** Les bases `[b0, b1)` en un tableau : p, t, d (9 flottants par base). */
  basePairs(b0: number, b1: number): Float32Array {
    const out = new Float32Array((b1 - b0) * 9);
    for (let b = b0; b < b1; b++) {
      const { p, t, d } = this.basePair(b);
      out.set(p, 9 * (b - b0));
      out.set(t, 9 * (b - b0) + 3);
      out.set(d, 9 * (b - b0) + 6);
    }
    return out;
  }

  /** L'axe de l'ADN de `[b0, b1)` échantillonné tous les `every` pas, pour le tube de 1 nm. */
  axis(b0: number, b1: number, every = 6): Float32Array {
    const pts: number[] = [];
    for (let b = b0; b < b1; b += every) pts.push(...this.basePoint(b));
    pts.push(...this.basePoint(b1 - 1));
    return Float32Array.from(pts);
  }
}
