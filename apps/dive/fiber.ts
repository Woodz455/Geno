/** Le niveau fin dessiné : la fibre de 2 kb, ses sites CTCF, ses cohésines, ses boucles.
 *
 *  Ce qu'on veut lire à l'échelle d'une boucle CTCF :
 *  - **la fibre**, un tube continu le long des 2 000 monomères ;
 *  - **les sites CTCF**, en flèches qui enveloppent la fibre et pointent dans le sens du
 *    motif : `+` vers les coordonnées croissantes, `−` vers les décroissantes. Deux flèches
 *    qui se font face, c'est une paire convergente — la règle de la semaine 8 se lit ;
 *  - **les cohésines tenues dans cette conformation**, en anneaux autour de leurs deux pieds ;
 *  - **les paires d'ancres de l'ensemble**, en arcs bombés hors de la chromatine, d'autant
 *    plus épais qu'une cohésine les tient souvent sur les 1 200 conformations. Un arc n'est
 *    pas une liaison de la structure affichée : c'est une fréquence de population, dessinée
 *    là où sont ses deux ancres dans cette structure-ci ;
 *  - **les nucléosomes**, illustratifs, quand on descend sous 250 nm de champ.
 *
 *  Tout est dans le repère du noyau (la pose de la semaine 9), en µm.
 */

import { type G3D, toNucleus } from "@geno/g3d";
import { EllipsoidRenderer, TubeRenderer } from "@geno/viewer";
import { nucleosomes } from "./nucleosomes.js";

export const IDS = {
  fiber: 2_000_000,
  glyph: 3_000_000,
  ring: 3_100_000,
  arc: 3_200_000,
  nucleosome: 4_000_000,
};

/** Seuil des arcs, fixé avant toute mesure : une paire tenue dans au moins 5 % des
 *  conformations est une boucle de la population. */
export const ARC_THRESHOLD = 0.05;

const COLORS = {
  fiber: [0.62, 0.66, 0.72],
  loop: [0.95, 0.62, 0.25],
  ctcf: [0.96, 0.82, 0.2],
  cohesin: [0.3, 0.82, 0.78],
  arc: [0.72, 0.5, 0.95],
  core: [0.86, 0.78, 0.68],
  dna: [0.55, 0.62, 0.7],
};

export interface Site {
  monomer: number;
  strand: 1 | -1;
  occupancy: number;
  /** Sens de la flèche, monde, unitaire. */
  dir: [number, number, number];
}

export interface Loop {
  i: number;
  j: number;
  kind: string;
  /** Centre des monomères de la boucle, µm, repère du noyau. */
  centroid: [number, number, number];
  /** Plus grande distance d'un monomère de la boucle à son centre, µm. */
  extent: number;
}

export interface FiberState {
  /** 0 → 1 : la fibre remplace les billes de 750 kb. */
  fiber: number;
  /** 0 → 1 : la boucle visée s'isole et ses annotations apparaissent. */
  loop: number;
  /** 0 → 1 : le tube cède aux nucléosomes. */
  nuc: number;
}

export interface Fiber {
  n: number;
  /** Positions des monomères, µm, repère du noyau. */
  points: Float32Array;
  radiusUm: number;
  sites: Site[];
  /** Boucles tenues par une cohésine dans la conformation affichée. */
  held: Loop[];
  /** Paires d'ancres au-dessus du seuil, avec leur fréquence. */
  arcs: { i: number; j: number; weight: number }[];
  /** La boucle visée : CTCF des deux côtés, la plus proche d'ACTB. */
  focus: Loop | null;
  draw(proj: Float32Array, view: Float32Array, height: number, s: FiberState): void;
  describe(id: number): Promise<string | null>;
  /** Nucléosomes générés (après la première descente sous 250 nm). */
  nucleosomeCount(): number;
}

const sub = (a: number[], b: number[]) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const add = (a: number[], b: number[]) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const mul = (a: number[], k: number) => [a[0] * k, a[1] * k, a[2] * k];
const norm = (a: number[]) => {
  const l = Math.hypot(a[0], a[1], a[2]) || 1;
  return [a[0] / l, a[1] / l, a[2] / l];
};
const cross = (a: number[], b: number[]) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];

export async function openFiber(
  gl: WebGL2RenderingContext,
  g3d: G3D,
  name = "fin",
  nearMonomer = 763,
  isolate = true,
): Promise<Fiber | null> {
  if (!g3d.header.levels.some((l) => l.name === name)) return null;
  const frame = await g3d.firstFrame(name);
  const lv = frame.level;
  const n = frame.count;
  const points = toNucleus(lv, frame.positions).map((v) => v / 1000);
  const radiusUm = (lv.radius_nm[0] / 1000) * 0.75;
  const P = (k: number) => [points[3 * k], points[3 * k + 1], points[3 * k + 2]];
  const tangent = (k: number) => norm(sub(P(Math.min(n - 1, k + 1)), P(Math.max(0, k - 1))));

  // — Sites CTCF —
  const fwd = await g3d.attribute(name, "ctcf_forward");
  const rev = await g3d.attribute(name, "ctcf_reverse");
  const sites: Site[] = [];
  for (let k = 0; k < n; k++) {
    if (fwd && Number.isFinite(fwd[k])) sites.push({ monomer: k, strand: 1, occupancy: fwd[k], dir: tangent(k) as Site["dir"] });
    if (rev && Number.isFinite(rev[k])) sites.push({ monomer: k, strand: -1, occupancy: rev[k], dir: mul(tangent(k), -1) as Site["dir"] });
  }

  // — Boucles et paires —
  const L = await g3d.links(name);
  const loopOf = (i: number, j: number, kind: string): Loop => {
    let c = [0, 0, 0];
    for (let k = i; k <= j; k++) c = add(c, P(k));
    c = mul(c, 1 / (j - i + 1));
    let extent = 0;
    for (let k = i; k <= j; k++) extent = Math.max(extent, Math.hypot(...sub(P(k), c)));
    return { i, j, kind, centroid: c as Loop["centroid"], extent };
  };
  const held: Loop[] = [];
  const arcs: Fiber["arcs"] = [];
  if (L) {
    for (let k = 0; k < L.i.length; k++) {
      const kind = L.kinds[L.kind[k]];
      const i = Math.min(L.i[k], L.j[k]);
      const j = Math.max(L.i[k], L.j[k]);
      if (kind === "ctcf_pair") {
        if (L.weight[k] >= ARC_THRESHOLD) arcs.push({ i, j, weight: L.weight[k] });
      } else held.push(loopOf(i, j, kind));
    }
  }
  const dist = (l: Loop) => (nearMonomer < l.i ? l.i - nearMonomer : nearMonomer > l.j ? nearMonomer - l.j : 0);
  // La boucle à montrer est une boucle CTCF **des deux façons** : tenue par une cohésine dont
  // les deux pieds sont sur des CTCF dans la conformation affichée, *et* fréquente dans la
  // population (son arc passe le seuil). La première version ne demandait que la première
  // condition et montrait 1341–1392 — un domaine convergent planté, tenu ici, mais tenu dans
  // 2,5 % des conformations seulement : une boucle de l'instant, pas de la population.
  const pairWeight = (l: Loop) =>
    Math.max(0, ...arcs.filter((a) => Math.abs(a.i - l.i) <= 2 && Math.abs(a.j - l.j) <= 2).map((a) => a.weight));
  const ctcfHeld = held.filter((l) => l.kind === "cohesin_ctcf").sort((a, b) => dist(a) - dist(b));
  const focus = ctcfHeld.find((l) => pairWeight(l) >= ARC_THRESHOLD) ?? ctcfHeld[0] ?? null;
  // La fenêtre de la boucle visée, avec 6 monomères (12 kb) de contexte de part et d'autre
  // pour que ses ancres ne flottent pas dans le vide — 15 masquaient la boucle de face. Fibre **et annotations** s'isolent ensemble :
  // la première version n'effaçait que la fibre, et les 40 flèches, 20 anneaux et 13 arcs de
  // toute la région flottaient autour de la boucle, sans chromatine sous eux.
  const MARGIN = 6;
  // `isolate = false` est le **témoin** de la mesure de lisibilité : mêmes annotations, rien
  // d'effacé autour.
  const inWindow = (k: number) => !isolate || !focus || (k >= focus.i - MARGIN && k <= focus.j + MARGIN);

  // — La fibre : deux palettes, territoire et « boucle visée » —
  const fiber = new TubeRenderer(gl);
  const c1 = new Float32Array(n * 3);
  const c2 = new Float32Array(n * 3);
  for (let k = 0; k < n; k++) {
    c1.set(COLORS.fiber, 3 * k);
    c2.set(focus && k >= focus.i && k <= focus.j ? COLORS.loop : COLORS.fiber, 3 * k);
  }
  fiber.upload({ points, radii: new Float32Array(n).fill(radiusUm), colors: c1, colors2: c2 });

  // — Flèches CTCF : un manchon qui enveloppe la fibre, puis une pointe —
  const glyphs = new TubeRenderer(gl);
  const glyphOwner: number[] = [];      // 1 : dans la fenêtre de la boucle visée
  {
    const pts: number[] = [];
    const radii: number[] = [];
    const ids: number[] = [];
    const breaks: number[] = [];
    const R = radiusUm * 1.45;
    sites.forEach((s, q) => {
      if (q) breaks.push(pts.length / 3);
      glyphOwner.push(...[0, 0, 0, 0].map(() => (inWindow(s.monomer) ? 1 : 0)));
      const p = P(s.monomer);
      const d = s.dir;
      // Avancée d'1,2 rayon dans le sens du motif : aux deux pieds d'une boucle tenue, les
      // sites se touchent, et deux flèches centrées sur eux se fondaient en une tache. Ainsi
      // décalées, chacune se pose sur le brin qu'elle désigne.
      for (const [off, r] of [[-0.4, R], [1.4, R], [1.45, R * 1.35], [3.1, 0]] as [number, number][]) {
        pts.push(...add(p, mul(d, off * radiusUm)));
        radii.push(Math.max(r, 1e-6));
        ids.push(q + 1);
      }
    });
    if (sites.length) {
      const m = radii.length;
      glyphs.upload({
        points: Float32Array.from(pts),
        radii: Float32Array.from(radii),
        colors: Float32Array.from(Array.from({ length: m }, () => COLORS.ctcf).flat()),
        ids: Uint32Array.from(ids),
        breaks,
      });
    }
  }

  // — Anneaux de cohésine autour des deux pieds des boucles tenues —
  const rings = new TubeRenderer(gl);
  const ringOwner: number[] = [];
  {
    const pts: number[] = [];
    const ids: number[] = [];
    const breaks: number[] = [];
    const SEG = 32;
    held.forEach((l, q) => {
      if (q) breaks.push(pts.length / 3);
      // Toutes les cohésines de la fenêtre : un même domaine peut en tenir deux à la fois.
      const keep = inWindow(l.i) && inWindow(l.j) ? 1 : 0;
      const a = P(l.i);
      const b = P(l.j);
      const c = mul(add(a, b), 0.5);
      // Les deux brins traversent l'anneau en sens contraires : son axe suit leur différence.
      let ax = norm(sub(tangent(l.i), tangent(l.j)));
      if (!Number.isFinite(ax[0])) ax = [0, 0, 1];
      const u0 = norm(Math.abs(ax[0]) < 0.9 ? cross(ax, [1, 0, 0]) : cross(ax, [0, 1, 0]));
      const v0 = cross(ax, u0);
      const r = Math.hypot(...sub(a, b)) / 2 + radiusUm * 1.6;
      for (let s = 0; s <= SEG; s++) {
        const th = (2 * Math.PI * s) / SEG;
        pts.push(...add(c, add(mul(u0, r * Math.cos(th)), mul(v0, r * Math.sin(th)))));
        ids.push(q + 1);
        ringOwner.push(keep);
      }
    });
    if (held.length) {
      const m = ids.length;
      rings.upload({
        points: Float32Array.from(pts),
        radii: new Float32Array(m).fill(radiusUm * 0.25),
        colors: Float32Array.from(Array.from({ length: m }, () => COLORS.cohesin).flat()),
        ids: Uint32Array.from(ids),
        breaks,
      });
    }
  }

  // — Arcs : Bézier quadratique bombée hors de la boucle, épaisseur ∝ fréquence —
  const arcTube = new TubeRenderer(gl);
  const arcOwner: number[] = [];
  {
    const pts: number[] = [];
    const radii: number[] = [];
    const ids: number[] = [];
    const breaks: number[] = [];
    const SEG = 28;
    const wmax = Math.max(...arcs.map((a) => a.weight), 1e-6);
    arcs.forEach((a, q) => {
      if (q) breaks.push(pts.length / 3);
      const p0 = P(a.i);
      const p1 = P(a.j);
      const mid = mul(add(p0, p1), 0.5);
      const g = loopOf(a.i, a.j, "").centroid;
      let out = sub(mid, g);
      if (Math.hypot(...out) < 1e-6) out = cross(sub(p1, p0), [0, 0, 1]);
      out = norm(out);
      // Bombé d'au moins la moitié de l'étendue de la boucle : les deux ancres d'une boucle
      // tenue se touchent, et un arc proportionnel à leur seul écart disparaissait entre elles.
      const lp = loopOf(a.i, a.j, "");
      const h = Math.max(0.6 * Math.hypot(...sub(p1, p0)), 0.5 * lp.extent) + 4 * radiusUm;
      const keep = inWindow(a.i) && inWindow(a.j) ? 1 : 0;
      const ctrl = add(mid, mul(out, h));
      const r = radiusUm * (0.25 + 0.75 * (a.weight / wmax));
      for (let s = 0; s <= SEG; s++) {
        const t = s / SEG;
        const b = add(add(mul(p0, (1 - t) * (1 - t)), mul(ctrl, 2 * t * (1 - t))), mul(p1, t * t));
        pts.push(...b);
        radii.push(r);
        ids.push(q + 1);
        arcOwner.push(keep);
      }
    });
    if (arcs.length) {
      const m = ids.length;
      arcTube.upload({
        points: Float32Array.from(pts),
        radii: Float32Array.from(radii),
        colors: Float32Array.from(Array.from({ length: m }, () => COLORS.arc).flat()),
        ids: Uint32Array.from(ids),
        breaks,
      });
    }
  }

  // — Nucléosomes et ADN de liaison, générés à la première descente —
  let nucs: EllipsoidRenderer | null = null;
  let linker: TubeRenderer | null = null;
  let nucMonomer: Uint32Array | null = null;
  let nucOwner: number[] = [];
  let nucLoop = -1;
  let nucCount = 0;
  const makeNucleosomes = () => {
    const nu = nucleosomes(points, lv.radius_nm[0] / 1000);
    nucCount = nu.count;
    nucMonomer = nu.monomer;
    const colors = new Float32Array(nu.count * 3);
    for (let q = 0; q < nu.count; q++) {
      const k = nu.monomer[q];
      colors.set(focus && k >= focus.i && k <= focus.j ? COLORS.loop : COLORS.core, 3 * q);
    }
    nucOwner = Array.from(nu.monomer, (k) => (inWindow(k) ? 1 : 0));
    nucs = new EllipsoidRenderer(gl);
    nucs.upload({
      centers: nu.centers,
      axes: nu.axes,
      radii: nu.radii,
      ids: Uint32Array.from({ length: nu.count }, (_, q) => IDS.nucleosome + q + 1),
      colors,
      count: nu.count,
    });
    linker = new TubeRenderer(gl);
    linker.upload({
      points: nu.centers,
      radii: new Float32Array(nu.count).fill(0.001),
      colors: Float32Array.from(Array.from({ length: nu.count }, () => COLORS.dna).flat()),
    });
  };

  const fade = new Float32Array(n);
  let lastLoop = -1;
  const ownerFade = (owner: number[], loop: number) => owner.map((o) => (o ? 1 : 1 - loop));

  return {
    n,
    points,
    radiusUm,
    sites,
    held,
    arcs,
    focus,
    nucleosomeCount: () => nucCount,

    draw(proj, view, height, s) {
      if (s.fiber <= 0) return;
      if (s.loop !== lastLoop) {
        for (let k = 0; k < n; k++) fade[k] = inWindow(k) ? 1 : 1 - s.loop;
        fiber.setFade(fade);
        if (glyphOwner.length) glyphs.setFade(ownerFade(glyphOwner, s.loop));
        if (ringOwner.length) rings.setFade(ownerFade(ringOwner, s.loop));
        if (arcOwner.length) arcTube.setFade(ownerFade(arcOwner, s.loop));
        lastLoop = s.loop;
      }
      // La fibre complète le noyau (qui garde h ≥ fiber) puis cède aux nucléosomes.
      const tubeFade = s.fiber * (1 - s.nuc);
      if (tubeFade > 0) fiber.draw(proj, { view, height, fade: tubeFade, mix: s.loop, idBase: IDS.fiber });
      // Les annotations apparaissent avec la boucle. Flèches et arcs sont des symboles à
      // l'échelle de la fibre : ils s'effacent quand on descend aux nucléosomes, où une flèche
      // de 29 nm écraserait des cœurs de 11 nm. L'anneau de cohésine, lui, est un objet réel
      // de cette taille-là : il reste.
      if (s.loop > 0) {
        const symbols = s.loop * (1 - s.nuc);
        if (symbols > 0) glyphs.draw(proj, { view, height, fade: symbols, idBase: IDS.glyph });
        rings.draw(proj, { view, height, fade: s.loop, idBase: IDS.ring });
        if (symbols > 0) arcTube.draw(proj, { view, height, fade: symbols, idBase: IDS.arc });
      }
      if (s.nuc > 0) {
        if (!nucs) makeNucleosomes();
        if (s.loop !== nucLoop) {
          const a = Float32Array.from(ownerFade(nucOwner, s.loop));
          nucs!.setAlpha(a);
          linker!.setFade(a);
          nucLoop = s.loop;
        }
        nucs!.draw(proj, { view, fade: s.nuc, invert: true });
        linker!.draw(proj, { view, height, fade: s.nuc, invert: true, idBase: IDS.nucleosome });
      }
    },

    async describe(id) {
      if (id >= IDS.nucleosome) {
        const q = id - IDS.nucleosome - 1;
        const k = nucMonomer ? nucMonomer[Math.max(0, q)] : 0;
        return `nucléosome illustratif sur le monomère ${k} — position tirée, pas mesurée`;
      }
      if (id >= IDS.arc) {
        const a = arcs[id - IDS.arc - 1];
        return a ? `paire d'ancres CTCF ${a.i}–${a.j} (${((a.j - a.i) * lv.bp_per_bead / 1000).toFixed(0)} kb) : tenue par une cohésine dans ${(a.weight * 100).toFixed(1)} % des conformations` : null;
      }
      if (id >= IDS.ring) {
        const l = held[id - IDS.ring - 1];
        const what = { cohesin_ctcf: "pieds arrêtés sur deux CTCF convergents", cohesin_one_side: "un pied arrêté", cohesin_free: "en cours d'extrusion" }[l?.kind ?? ""] ?? "";
        return l ? `cohésine tenant ${l.i}–${l.j} dans cette conformation (${what})` : null;
      }
      if (id >= IDS.glyph) {
        const st = sites[id - IDS.glyph - 1];
        return st ? `site CTCF ${st.strand > 0 ? "+" : "−"} sur le monomère ${st.monomer}, occupation ${st.occupancy.toFixed(2)} (planté)` : null;
      }
      if (id >= IDS.fiber) {
        const k = id - IDS.fiber - 1;
        const b = await g3d.identify(name, k + 1);
        return `${b.copy}:${b.start.toLocaleString("fr")}–${b.end.toLocaleString("fr")} (fibre, 2 kb) · dispersion ${b.variabilityNm.toFixed(0)} nm`;
      }
      return null;
    },
  };
}
