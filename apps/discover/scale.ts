/** Les échelles comparées du mode Découverte (semaine 19).
 *
 *  Chaque analogie est **calculée** depuis des nombres sourcés, jamais écrite à la main : le texte
 *  porte des marques `{clé}` que `fill()` remplace, et `tools/test-discover.mjs` recalcule chaque
 *  valeur par un autre chemin. Changer un nombre ici change le texte partout.
 */

/** Les nombres de départ, avec leur source. Longueurs en mètres. */
export const FACTS = {
  /** Génome diploïde de GM12878 (46,XX) : 2 × ~3,1 Gb (hg38 sans chrY). */
  diploidBp: 2 * 3.1e9,
  /** Montée de la double hélice B, par paire de bases (Watson et Crick ; Dickerson 1982). */
  risePerBp: 0.34e-9,
  /** Diamètre du noyau du modèle (semaine 6) et d'un noyau de lymphocyte, ordre de grandeur. */
  nucleus: 10e-6,
  /** Diamètre de la double hélice B. */
  helix: 2e-9,
  /** Une bille du noyau (semaine 6). */
  beadBp: 750_000,
  /** Un monomère de la fibre (semaine 8). */
  monomerBp: 2_000,
  /** ADN enroulé sur un nucléosome, et ses tours de superhélice gauche (Luger et al. 1997). */
  nucleosomeBp: 147,
  nucleosomeTurns: 1.65,
  /** Taille médiane d'un domaine de contact dans GM12878 (Rao et al. 2014). */
  tadBp: 185_000,
  /** Une balle de tennis (règlement ITF : 6,54–6,86 cm), une page de roman (~1 800 signes),
   *  un roman (~300 pages). */
  tennisBall: 0.067,
  charsPerPage: 1_800,
  pagesPerNovel: 300,
} as const;

const fmt = (x: number, digits = 0): string =>
  x.toLocaleString("fr-FR", { maximumFractionDigits: digits, minimumFractionDigits: digits });

/** La longueur de l'ADN diploïde déplié, en mètres. */
export const dnaLength = (): number => FACTS.diploidBp * FACTS.risePerBp;

/** Toutes les valeurs que les textes peuvent citer, déjà mises en forme. */
export function values(): Record<string, string> {
  const L = dnaLength();
  const ball = FACTS.tennisBall / FACTS.nucleus;           // grossissement noyau → balle
  const novels = (bp: number) => bp / (FACTS.charsPerPage * FACTS.pagesPerNovel);
  return {
    dna_m: fmt(L, 1),
    nucleus_um: fmt(FACTS.nucleus * 1e6),
    fold: fmt(L / FACTS.nucleus),                          // fois plus long que le noyau est large
    ball_dna_km: fmt((L * ball) / 1000, 0),                // l'ADN, si le noyau était une balle de tennis
    ball_helix_um: fmt(FACTS.helix * ball * 1e6, 1),       // l'épaisseur du fil, à cette échelle
    bead_kb: fmt(FACTS.beadBp / 1000),
    bead_novels: fmt(novels(FACTS.beadBp), 1),
    genome_novels: fmt(novels(FACTS.diploidBp / 2)),
    tad_kb: fmt(FACTS.tadBp / 1000),
    tad_um: fmt(FACTS.tadBp * FACTS.risePerBp * 1e6, 0),   // un TAD, déplié
    monomer_kb: fmt(FACTS.monomerBp / 1000),
    nuc_bp: fmt(FACTS.nucleosomeBp),
    nuc_turns: fmt(FACTS.nucleosomeTurns, 2),
    rise_nm: fmt(FACTS.risePerBp * 1e9, 2),
    helix_nm: fmt(FACTS.helix * 1e9),
  };
}

/** Remplace chaque `{clé}` par sa valeur ; une clé inconnue est une erreur, pas un trou. */
export function fill(text: string, v: Record<string, string> = values()): string {
  return text.replace(/\{([a-z_]+)\}/g, (_, k: string) => {
    if (!(k in v)) throw new Error(`valeur d'échelle inconnue : {${k}}`);
    return v[k];
  });
}

/** Le champ de la caméra, lisible : « 2,6 µm », « 120 nm », « 4,0 nm ». */
export function fieldText(widthUm: number): string {
  if (widthUm >= 1) return `${fmt(widthUm, widthUm < 10 ? 1 : 0)} µm`;
  const nm = widthUm * 1000;
  return `${fmt(nm, nm < 10 ? 1 : 0)} nm`;
}

/** Ce que mesure le champ, en objets familiers : combien de largeurs de noyau, de bille, d'hélice. */
export function fieldCompare(widthUm: number): string {
  const m = widthUm * 1e-6;
  if (m >= FACTS.nucleus / 2) return `${fmt(m / FACTS.nucleus, 1)} fois la largeur du noyau`;
  if (m >= 50e-9) return `${fmt(m / FACTS.nucleus * 100, 1)} % de la largeur du noyau`;
  return `${fmt(m / FACTS.helix, 0)} fois l'épaisseur de la double hélice`;
}
