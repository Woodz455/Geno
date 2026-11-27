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

/** Les langues du parcours (semaine 20). Les nombres sont calculés une fois, mis en forme par langue. */
export type Locale = "fr-FR" | "en-GB";

const fmt = (x: number, digits = 0, locale: Locale = "fr-FR"): string =>
  x.toLocaleString(locale, { maximumFractionDigits: digits, minimumFractionDigits: digits });

/** La longueur de l'ADN diploïde déplié, en mètres. */
export const dnaLength = (): number => FACTS.diploidBp * FACTS.risePerBp;

/** Toutes les valeurs que les textes peuvent citer, déjà mises en forme. */
export function values(locale: Locale = "fr-FR"): Record<string, string> {
  const f = (x: number, d = 0) => fmt(x, d, locale);
  const L = dnaLength();
  const ball = FACTS.tennisBall / FACTS.nucleus;           // grossissement noyau → balle
  const novels = (bp: number) => bp / (FACTS.charsPerPage * FACTS.pagesPerNovel);
  return {
    dna_m: f(L, 1),
    nucleus_um: f(FACTS.nucleus * 1e6),
    fold: f(L / FACTS.nucleus),                          // fois plus long que le noyau est large
    ball_dna_km: f((L * ball) / 1000, 0),                // l'ADN, si le noyau était une balle de tennis
    ball_helix_um: f(FACTS.helix * ball * 1e6, 1),       // l'épaisseur du fil, à cette échelle
    bead_kb: f(FACTS.beadBp / 1000),
    bead_novels: f(novels(FACTS.beadBp), 1),
    genome_novels: f(novels(FACTS.diploidBp / 2)),
    tad_kb: f(FACTS.tadBp / 1000),
    tad_um: f(FACTS.tadBp * FACTS.risePerBp * 1e6, 0),   // un TAD, déplié
    monomer_kb: f(FACTS.monomerBp / 1000),
    nuc_bp: f(FACTS.nucleosomeBp),
    nuc_turns: f(FACTS.nucleosomeTurns, 2),
    rise_nm: f(FACTS.risePerBp * 1e9, 2),
    helix_nm: f(FACTS.helix * 1e9),
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
export function fieldText(widthUm: number, locale: Locale = "fr-FR"): string {
  if (widthUm >= 1) return `${fmt(widthUm, widthUm < 10 ? 1 : 0, locale)} µm`;
  const nm = widthUm * 1000;
  return `${fmt(nm, nm < 10 ? 1 : 0, locale)} nm`;
}

const COMPARE: Record<Locale, [string, string, string]> = {
  "fr-FR": ["{x} fois la largeur du noyau", "{x} % de la largeur du noyau", "{x} fois l'épaisseur de la double hélice"],
  "en-GB": ["{x} times the width of the nucleus", "{x} % of the width of the nucleus", "{x} times the thickness of the double helix"],
};

/** Ce que mesure le champ, en objets familiers : combien de largeurs de noyau, de bille, d'hélice. */
export function fieldCompare(widthUm: number, locale: Locale = "fr-FR"): string {
  const m = widthUm * 1e-6;
  const [nuc, pct, hel] = COMPARE[locale];
  if (m >= FACTS.nucleus / 2) return nuc.replace("{x}", fmt(m / FACTS.nucleus, 1, locale));
  if (m >= 50e-9) return pct.replace("{x}", fmt(m / FACTS.nucleus * 100, 1, locale));
  return hel.replace("{x}", fmt(m / FACTS.helix, 0, locale));
}
