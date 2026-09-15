"""Ce qu'on a le droit de mesurer sur un ensemble, et ce qui n'a pas de sens.

Le piège de la semaine est géométrique. Deux noyaux recuits indépendamment n'ont
**aucun repère commun** : ni orientation — rien ne distingue un axe dans une
sphère — ni même attribution des territoires, puisque chr7:a atterrit ailleurs à
chaque tirage. « L'écart-type de la position de la bille i » en x, y, z ne veut
donc rien dire, et un alignement de Procruste ne le sauve pas : il n'y a rien à
aligner. C'est mesuré ici plutôt qu'affirmé (`frame`).

Ce qui reste, et qui est invariant par rotation :

- la **position radiale** de chaque bille ;
- les **distances par paires**, donc la carte de contacts et P(s) ;
- les grandeurs par copie — rayon de giration, position radiale du centroïde.

Tout le reste du module ne manipule que ça.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree


# --------------------------------------------------------------------------
# Le repère n'existe pas
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Frame:
    """Ce que l'alignement optimal de deux structures rattrape. Réponse : presque rien.

    Trois écarts, tous après rotation optimale :

    - `aligned` — la vraie correspondance, bille i contre bille i ;
    - `radial_shuffled` — les billes permutées **à l'intérieur de leur couche
      radiale**, ce qui détruit leur identité mais garde exactement leur profondeur ;
    - `shuffled` — les billes permutées n'importe comment, la référence sans information.

    `kept` est la part de l'écart de référence qui **survit** au meilleur
    alignement possible. Proche de 1, il n'y a pas de repère commun.
    """

    aligned: float            # RMSD après alignement optimal, µm
    radial_shuffled: float    # … en permutant à l'intérieur des couches radiales
    shuffled: float           # … en permutant librement
    nuclear_radius: float
    pairs: int

    @property
    def kept(self) -> float:
        return self.aligned / self.shuffled

    @property
    def radial_explains(self) -> float:
        """Part du gain de l'alignement qu'on retrouve en ne gardant que les rayons."""
        gain = self.shuffled - self.aligned
        return (self.shuffled - self.radial_shuffled) / gain if gain > 1e-12 else float("nan")

    def __str__(self) -> str:
        return (
            f"RMSD aligné {self.aligned:.2f} µm dans un noyau de rayon "
            f"{self.nuclear_radius:.1f} µm · référence sans information {self.shuffled:.2f} µm "
            f"→ {self.kept:.0%} de l'écart survit au meilleur alignement · en ne gardant "
            f"que les rayons on en retrouve {self.radial_explains:.0%}"
        )

def _kabsch(a: np.ndarray, b: np.ndarray) -> float:
    """RMSD de `b` sur `a` après la rotation optimale. La réflexion est permise —
    une matrice de distances ne détermine la structure qu'à une isométrie près,
    et la semaine 5 avait déjà tranché ce point."""
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    u, _, vt = np.linalg.svd(b.T @ a)
    return float(np.sqrt((((b @ (u @ vt)) - a) ** 2).sum(axis=1).mean()))


def _shuffle_within_shells(
    x: np.ndarray, rng: np.random.Generator, shells: int
) -> np.ndarray:
    """Permute les billes sans changer l'ensemble des profondeurs occupées."""
    order = np.argsort(np.linalg.norm(x, axis=1))
    out = np.arange(len(x))
    for block in np.array_split(order, shells):
        out[block] = block[rng.permutation(len(block))]
    return x[out]


def frame(
    coords: np.ndarray,
    nuclear_radius: float,
    pairs: int = 20,
    seed: int = 0,
    shells: int = 20,
) -> Frame:
    """Mesure, plutôt que d'affirmer, qu'il n'y a pas de repère commun.

    Le témoin « couches radiales » a été ajouté sur une hypothèse — que le peu
    que l'alignement rattrape serait la stratification radiale partagée — et
    **cette hypothèse est fausse** : mêler les billes à profondeur constante
    ramène pratiquement à la référence sans information. Ce que la correspondance
    porte est donc ailleurs, dans la compacité des territoires : deux structures
    ont chacune 46 blobs, et une rotation bien choisie en superpose quelques-uns.
    Le témoin est resté parce qu'il dit ça.
    """
    rng = np.random.default_rng(seed)
    n, k = coords.shape[0], coords.shape[1]
    aligned, radial, free = [], [], []
    for _ in range(pairs):
        s, t = rng.choice(n, 2, replace=False)
        aligned.append(_kabsch(coords[s], coords[t]))
        radial.append(_kabsch(coords[s], _shuffle_within_shells(coords[t], rng, shells)))
        free.append(_kabsch(coords[s], coords[t][rng.permutation(k)]))
    return Frame(
        aligned=float(np.mean(aligned)),
        radial_shuffled=float(np.mean(radial)),
        shuffled=float(np.mean(free)),
        nuclear_radius=nuclear_radius,
        pairs=pairs,
    )


# --------------------------------------------------------------------------
# Reproductibilité de la position radiale
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Reproducibility:
    """Quelle part de la position radiale d'une bille tient à la bille, et non au tirage."""

    icc: float              # part de variance attribuable à la bille
    between_nm: float       # écart-type des moyennes par bille
    within_nm: float        # écart-type moyen d'une bille sur l'ensemble
    pairwise_r: float       # corrélation moyenne entre profils de deux structures
    n_structures: int
    n_beads: int

    def __str__(self) -> str:
        return (
            f"ICC {self.icc:.3f} · dispersion entre billes {self.between_nm:.0f} nm "
            f"contre {self.within_nm:.0f} nm au sein d'une bille · corrélation entre "
            f"deux structures r = {self.pairwise_r:+.3f}"
        )


def reproducibility(
    radial: np.ndarray, nuclear_radius: float, *, pairs: int = 200, seed: int = 0
) -> Reproducibility:
    """Décomposition de variance sur la position radiale, façon ANOVA à un facteur.

    La bille est le sujet, la structure le juge. L'ICC répond à la question qui
    fait le principe n° 1 du projet : **« la bille i est à telle profondeur » est-il
    un énoncé sur la bille, ou sur le tirage ?**

    On donne aussi la corrélation entre les profils radiaux de deux structures
    prises au hasard, qui dit la même chose plus directement et se lit sans
    connaître l'ICC.
    """
    m, n = radial.shape                       # m structures, n billes
    per_bead = radial.mean(axis=0)
    grand = float(radial.mean())

    ms_between = m * ((per_bead - grand) ** 2).sum() / (n - 1)
    ms_within = ((radial - per_bead) ** 2).sum() / (n * (m - 1))
    icc = (ms_between - ms_within) / (ms_between + (m - 1) * ms_within)

    rng = np.random.default_rng(seed)
    rs = []
    for _ in range(pairs):
        s, t = rng.choice(m, 2, replace=False)
        rs.append(np.corrcoef(radial[s], radial[t])[0, 1])

    scale = nuclear_radius * 1000.0
    return Reproducibility(
        icc=float(icc),
        between_nm=float(per_bead.std()) * scale,
        within_nm=float(radial.std(axis=0).mean()) * scale,
        pairwise_r=float(np.mean(rs)),
        n_structures=m,
        n_beads=n,
    )


def medoid(radial: np.ndarray) -> tuple[int, np.ndarray]:
    """Structure la plus centrale **pour la distance entre profils radiaux**.

    La précision compte : il n'existe pas de médoïde « de l'ensemble » dans
    l'absolu, seulement un médoïde pour une distance donnée, et celle-ci ne voit
    que la profondeur des billes. Deux structures aux territoires disposés tout
    autrement peuvent avoir le même profil radial et ce médoïde ne les
    distinguera pas. C'est assumé : la profondeur est ce que l'ensemble mesure.
    """
    centred = radial - radial.mean(axis=0)
    gram = centred @ centred.T
    sq = np.diag(gram)
    d2 = np.maximum(sq[:, None] + sq[None, :] - 2.0 * gram, 0.0)
    total = np.sqrt(d2).sum(axis=1)
    return int(total.argmin()), total


# --------------------------------------------------------------------------
# Contacts
# --------------------------------------------------------------------------


def contact_pairs(x: np.ndarray, radius: np.ndarray, cutoff: float) -> np.ndarray:
    """Paires de billes dont les centres sont à moins de `cutoff × (r_i + r_j)`.

    `cutoff = 1` veut dire « en contact » au sens du volume exclu. Au-delà, on
    élargit le voisinage : une ligature Hi-C n'exige pas que deux nucléosomes se
    touchent, seulement qu'ils soient assez proches pour se retrouver dans le même
    complexe. Le seuil est donc un paramètre, pas une constante, et le rapport le
    balaie.
    """
    pairs = cKDTree(x).query_pairs(
        cutoff * 2.0 * float(radius.max()), output_type="ndarray"
    )
    if len(pairs) == 0:
        return pairs
    i, j = pairs[:, 0], pairs[:, 1]
    d = np.linalg.norm(x[i] - x[j], axis=1)
    return pairs[d < cutoff * (radius[i] + radius[j])]


def slope(separation: np.ndarray, probability: np.ndarray, lo: float, hi: float) -> float:
    """Pente de log P contre log s sur la fenêtre `[lo, hi]` en paires de bases."""
    window = (separation >= lo) & (separation <= hi) & (probability > 0)
    if window.sum() < 3:
        return float("nan")
    return float(np.polyfit(np.log10(separation[window]), np.log10(probability[window]), 1)[0])


# P(s) du modèle n'a pas une pente, elle en a trois. Les fenêtres portent des noms
# parce que chacune dit autre chose (§ VALIDATION S7).
REGIMES: tuple[tuple[str, float, float], ...] = (
    ("chaîne", 1.4e6, 3.8e6),
    ("polymère", 2.2e6, 9.0e6),
    ("territoire", 1.5e7, 1.5e8),
)


@dataclass(frozen=True)
class Contacts:
    """Ce qu'un ensemble prédit d'une expérience Hi-C."""

    cis: int                 # contacts intra-copie
    trans: int               # entre copies différentes
    homolog: int             # entre les deux copies d'un même chromosome
    per_structure: float
    cutoff: float
    separation: np.ndarray   # (k,) séparation génomique médiane du bin, pb
    probability: np.ndarray  # (k,) contacts observés / paires possibles / structures
    exponent: float          # pente sur le régime « polymère », celle qui se compare au Hi-C
    fit_range: tuple[float, float]
    regimes: tuple[tuple[str, float], ...]
    plateau: float           # P(s) moyenne sur le régime « territoire »

    @property
    def trans_fraction(self) -> float:
        return self.trans / max(self.cis + self.trans, 1)

    def __str__(self) -> str:
        return (
            f"{self.per_structure:,.0f} contacts par structure au seuil {self.cutoff:g} · "
            f"trans {self.trans_fraction:.1%} · homologues {self.homolog / max(self.trans, 1):.1%} "
            f"des trans · P(s) ∝ s^{self.exponent:+.2f} entre "
            f"{self.fit_range[0] / 1e6:.1f} et {self.fit_range[1] / 1e6:.0f} Mb, "
            f"puis plateau à {self.plateau:.3f}"
        )


def homolog_map(labels: tuple[str, ...]) -> np.ndarray:
    """Pour chaque copie, l'index de son homologue — déduit des étiquettes.

    Les copies sortent du caryotype par paires consécutives, si bien qu'un simple
    `index ^ 1` donnerait la bonne réponse aujourd'hui. Il la donnerait encore
    sur un caryotype aneuploïde, où elle serait fausse : une trisomie 21 casse
    l'appariement sans casser l'ordre. On lit donc le nom du chromosome.
    """
    chroms = [label.split(":", 1)[0] for label in labels]
    partner = np.arange(len(labels))
    for chrom in set(chroms):
        same = [k for k, c in enumerate(chroms) if c == chrom]
        if len(same) == 2:
            partner[same[0]], partner[same[1]] = same[1], same[0]
    return partner


def contacts(
    coords: np.ndarray,
    radius: np.ndarray,
    copy_id: np.ndarray,
    labels: tuple[str, ...],
    bp_per_bead: float,
    *,
    cutoff: float = 1.5,
    fit_from: float = 2.2e6,
    fit_to: float = 9.0e6,
) -> Contacts:
    """Carte de contacts agrégée sur l'ensemble, et la courbe P(s) qui en sort.

    P(s) n'a été demandée nulle part au modèle : elle tombe de la géométrie —
    chaîne de sphères tangentes, volume exclu, confinement, territoires. C'est
    donc une **prédiction**, et la seule de cet ensemble qu'une expérience Hi-C
    puisse contredire.

    La normalisation compte : à la séparation `s` en billes, le nombre de paires
    *possibles* vaut `Σ_c (k_c − s)` sur les copies. Diviser par le nombre de
    contacts observés sans ça donnerait une décroissance en partie factice, due à
    la seule raréfaction des paires disponibles aux grandes séparations.
    """
    n_struct, n_beads = coords.shape[:2]
    counts = np.bincount(copy_id)
    max_sep = int(counts.max())

    observed = np.zeros(max_sep + 1, dtype=np.int64)
    cis = trans = homolog = 0
    partner = homolog_map(labels)

    for s in range(n_struct):
        pairs = contact_pairs(coords[s], radius, cutoff)
        if len(pairs) == 0:
            continue
        ci, cj = copy_id[pairs[:, 0]], copy_id[pairs[:, 1]]
        same = ci == cj
        cis += int(same.sum())
        trans += int((~same).sum())
        # `& ~same` n'est pas une précaution : une copie sans homologue est
        # son propre partenaire dans la table, si bien que sans ce masque tous
        # ses contacts *cis* seraient comptés comme des contacts d'homologues.
        # Un contact d'homologues est trans par définition.
        homolog += int(((partner[ci] == cj) & ~same).sum())
        sep = pairs[same, 1] - pairs[same, 0]
        observed += np.bincount(sep, minlength=max_sep + 1)

    possible = np.array(
        [np.maximum(counts - s, 0).sum() for s in range(max_sep + 1)], dtype=np.float64
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        prob = observed / (possible * n_struct)

    sep_bp = np.arange(max_sep + 1) * bp_per_bead
    ok = (sep_bp > 0) & np.isfinite(prob) & (prob > 0)
    sep_ok, prob_ok = sep_bp[ok], prob[ok]

    # La séparation 1 vaut exactement P = 1 par construction — deux billes liées
    # sont toujours à portée dès que le seuil dépasse l'allongement maximal de la
    # liaison. L'inclure dans un ajustement mesurerait le modèle de liaison, pas
    # le repliement ; les fenêtres commencent donc au-delà.
    plateau_window = (sep_ok >= REGIMES[2][1]) & (sep_ok <= REGIMES[2][2])

    return Contacts(
        cis=cis,
        trans=trans,
        homolog=homolog,
        per_structure=(cis + trans) / n_struct,
        cutoff=cutoff,
        separation=sep_ok,
        probability=prob_ok,
        exponent=slope(sep_ok, prob_ok, fit_from, fit_to),
        fit_range=(fit_from, fit_to),
        regimes=tuple(
            (name, slope(sep_ok, prob_ok, lo, hi)) for name, lo, hi in REGIMES
        ),
        plateau=float(prob_ok[plateau_window].mean()) if plateau_window.any() else float("nan"),
    )


# --------------------------------------------------------------------------
# Profils radiaux, et la corrélation qu'on ne peut pas encore faire
# --------------------------------------------------------------------------


def radial_by_quartile(
    radial: np.ndarray, lad: np.ndarray, groups: int = 4
) -> list[tuple[float, float, float, float]]:
    """Position radiale par quartile de contenu LAD, découpé **par rang**.

    Découper sur les valeurs ne marche pas : une bille de 750 kb sans le moindre
    LAD est fréquente, si bien que le premier et le deuxième quartile de la
    distribution valent tous deux zéro et que le groupe Q1 ressort vide. Mesuré :
    le rapport affichait `nan`. Le rang, lui, est toujours défini.

    Renvoie par groupe `(rayon moyen, dispersion entre billes, LAD min, LAD max)` —
    les bornes LAD sont là pour que les ex æquo se voient.
    """
    order = np.argsort(lad, kind="stable")
    out = []
    for block in np.array_split(order, groups):
        values = radial[:, block]
        out.append(
            (
                float(values.mean()),
                float(values.mean(axis=0).std()),
                float(lad[block].min()),
                float(lad[block].max()),
            )
        )
    return out


def self_consistency(radial: np.ndarray, lad: np.ndarray) -> float:
    """Corrélation entre la position radiale moyenne et la piste LAD **d'entrée**.

    Ce n'est **pas** une validation : la piste LAD est ce qui a fixé les rayons
    visés du modèle, donc cette corrélation ne teste qu'une chose, que le solveur
    a fait ce qu'on lui a demandé. Elle est ici pour ça et pour rien d'autre.
    La validation, elle, corrèle contre un DamID *mesuré* — voir `damid`.
    """
    return float(np.corrcoef(radial.mean(axis=0), lad)[0, 1])


def read_bedgraph(path) -> tuple[dict[str, list[tuple[int, int, float]]], int]:
    """Lit un bedGraph `chrom start end valeur`. Renvoie la piste et le nombre de lignes.

    Les lignes `track`, `browser` et les commentaires sont sautés ; une valeur non
    numérique fait échouer la lecture plutôt que d'être silencieusement ignorée —
    un DamID à demi lu produirait une corrélation à demi fausse, ce qui est pire
    qu'une erreur.
    """
    track: dict[str, list[tuple[int, int, float]]] = {}
    rows = 0
    with open(path, encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_no} — quatre colonnes attendues, {len(fields)}")
            try:
                start, end, value = int(fields[1]), int(fields[2]), float(fields[3])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no} — {exc}") from exc
            track.setdefault(fields[0], []).append((start, end, value))
            rows += 1
    return track, rows


@dataclass(frozen=True)
class DamID:
    """Le livrable de la semaine 7 : position radiale modélisée contre DamID mesuré."""

    correlation: float
    covered: int             # billes recouvertes par la piste
    total: int
    source: str

    def __str__(self) -> str:
        return (
            f"r = {self.correlation:+.3f} sur {self.covered:,}/{self.total:,} billes "
            f"couvertes · source {self.source}"
        )


def damid(
    radial: np.ndarray,
    starts: np.ndarray,
    ends: np.ndarray,
    copy_id: np.ndarray,
    labels: tuple[str, ...],
    track: dict[str, list[tuple[int, int, float]]],
    source: str,
) -> DamID:
    """Corrèle la position radiale moyenne à une piste DamID mesurée.

    `track` associe à chaque chromosome des intervalles `(début, fin, score)`.
    Chaque bille reçoit la moyenne des scores pondérée par le recouvrement ; les
    billes non couvertes — bras acrocentriques, trous d'assemblage — sont exclues
    plutôt que mises à zéro, un zéro étant une valeur et non une absence.

    **Cette fonction n'a encore jamais tourné sur des données réelles.** Aucune
    entrée DamID du manifeste n'a pu être récupérée (`DATA_SOURCES.md` § 8). Elle
    est écrite, testée sur une piste construite exprès, et attend son fichier.
    """
    chroms = [label.split(":", 1)[0] for label in labels]
    mean_r = radial.mean(axis=0)

    score = np.full(len(starts), np.nan)
    for k in range(len(starts)):
        spans = track.get(chroms[copy_id[k]])
        if not spans:
            continue
        s0, e0 = int(starts[k]), int(ends[k])
        weight = value = 0.0
        for s1, e1, v in spans:
            overlap = min(e0, e1) - max(s0, s1)
            if overlap > 0:
                weight += overlap
                value += overlap * v
        if weight > 0:
            score[k] = value / weight

    ok = np.isfinite(score)
    r = float(np.corrcoef(mean_r[ok], score[ok])[0, 1]) if ok.sum() > 2 else float("nan")
    return DamID(correlation=r, covered=int(ok.sum()), total=len(starts), source=source)
