"""Recuit : place des sphères sans qu'elles se traversent, dans un noyau, sur des chaînes.

Trois contraintes, et un ordre qui n'est pas négociable :

1. **volume exclu** — deux billes ne se traversent pas. `d ≥ r_i + r_j`.
2. **connectivité** — la chaîne ne casse pas. `d ≤ stretch · (r_i + r_j)`, une
   longueur *maximale* et rien d'autre : entre le contact et cette limite, la
   distance est libre (§ `_tighten`).
3. **confinement** — tout tient dans le noyau. `|x| ≤ R - r`.

Plus un rappel radial facultatif vers un rayon visé — une force de rappel, pas une
poussée, et la distinction n'est pas cosmétique (§ `_radial_bias`).

La méthode est une descente sous contraintes avec bruit recuit : à chaque pas on
injecte un déplacement gaussien d'amplitude `T`, on projette les contraintes, et
`T` décroît géométriquement. C'est du recuit simulé au sens continu — le bruit
permet de sortir des configurations bloquées tôt, et sa décroissance fige la
solution. Ce n'est **pas** du Metropolis : aucun pas n'est rejeté. Dit autrement,
on ne prétend pas échantillonner une distribution de Boltzmann, seulement trouver
une configuration admissible ; c'est exactement ce dont la semaine 6 a besoin, et
la semaine 7 (ensembles) tirera sa variabilité de graines différentes, pas d'un
équilibre thermodynamique.

Les contraintes sont projetées « à la Jacobi » : chaque bille encaisse la
**moyenne** des corrections que ses paires lui demandent, pas leur somme. Une
bille prise entre dix voisines recevrait sinon dix fois la correction et
partirait à l'autre bout du noyau.

Deux points fixes guettent ce genre de solveur, tous deux mesurés ici et traités
dans `relax` : une projection de Jacobi à gain 1 qui stagne, et des cages
octaédriques où six billes en contact mutuel verrouillent trois liaisons.

Une limite structurelle à dire tout de suite : la relaxation ne fait jamais se
croiser deux chaînes. La topologie initiale est donc conservée — le nombre
d'enlacements entre chromosomes est fixé par l'initialisation, pas par le recuit.
C'est pour ça que l'initialisation est motivée biologiquement (§ `build`) et que
le rapport mesure la territorialité **avant et après** : si le recuit la
détruisait, le modèle serait faux ; s'il la crée de rien, c'est qu'on l'a
imposée.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.spatial import cKDTree

EMPTY_PAIRS = np.empty((0, 2), dtype=np.int64)


@dataclass
class System:
    """Des sphères, éventuellement enchaînées, dans une sphère."""

    x: np.ndarray                                  # (n, 3) centres, µm
    r: np.ndarray                                  # (n,) rayons, µm
    R: float                                       # rayon de confinement, µm
    bonds: np.ndarray = field(default_factory=lambda: EMPTY_PAIRS)
    bond_len: np.ndarray = field(default_factory=lambda: np.empty(0))   # longueur MAX
    outward: np.ndarray | None = None              # (n,) rayon visé, en fraction de (R - r)

    @property
    def n(self) -> int:
        return len(self.x)


@dataclass(frozen=True)
class Quality:
    """Ce qui décide si une configuration est admissible, en nombres."""

    max_overlap: float      # chevauchement max, en fraction de (r_i + r_j)
    mean_overlap: float     # moyenne sur les paires en contact
    n_overlapping: int
    outside: int            # billes dont la surface sort du noyau
    bond_stretch: float     # allongement max au-delà de la limite, relatif
    shakes: int = 0         # secousses locales qu'il a fallu pour y arriver

    def acceptable(self, tol: float = 0.01) -> bool:
        """Les trois conditions dures, à la même tolérance relative.

        La tension de liaison en fait partie : une chaîne tendue 20 % au-delà de
        sa limite viole une contrainte dure autant qu'une interpénétration, et
        c'est justement la signature des cages octaédriques (§ `relax`). Un
        critère qui ne regarderait que le chevauchement les laisserait passer.
        """
        return (
            self.max_overlap <= tol
            and self.bond_stretch <= tol
            and self.outside == 0
        )

    def __str__(self) -> str:
        shaken = f" · {self.shakes} secousse(s)" if self.shakes else ""
        return (
            f"chevauchement max {self.max_overlap:.4%} · {self.n_overlapping} paires · "
            f"{self.outside} hors noyau · liaison la plus tendue "
            f"+{self.bond_stretch:.2%}{shaken}"
        )


def _pairs(x: np.ndarray, r: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Paires qui se chevauchent, avec leur direction et leur profondeur.

    L'arbre cherche à `2·r_max`, ce qui est une sur-recherche d'autant plus faible
    que les rayons sont homogènes — ici ±0,5 %, les billes ayant toutes la même
    taille génomique à l'arrondi près.
    """
    pairs = cKDTree(x).query_pairs(2.0 * float(r.max()), output_type="ndarray")
    if len(pairs) == 0:
        return EMPTY_PAIRS, np.empty((0, 3)), np.empty(0)

    i, j = pairs[:, 0], pairs[:, 1]
    delta = x[i] - x[j]
    d = np.linalg.norm(delta, axis=1)
    contact = r[i] + r[j]
    hit = d < contact
    if not hit.any():
        return EMPTY_PAIRS, np.empty((0, 3)), np.empty(0)

    pairs, delta, d, contact = pairs[hit], delta[hit], d[hit], contact[hit]
    # Billes confondues : aucune direction à lire dans le vecteur nul.
    dead = d < 1e-12
    if dead.any():
        delta = delta.copy()
        delta[dead] = np.array([1.0, 0.0, 0.0])
        d = d.copy()
        d[dead] = 1.0
    return pairs, delta / d[:, None], (contact - d) / contact


def quality(sys: System) -> Quality:
    """Toujours mesurée à taille pleine : l'inflation est un moyen, pas un résultat."""
    pairs, _, rel = _pairs(sys.x, sys.r)
    outside = int((np.linalg.norm(sys.x, axis=1) + sys.r > sys.R + 1e-9).sum())

    stretch = 0.0
    if len(sys.bonds):
        b = sys.bonds
        d = np.linalg.norm(sys.x[b[:, 0]] - sys.x[b[:, 1]], axis=1)
        stretch = float(max((d / sys.bond_len).max() - 1.0, 0.0))

    return Quality(
        max_overlap=float(rel.max()) if len(rel) else 0.0,
        mean_overlap=float(rel.mean()) if len(rel) else 0.0,
        n_overlapping=len(pairs),
        outside=outside,
        bond_stretch=stretch,
    )


def _separate(sys: System, gain: float, scale: float = 1.0) -> int:
    r = sys.r * scale
    pairs, u, rel = _pairs(sys.x, r)
    if len(pairs) == 0:
        return 0
    i, j = pairs[:, 0], pairs[:, 1]
    depth = rel * (r[i] + r[j])                  # profondeur absolue, µm

    disp = np.zeros_like(sys.x)
    count = np.zeros(sys.n)
    half = 0.5 * depth[:, None] * u
    np.add.at(disp, i, half)
    np.add.at(disp, j, -half)
    np.add.at(count, i, 1.0)
    np.add.at(count, j, 1.0)

    sys.x += gain * disp / np.maximum(count, 1.0)[:, None]
    return len(pairs)


def _tighten(sys: System, gain: float) -> None:
    """Empêche une liaison de s'allonger au-delà de sa limite. **Rien de plus.**

    Une liaison à longueur imposée se bat contre le volume exclu : dans un repli
    serré, la liaison rapproche deux billes et en écrase une troisième, et comme
    la liaison gagne (mesuré : écart aux liaisons 0,3 %, chevauchement résiduel
    6 %), c'est le volume exclu qui cède. Or c'est lui la condition.

    Une contrainte **unilatérale** lève la frustration : la chaîne ne peut pas se
    rompre, les billes ne peuvent pas se traverser, et entre les deux la distance
    est libre. La chromatine ne dit pas autre chose — un lien s'étire, il ne fixe
    pas une distance.
    """
    if not len(sys.bonds):
        return
    i, j = sys.bonds[:, 0], sys.bonds[:, 1]
    delta = sys.x[i] - sys.x[j]
    d = np.maximum(np.linalg.norm(delta, axis=1), 1e-12)
    over = d > sys.bond_len
    if not over.any():
        return
    i, j, delta, d = i[over], j[over], delta[over], d[over]
    pull = (0.5 * gain * (d - sys.bond_len[over]) / d)[:, None] * delta
    np.add.at(sys.x, i, -pull)
    np.add.at(sys.x, j, pull)


def _confine(sys: System, scale: float = 1.0) -> None:
    d = np.linalg.norm(sys.x, axis=1)
    limit = sys.R - sys.r * scale
    over = d > limit
    if over.any():
        sys.x[over] *= (limit[over] / np.maximum(d[over], 1e-12))[:, None]


def _radial_bias(sys: System, gain: float) -> None:
    """Rappel vers un rayon cible, proportionnel au poids LAD de chaque bille.

    C'est une **force de rappel**, pas une poussée. La distinction n'est pas
    cosmétique et la semaine 5 l'avait déjà payée : une poussée vers l'extérieur
    s'accumule sur tous les pas du recuit, si bien que son effet dépend du nombre
    de pas et non du modèle. Mesuré ici, 500 pas d'une poussée même faible (gain
    0,03) plaquaient toutes les billes LAD contre l'enveloppe et y créaient une
    croûte bloquée dont le polissage ne sortait plus : chevauchement résiduel
    21 %, contre 5 % sans poussée du tout. Un rappel, lui, converge vers sa cible
    et peut ramener une bille vers l'intérieur.

    `outward` est le rayon visé en fraction de `R - r`. Ce n'est pas une position
    imposée — le volume exclu interdit à toutes les billes d'y arriver, et c'est
    justement ce que la semaine mesure.
    """
    if sys.outward is None:
        return
    d = np.maximum(np.linalg.norm(sys.x, axis=1), 1e-12)
    target = (sys.R - sys.r) * sys.outward
    sys.x += (gain * (target - d) / d)[:, None] * sys.x


def _offenders(sys: System, tol: float) -> np.ndarray:
    """Billes qui violent une condition dure au-delà de `tol`, chevauchement ou liaison.

    Exactement le critère de `Quality.acceptable`, pour que la secousse porte sur ce
    qui bloque et sur rien d'autre.
    """
    pairs, _, rel = _pairs(sys.x, sys.r)
    bad = pairs[rel > tol].ravel() if len(pairs) else np.empty(0, dtype=np.int64)
    if len(sys.bonds):
        d = np.linalg.norm(sys.x[sys.bonds[:, 0]] - sys.x[sys.bonds[:, 1]], axis=1)
        bad = np.concatenate([bad, sys.bonds[d > sys.bond_len * (1.0 + tol)].ravel()])
    return np.unique(bad)


def _shake(
    sys: System,
    rng: np.random.Generator,
    amplitude: float,
    who: np.ndarray,
    reach: float,
) -> int:
    """Secoue fort, mais **seulement** autour des fautifs et de leurs voisins.

    Un réchauffage global défait tout le noyau pour réparer six billes, et le
    paie : sur la graine 3, quatre réchauffages à 0,10 faisaient tomber la
    corrélation LAD-rayon de +0,665 à +0,596 sans même régler le problème.

    L'amplitude doit être de l'ordre du rayon d'une bille. En dessous, une cage
    de contacts ne s'ouvre pas : ce n'est pas un puits peu profond qu'on quitte
    par agitation, c'est une structure rigide qu'il faut casser.
    """
    if len(who) == 0:
        return 0
    around = cKDTree(sys.x).query_ball_point(sys.x[who], reach * float(sys.r.max()))
    idx = np.unique(np.concatenate([np.asarray(a, dtype=np.int64) for a in around]))
    sys.x[idx] += rng.normal(0.0, amplitude, (len(idx), 3)) * sys.r[idx, None]
    return len(idx)


def sealed(inflate_from: float, bond_stretch: float) -> bool:
    """Le tube de la chaîne reste-t-il étanche pendant toute l'inflation ?

    À l'échelle `s`, deux billes liées s'excluent au rayon `s·r` mais restent
    séparées d'au plus `bond_stretch·(r_i + r_j)`. Une troisième bille, de rayon
    `s·r` elle aussi, se glisse entre elles dès que l'écart le permet :

        bond_stretch · 2r − 2·s·r  ≥  2·s·r     ⟺     s ≤ bond_stretch / 2

    En dessous du seuil, une chaîne étrangère peut traverser une liaison. Mesuré à
    s = 0,50 avec stretch = 1,15 (seuil 0,575) : une liaison bloquée à 20 % au-delà
    de sa limite et 5 % de chevauchement résiduel, quand tout le reste était déjà à
    0,04 %.

    Ce n'est pas pour autant définitif, et il faut le dire : la même configuration,
    polie avec un plancher de bruit de 0,05 au lieu de 0,02 et un gain de liaison
    de 1,0 au lieu de 0,5, retombe à 0,35 %. Le passage crée un défaut *difficile*
    à défaire, pas un défaut indéfaisable. D'où le choix de rester au-dessus du
    seuil par défaut — c'est gratuit — plutôt que d'interdire en dessous.
    """
    return inflate_from > bond_stretch / 2.0


def relax(
    sys: System,
    *,
    steps: int = 1_500,
    t0: float = 0.45,
    t1: float = 0.004,
    inflate_from: float = 1.0,
    outward_gain: float = 0.05,
    separate_gain: float = 1.6,
    bond_gain: float = 0.5,
    polish: int = 800,
    polish_t: float = 0.02,
    shake_rounds: int = 6,
    shake: float = 0.8,
    shake_reach: float = 3.0,
    tol: float = 0.01,
    seed: int = 0,
    trace_every: int = 0,
) -> tuple[Quality, list[tuple[int, float, float]]]:
    """Recuit puis polissage. Renvoie la qualité finale et la trace éventuelle.

    **Inflation.** Placer d'emblée des billes à leur taille finale ne marche pas :
    une chaîne auto-évitante construite gloutonnement se piège, et on démarre à
    99 % de chevauchement — deux billes confondues — dont la relaxation ne sort
    plus. On démarre donc les rayons à `inflate_from` de leur valeur et on les
    fait croître pendant le recuit. C'est l'idée de Lubachevsky–Stillinger, et
    elle change la nature du problème : au lieu de *défaire* un empilement
    impossible, on en *fait croître* un possible.

    Les liaisons, elles, ne sont **pas** mises à l'échelle. L'inflation ne doit
    porter que sur le volume exclu ; contracter aussi les chaînes ferait
    s'effondrer les territoires que l'initialisation vient d'établir, et il
    faudrait ensuite les rouvrir contre des voisins déjà en place.

    `separate_gain` dépasse 1 : c'est de la sur-relaxation. Une projection de
    Jacobi à gain 1 stagne dès que les corrections demandées à une bille se
    compensent — mesuré, elle plafonnait à 21 % de chevauchement.

    Le polissage tourne à taille pleine et **sans biais radial** : la périphérie
    est une préférence, l'absence d'interpénétration est une condition. Quand les
    deux s'opposent, c'est la condition qui gagne — et le rapport mesure ensuite ce
    qu'il reste de la préférence, ce qui est précisément le résultat intéressant
    de la semaine.

    Il garde en revanche un **plancher de bruit** qui décroît linéairement jusqu'à
    zéro. À température strictement nulle, le polissage laissait le fond de la
    distribution s'installer : ~11 600 paires en chevauchement résiduel, contre
    **39** avec un plancher à 0,02. Une agitation résiduelle permet aux voisins de
    bouger ensemble, ce qu'une projection locale ne fait jamais. Le plancher
    atteint exactement zéro à la fin, donc la configuration rendue n'est jamais
    bruitée.

    Ce qu'il ne règle pas, ce sont les cages octaédriques traitées plus bas : elles
    survivaient aux deux réglages. Le bruit nettoie le fond, la secousse casse les
    structures rigides — deux problèmes, deux outils.
    """
    rng = np.random.default_rng(seed)
    trace: list[tuple[int, float, float]] = []
    clock = 0

    def anneal(n: int, hot: float, cold: float, inflate: bool, bias: float) -> None:
        nonlocal clock
        decay = (cold / hot) ** (1.0 / max(n - 1, 1))
        temperature = hot
        for step in range(n):
            scale = (
                inflate_from + (1.0 - inflate_from) * (step / max(n - 1, 1))
                if inflate
                else 1.0
            )
            sys.x += rng.normal(0.0, temperature, sys.x.shape) * sys.r[:, None]
            _radial_bias(sys, bias)
            _separate(sys, separate_gain, scale)
            _tighten(sys, bond_gain)
            _confine(sys, scale)
            temperature *= decay
            if trace_every and step % trace_every == 0:
                trace.append((clock + step, temperature, quality(sys).max_overlap))
        clock += n

    def polish_down(n: int) -> Quality:
        nonlocal clock
        for step in range(n):
            floor = polish_t * (1.0 - step / max(n - 1, 1))
            if floor > 0.0:
                sys.x += rng.normal(0.0, floor, sys.x.shape) * sys.r[:, None]
            _tighten(sys, bond_gain)
            _confine(sys)
            _separate(sys, separate_gain)
            _confine(sys)
            if step % 25 == 24:
                q = quality(sys)
                if trace_every:
                    trace.append((clock + step, floor, q.max_overlap))
                if q.acceptable(tol):
                    clock += step
                    return q
        clock += n
        return quality(sys)

    anneal(steps, t0, t1, inflate=True, bias=outward_gain)
    q = polish_down(polish)

    # Cages octaédriques. Le polissage a des **points fixes** : six billes en contact
    # mutuel forment un octaèdre rigide, et trois liaisons tombées sur ses trois
    # diagonales ne peuvent plus se raccourcir — raccourcir une diagonale écarterait
    # les quatre billes de l'équateur, ce que les deux autres liaisons interdisent.
    #
    # La signature est arithmétique et ne trompe pas : la diagonale d'un octaèdre
    # vaut √2 fois son arête, donc la liaison se fige à 1,3849 × contact pour une
    # arête à 0,98 × contact, soit **exactement +20,42 %** au-delà d'une limite à
    # 1,15. Ce nombre est ressorti à l'identique sur cinq configurations et trois
    # graines différentes — c'est √2, pas un hasard local.
    #
    # Une cage de contacts ne se quitte pas par agitation : il faut la casser. D'où
    # une secousse d'amplitude comparable au rayon d'une bille, appliquée aux seuls
    # fautifs et à leur voisinage, suivie d'un polissage.
    for attempt in range(1, shake_rounds + 1):
        if q.acceptable(tol):
            break
        _shake(sys, rng, shake, _offenders(sys, tol), shake_reach)
        q = replace(polish_down(polish // 4), shakes=attempt)

    return q, trace
