"""Découpage du génome en billes d'échelle TAD, et ce que chaque bille porte.

Une bille = un TAD. Reste à dire *quel* TAD : Dixon 2012 en compte ~2 200 sur le
génome haploïde, de taille moyenne ~880 kb ; Rao 2014 en compte 9 274, médiane
185 kb. Ce ne sont pas deux mesures du même objet à quatre près, ce sont deux
définitions. À 750 kb par bille on est à l'échelle Dixon, et le noyau diploïde
tient en ~8 000 billes ; à l'échelle Rao il en faudrait ~37 000. La feuille de
route demande 6 000–10 000 pour la semaine 6, donc **échelle Dixon**, et le
rapport le dit plutôt que de laisser croire qu'« un TAD » est une unité.

Rayons. La chromatine a une densité locale à peu près constante, donc le volume
d'une bille suit sa longueur en paires de bases et son rayon suit la racine
cubique : `r ∝ L^(1/3)`. La constante est fixée par une seule quantité physique,
la **fraction volumique** `phi` occupée par la chromatine dans le noyau — mesurée
entre 12 % et 52 % selon les régions par ChromEMT (Ou 2017), ~30 % en moyenne.

Conséquence à retenir : une fois `phi` et le nombre de billes fixés, `r/R` est
déterminé. Le rayon nucléaire `R` n'est donc **qu'une unité** — il change les
micromètres affichés, rien d'autre dans la géométrie.

LADs. À 750 kb par bille, un LAD médian (~500 kb) est plus petit qu'une bille :
un drapeau booléen par bille serait un mensonge de résolution. Chaque bille porte
donc une **fraction LAD** continue dans [0, 1], et l'attraction vers la lamina en
est proportionnelle. La piste elle-même est pour l'instant synthétique — le
réseau est fermé — et c'est estampillé partout où ça sort.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .karyotype import Karyotype

PACKING_2D = np.pi / np.sqrt(12.0)   # 0,9069 — densité hexagonale, borne 2D


@dataclass(frozen=True)
class Beads:
    """Le génome diploïde découpé en billes, prêt à être placé dans l'espace."""

    copy_id: np.ndarray       # (n,) index de la copie chromosomique
    start: np.ndarray         # (n,) début sur le chromosome, 0-based
    end: np.ndarray           # (n,) fin, exclue
    radius: np.ndarray        # (n,) rayon en µm
    lad: np.ndarray           # (n,) fraction LAD dans [0, 1]
    labels: tuple[str, ...]   # étiquette par copie, indexée par copy_id
    acrocentric: np.ndarray   # (n,) bille sur un bras court acrocentrique
    nuclear_radius: float     # µm
    phi: float                # fraction volumique de chromatine
    lad_source: str           # "synthetic" ou le chemin de la piste réelle

    @property
    def n(self) -> int:
        return len(self.start)

    @property
    def bonds(self) -> np.ndarray:
        """Paires (i, i+1) consécutives *à l'intérieur* d'une même copie.

        Le dernier bille d'un chromosome n'est pas liée à la première du suivant :
        c'est tout l'intérêt d'un noyau à 46 polymères plutôt qu'un seul.
        """
        same = self.copy_id[:-1] == self.copy_id[1:]
        i = np.flatnonzero(same)
        return np.c_[i, i + 1]

    def bond_limits(self, stretch: float = 1.15) -> np.ndarray:
        """Longueur **maximale** d'une liaison : au-delà, la chaîne casserait.

        Pas une longueur visée. Deux billes liées vivent entre `r_i + r_j`, que le
        volume exclu impose, et `stretch × (r_i + r_j)`, que la liaison impose ;
        entre les deux, la distance est libre. Un lien de chromatine s'étire, il
        ne fixe pas une distance, et une contrainte bilatérale se battrait contre
        le volume exclu dans les replis serrés.
        """
        b = self.bonds
        return stretch * (self.radius[b[:, 0]] + self.radius[b[:, 1]])

    def radial_targets(self, seed: int = 0) -> np.ndarray:
        """Rayon visé par bille, en fraction de `R - r`, ordonné par contenu LAD.

        Le rang compte, pas la valeur. Prendre la fraction LAD elle-même comme
        rayon visé placerait la bille moyenne — 34 % de LAD — à 0,34 R, donc
        **comprimerait** tout le noyau vers son centre : mesuré, le rayon moyen
        passait de 0,67 à 0,56 et plus aucune bille n'atteignait la lamina. Le
        biais doit trier, pas tasser.

        D'où la cible par quantile : la bille de rang `q` vise `q^(1/3)`, qui est
        exactement la loi des rayons d'une sphère de densité uniforme. Le profil
        de densité visé est donc celui qu'on aurait sans biais du tout, et seule
        l'attribution des rayons change. Les ex æquo — nombreux, une bille de
        750 kb sans aucun LAD est fréquente — sont départagés au hasard, ce qui
        est le traitement correct : entre deux billes de même contenu LAD, le
        modèle n'a pas de raison d'en préférer une.
        """
        rng = np.random.default_rng(seed)
        order = np.lexsort((rng.random(self.n), self.lad))
        rank = np.empty(self.n)
        rank[order] = np.arange(self.n)
        return np.cbrt((rank + 0.5) / self.n)

    def lamina_capacity(self) -> tuple[int, float]:
        """Combien de billes peuvent toucher l'enveloppe **en même temps**.

        Borne géométrique, indépendante du modèle : les billes au contact ont leur
        centre sur une sphère de rayon `R - r`, et leurs disques de rayon `r` y
        pavent une surface au mieux à la densité hexagonale.

            N1 = eta · 4·pi·(R - r)² / (pi · r²)

        Rapportée au total `N = phi · (R/r)³`, la fraction du génome qui peut être
        à la lamina vaut `4·eta·r·(R - r)² / (phi·R³)` — **linéaire en r**. Un
        modèle plus fin a proportionnellement moins de place à la périphérie.
        """
        r = float(self.radius.mean())
        R = self.nuclear_radius
        n1 = int(PACKING_2D * 4.0 * (R - r) ** 2 / r**2)
        return n1, n1 / self.n


def capacity_at(
    total_bp: int,
    bp_per_bead: int,
    *,
    nuclear_radius: float = 5.0,
    phi: float = 0.30,
) -> tuple[int, float, float]:
    """Borne de périphérie à une résolution donnée, sans construire le modèle.

    Renvoie `(nombre de billes, rayon en nm, fraction du génome qui peut toucher)`.
    Utile parce que la réponse ne dépend que de la résolution : `f = 4·eta·r/(phi·R)`
    au premier ordre. Un modèle deux fois plus fin a deux fois moins de place à la
    lamina — donc « fraction des LADs à la lamina » n'est pas une quantité
    comparable entre deux modèles de granularité différente.
    """
    n = max(1, int(round(total_bp / bp_per_bead)))
    r = nuclear_radius * (phi / n) ** (1 / 3)
    n1 = PACKING_2D * 4.0 * (nuclear_radius - r) ** 2 / r**2
    return n, r * 1000.0, min(n1 / n, 1.0)


def _lad_fraction(
    length: int,
    edges: np.ndarray,
    rng: np.random.Generator,
    *,
    resolution: int = 50_000,
    mean_lad: float = 1_200_000,
    mean_ilad: float = 2_200_000,
) -> np.ndarray:
    """Piste LAD synthétique : chaîne de Markov à deux états, moyennée par bille.

    Les longueurs de séjour reproduisent ce que DamID rapporte — LADs de ~0,1 à
    10 Mb, médiane de l'ordre du demi-mégabase, couverture ~35 % du génome — et
    surtout leur **autocorrélation** : un LAD est un bloc, pas un tirage bille par
    bille. C'est cette autocorrélation qui décide si la contrainte de périphérie
    est satisfaisable, parce qu'elle décide si les billes LAD d'une même copie
    peuvent atteindre la lamina ensemble.

    Reste que c'est une piste **fabriquée**. Ce qu'on peut en tirer, c'est la
    géométrie ; ce qu'on ne peut pas, c'est la corrélation avec le DamID publié.
    """
    n_steps = max(1, length // resolution)
    p_leave = resolution / mean_lad      # LAD → inter-LAD
    p_enter = resolution / mean_ilad     # inter-LAD → LAD

    state = np.empty(n_steps, dtype=bool)
    state[0] = rng.random() < mean_lad / (mean_lad + mean_ilad)
    draws = rng.random(n_steps)
    cur = state[0]
    for k in range(1, n_steps):
        cur = (not cur) if draws[k] < (p_leave if cur else p_enter) else cur
        state[k] = cur

    # Moyenne par bille. `edges` est en pb ; on convertit en pas de résolution.
    cuts = np.clip(edges // resolution, 0, n_steps)
    cum = np.concatenate([[0.0], np.cumsum(state, dtype=np.float64)])
    width = np.maximum(np.diff(cuts), 1)
    return (cum[cuts[1:]] - cum[cuts[:-1]]) / width


def segment(
    karyotype: Karyotype,
    *,
    bp_per_bead: int = 750_000,
    nuclear_radius: float = 5.0,
    phi: float = 0.30,
    seed: int = 0,
) -> Beads:
    """Découpe chaque copie en billes d'échelle TAD et leur donne un rayon.

    Le nombre de billes d'une copie est arrondi, donc leur taille exacte varie
    d'un chromosome à l'autre — jamais à l'intérieur d'un chromosome, pour que le
    long des chaînes les rayons soient homogènes.
    """
    rng = np.random.default_rng(seed)

    copy_id, starts, ends, lads, acro = [], [], [], [], []
    for ci, cp in enumerate(karyotype.copies):
        k = max(1, int(round(cp.length / bp_per_bead)))
        edges = np.linspace(0, cp.length, k + 1).round().astype(np.int64)
        copy_id.append(np.full(k, ci, dtype=np.int32))
        starts.append(edges[:-1])
        ends.append(edges[1:])
        lads.append(_lad_fraction(cp.length, edges, rng))
        # Le bras court acrocentrique fait grossièrement les 10 premiers Mb.
        acro.append((edges[:-1] < 10_000_000) & cp.acrocentric)

    copy_id = np.concatenate(copy_id)
    start = np.concatenate(starts)
    end = np.concatenate(ends)
    lad = np.concatenate(lads)
    acrocentric = np.concatenate(acro)

    # r ∝ L^(1/3), normalisé pour que le volume total des billes fasse
    # exactement `phi` fois le volume nucléaire.
    span = (end - start).astype(np.float64)
    shape = np.cbrt(span)
    scale = nuclear_radius * np.cbrt(phi / (shape**3).sum())
    radius = scale * shape

    return Beads(
        copy_id=copy_id,
        start=start,
        end=end,
        radius=radius,
        lad=lad,
        labels=tuple(c.label for c in karyotype.copies),
        acrocentric=acrocentric,
        nuclear_radius=nuclear_radius,
        phi=phi,
        lad_source="synthetic",
    )
