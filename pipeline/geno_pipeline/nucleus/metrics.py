"""Ce qu'on mesure sur un noyau une fois qu'il est placé.

Trois familles, qui ne se valent pas et qu'il ne faut pas confondre :

- **Conditions** — interpénétration, confinement, liaisons. Elles sont imposées
  par le solveur ; les mesurer vérifie le solveur, pas la biologie.
- **Propriétés imposées** — la stratification radiale des LADs. Elle est demandée
  par un terme du modèle ; la mesurer dit dans quelle mesure la demande a été
  *satisfaite*, ce qui est intéressant seulement parce qu'elle peut ne pas l'être.
- **Propriétés héritées** — la territorialité. Elle vient de l'initialisation, pas
  du recuit, parce qu'une relaxation ne fait jamais se croiser deux chaînes. La
  mesurer avant *et* après dit si le recuit la préserve, ce qui est la seule
  question honnête qu'on puisse lui poser.

Rien ici n'est une validation contre des données. La validation contre le DamID
publié est le livrable de la semaine 7, et il faut le réseau.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree


@dataclass(frozen=True)
class Territoriality:
    """À quel point une bille est entourée de billes de son propre chromosome."""

    same: float            # fraction moyenne de voisins de la même copie
    expected: float        # la même chose si tout était mélangé au hasard
    index: float           # same / expected — 1 = mélangé, >> 1 = territoires
    intermingling: float   # 1 - same
    k: int

    def __str__(self) -> str:
        return (
            f"voisins de la même copie {self.same:.1%} (hasard {self.expected:.1%}) "
            f"→ indice {self.index:.1f}× · entremêlement {self.intermingling:.1%}"
        )


def territoriality(x: np.ndarray, copy_id: np.ndarray, k: int = 20) -> Territoriality:
    """Fraction des `k` plus proches voisins appartenant à la même copie.

    Les voisins de chaîne immédiats (±2 indices) sont exclus : ils sont voisins
    par construction — ils sont *liés* — et les compter mesurerait la liaison, pas
    le territoire. Aux 46 jointures entre copies, cette exclusion retire aussi
    quelques billes d'une autre copie ; sur 8 000 billes l'effet est sous le
    dixième de pour-cent.
    """
    n = len(x)
    k = min(k, n - 6)
    _, nbr = cKDTree(x).query(x, k=k + 6)
    rows = np.arange(n)[:, None]

    eligible = np.abs(nbr - rows) > 2
    rank = np.cumsum(eligible, axis=1) - 1
    take = eligible & (rank < k)
    kept = np.maximum(take.sum(axis=1), 1)

    same = ((copy_id[nbr] == copy_id[:, None]) & take).sum(axis=1) / kept

    counts = np.bincount(copy_id, minlength=int(copy_id.max()) + 1)
    expected = (counts[copy_id] - 1) / (n - 1)

    s, e = float(same.mean()), float(expected.mean())
    return Territoriality(same=s, expected=e, index=s / e, intermingling=1.0 - s, k=k)


@dataclass(frozen=True)
class Periphery:
    """Ce qui touche la lamina, ce qui pourrait y être, et ce qui y serait « normalement ».

    Trois fractions qu'il ne faut surtout pas confondre, et c'est le résultat de
    la semaine :

    - `share` — ce que le modèle obtient.
    - `uniform_share` — la part des billes qui serait dans la coquille de contact
      si la densité était uniforme. C'est la part de **volume** de cette coquille
      dans le volume accessible aux centres, rien d'autre, et elle vaut environ
      `3·gap·r/(R - r)`.
    - `capacity_share` — la borne d'empilement monocouche, `4·eta·(R-r)²/r²`
      rapportée au total. Elle est bien plus grande, mais l'atteindre suppose un
      noyau à croûte dense et intérieur creux.

    « 35 % du génome est en LAD » et « 35 % du génome touche la lamina » ne sont
    donc pas la même phrase, et la seconde ne découle pas de la première.
    """

    at_lamina: int             # billes en contact avec l'enveloppe
    share: float               # at_lamina / n
    uniform_share: float       # part de volume de la coquille de contact
    capacity_share: float      # borne monocouche, rapportée au total
    enrichment: float          # share / uniform_share
    lad_at_lamina: float       # part du contenu LAD qui touche, pondérée par la longueur
    lad_coverage: float        # part LAD du génome modélisé
    r_lad: float               # position radiale moyenne du quartile le plus LAD
    r_ilad: float              # … et du quartile le moins LAD
    correlation: float         # corrélation de Pearson entre fraction LAD et rayon

    def __str__(self) -> str:
        return (
            f"{self.at_lamina} billes à la lamina ({self.share:.1%}) · "
            f"densité uniforme {self.uniform_share:.1%} → enrichissement "
            f"{self.enrichment:.1f}× · borne d'empilement {self.capacity_share:.0%} · "
            f"{self.lad_at_lamina:.0%} du contenu LAD y touche · radial LAD {self.r_lad:.3f} "
            f"vs iLAD {self.r_ilad:.3f} · r = {self.correlation:+.3f}"
        )


def radial(x: np.ndarray, nuclear_radius: float) -> np.ndarray:
    """Position radiale normalisée : 0 au centre, 1 à l'enveloppe."""
    return np.linalg.norm(x, axis=1) / nuclear_radius


def periphery(
    x: np.ndarray,
    r: np.ndarray,
    lad: np.ndarray,
    nuclear_radius: float,
    capacity: int,
    *,
    gap: float = 0.5,
) -> Periphery:
    """Contact avec la lamina, et les deux façons de dire ce qu'on en attendait.

    « Contact » = l'écart entre la surface de la bille et l'enveloppe est
    inférieur à `gap` rayons. Cette bande est plus fine qu'un empilement de deux
    couches (~1,63 r entre couches voisines), donc elle ne contient qu'une
    monocouche — c'est ce qui rend la comparaison avec la borne légitime.
    """
    d = np.linalg.norm(x, axis=1)
    inner = nuclear_radius - (1.0 + gap) * r
    touching = d > inner

    # Le volume accessible à un *centre* de bille est la boule de rayon R - r,
    # pas R : le confinement porte sur la surface. Normaliser par R³ sous-estimait
    # la part uniforme de 10 % et faisait sortir un enrichissement de 1,10 sur des
    # points pourtant tirés uniformément.
    mean_r, R = float(r.mean()), nuclear_radius
    free = (R - mean_r) ** 3
    uniform = (free - (R - (1.0 + gap) * mean_r) ** 3) / free

    q = np.quantile(lad, [0.25, 0.75])
    rad = d / R
    span = lad.max() - lad.min()
    corr = float(np.corrcoef(lad, rad)[0, 1]) if span > 1e-12 else float("nan")

    weight = lad * (r**3)      # pondérer par la longueur génomique, pas par bille
    share = float(touching.mean())
    return Periphery(
        at_lamina=int(touching.sum()),
        share=share,
        uniform_share=uniform,
        capacity_share=min(capacity / len(x), 1.0),
        enrichment=share / uniform if uniform > 1e-12 else float("nan"),
        lad_at_lamina=float(weight[touching].sum() / weight.sum()),
        lad_coverage=float((lad * r**3).sum() / (r**3).sum()),
        r_lad=float(rad[lad >= q[1]].mean()),
        r_ilad=float(rad[lad <= q[0]].mean()),
        correlation=corr,
    )


def gyration(x: np.ndarray, copy_id: np.ndarray) -> np.ndarray:
    """Rayon de giration de chaque copie chromosomique."""
    out = np.empty(int(copy_id.max()) + 1)
    for c in range(len(out)):
        pts = x[copy_id == c]
        out[c] = np.sqrt(((pts - pts.mean(axis=0)) ** 2).sum(axis=1).mean())
    return out
