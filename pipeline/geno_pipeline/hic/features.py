"""Extraction des features de conformation depuis une matrice de contacts.

On n'écrit pas ces algorithmes nous-mêmes : `cooler` et `cooltools` sont le
standard de fait du domaine (consortium Open2C), et réimplémenter ICE ou
l'insulation score serait se donner des bugs que personne d'autre n'a. Ce module
est la couche mince qui les enchaîne, fixe nos conventions et rend des tables
stables pour la suite du pipeline.

Enchaînement : équilibrage → compartiments → frontières de TAD → boucles.
Chaque étape suppose la précédente faite.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


def balance(clr, *, ignore_diags: int = 2, store: bool = True, **kw) -> np.ndarray:
    """Équilibrage ICE : retire les biais de couverture par bin.

    Le Hi-C brut porte des biais par bin — accessibilité, contenu GC, longueur
    de fragment, mappabilité. ICE les absorbe dans un poids par bin en imposant
    des marginales uniformes. Sans ça, un bin très couvert ressemble à un bin
    très connecté, et le vecteur propre décrit la couverture au lieu du
    compartiment.
    """
    import cooler

    weights, _stats = cooler.balance_cooler(
        clr, ignore_diags=ignore_diags, store=store, **kw
    )
    return weights


def compartments(clr, phasing_track: pd.DataFrame | None = None, **kw) -> pd.DataFrame:
    """Compartiments A/B par décomposition en vecteurs propres.

    **Le signe de E1 est arbitraire.** La décomposition rend une partition, pas
    une étiquette : rien dans la matrice ne dit quel côté est « actif ». Il faut
    une information extérieure — contenu GC ou densité génique — pour orienter.
    C'est le rôle de `phasing_track`. Sans elle, on rend la partition brute et
    l'appelant ne doit pas interpréter le signe.
    """
    import cooltools

    view = _whole_chrom_view(clr)
    _eigvals, eigvecs = cooltools.eigs_cis(
        clr, phasing_track=phasing_track, view_df=view, n_eigs=3, **kw
    )
    out = eigvecs[["chrom", "start", "end", "E1"]].copy()
    out["compartment"] = np.where(out["E1"] > 0, "A", "B")
    out.loc[out["E1"].isna(), "compartment"] = None
    out.attrs["phased"] = phasing_track is not None
    return out


def tad_boundaries(clr, window_bp: int, *, threshold: str = "Li", **kw) -> pd.DataFrame:
    """Frontières de TAD par insulation score.

    Le score mesure, pour chaque bin, le flux de contacts qui traverse une
    fenêtre glissante centrée sur lui. Un minimum local = peu de contacts
    traversants = une frontière.

    **La fenêtre doit être nettement plus petite que le TAD cherché**, pas de sa
    taille. Mesuré sur structure plantée (TADs de 490 kb en médiane), le rappel
    passe de 92 % à 100–150 kb à 58 % à 600 kb : une fenêtre trop large enjambe
    la frontière et lisse le minimum qu'on cherche. Voir `make hic-validate`.
    """
    import cooltools

    table = cooltools.insulation(
        clr, [window_bp], view_df=_whole_chrom_view(clr), threshold=threshold, **kw
    )
    score_col = f"log2_insulation_score_{window_bp}"
    flag_col = f"is_boundary_{window_bp}"
    out = table[["chrom", "start", "end", score_col, flag_col]].rename(
        columns={score_col: "insulation", flag_col: "is_boundary"}
    )
    out.attrs["window_bp"] = window_bp
    return out


def loops(
    clr,
    *,
    max_loci_separation: int = 2_000_000,
    lambda_bin_fdr: float = 0.1,
    clustering_radius: int | None = None,
    **kw,
) -> pd.DataFrame:
    """Boucles par détection de pics locaux (approche HiCCUPS).

    Un pixel est une boucle s'il dépasse nettement son voisinage local dans
    plusieurs directions — pas seulement la moyenne à cette distance. C'est ce
    qui distingue une vraie ancre d'une simple bande dense.
    """
    import cooltools

    view = _whole_chrom_view(clr)
    expected = cooltools.expected_cis(clr, view_df=view, nproc=1)
    if clustering_radius is None:
        clustering_radius = 2 * clr.binsize
    return cooltools.dots(
        clr,
        expected,
        view_df=view,
        max_loci_separation=max_loci_separation,
        lambda_bin_fdr=lambda_bin_fdr,
        clustering_radius=clustering_radius,
        nproc=1,
        **kw,
    )


def _whole_chrom_view(clr) -> pd.DataFrame:
    """Une région par chromosome — la vue par défaut de cooltools."""
    return pd.DataFrame(
        {
            "chrom": list(clr.chromnames),
            "start": 0,
            "end": [clr.chromsizes[c] for c in clr.chromnames],
            "name": list(clr.chromnames),
        }
    )


# --------------------------------------------------------------------------
# Évaluation contre une vérité terrain plantée
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Agreement:
    """Accord entre un appel et la vérité. Une seule façon de compter, partout."""

    recall: float
    precision: float
    n_true: int
    n_called: int
    n_matched: int

    @property
    def f1(self) -> float:
        if self.recall + self.precision == 0:
            return 0.0
        return 2 * self.recall * self.precision / (self.recall + self.precision)

    def __str__(self) -> str:
        return (
            f"rappel {self.recall:.0%}  précision {self.precision:.0%}  "
            f"F1 {self.f1:.0%}  ({self.n_matched}/{self.n_true} trouvés, "
            f"{self.n_called} appelés)"
        )


def match_positions(called: np.ndarray, true: np.ndarray, tol: int = 1) -> Agreement:
    """Apparie deux jeux de positions 1D à ± tol bins."""
    called = np.unique(np.asarray(called, dtype=np.int64))
    true = np.unique(np.asarray(true, dtype=np.int64))
    if true.size == 0:
        return Agreement(0.0, 0.0, 0, called.size, 0)

    hit_true = sum(1 for t in true if np.any(np.abs(called - t) <= tol)) if called.size else 0
    hit_called = sum(1 for c in called if np.any(np.abs(true - c) <= tol)) if called.size else 0
    return Agreement(
        recall=hit_true / true.size,
        precision=hit_called / called.size if called.size else 0.0,
        n_true=int(true.size),
        n_called=int(called.size),
        n_matched=int(hit_true),
    )


def match_pairs(called: np.ndarray, true: np.ndarray, tol: int = 2) -> Agreement:
    """Apparie deux jeux de paires 2D (boucles) à ± tol bins sur chaque axe."""
    called = np.atleast_2d(np.asarray(called, dtype=np.int64))
    true = np.atleast_2d(np.asarray(true, dtype=np.int64))
    if true.size == 0:
        return Agreement(0.0, 0.0, 0, len(called), 0)
    if called.size == 0:
        return Agreement(0.0, 0.0, len(true), 0, 0)

    close = (np.abs(called[:, None, :] - true[None, :, :]) <= tol).all(axis=2)
    return Agreement(
        recall=float(close.any(axis=0).mean()),
        precision=float(close.any(axis=1).mean()),
        n_true=len(true),
        n_called=len(called),
        n_matched=int(close.any(axis=0).sum()),
    )
