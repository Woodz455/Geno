"""Matrice Hi-C synthétique à structure **plantée**.

Pourquoi ce module existe alors qu'on veut des données réelles.

Le Hi-C réel n'a pas de vérité terrain. Les appels publiés de Rao 2014 sont un
*point de comparaison* entre deux méthodes, pas une vérité : quand notre caller
diverge, rien ne dit lequel des deux a raison. On ne peut donc pas prouver qu'un
caller est correct sur des données réelles — seulement qu'il est d'accord avec
un autre.

Sur du synthétique, on connaît la réponse parce qu'on l'a écrite. On plante des
compartiments, des TADs, des boucles et des biais de couverture, puis on vérifie
que chaque caller les retrouve. C'est la seule étape où « correct » a un sens.

L'ordre est donc : prouver la correction ici, puis comparer à Rao sur les vraies
données. Ce module ne remplace pas les données réelles, il les précède.

Modèle génératif, pour chaque paire de bins (i, j) :

    E[i,j] = profondeur
           × (|i-j| + 1)^(-alpha)          décroissance avec la distance, P(s)
           × (1 ± comp_gain)               compartiment : même type ou non
           × (1 + tad_gain)                si i et j sont dans le même TAD
           × (1 + bosse de boucle)         aux coins des TADs, là où CTCF agit
           × biais_i × biais_j             biais de couverture, ce qu'ICE doit retirer

    comptage[i,j] ~ Poisson(E[i,j])
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Truth:
    """Ce qui a été planté. Les tests comparent les appels à ça, rien d'autre."""

    chrom: str
    resolution: int
    n_bins: int
    compartment: np.ndarray   # +1 = A, -1 = B, par bin
    domain: np.ndarray        # identifiant de TAD, par bin
    boundaries: np.ndarray    # indices de bins où un TAD commence
    loops: np.ndarray         # (n, 2) paires de bins, coins de TAD
    bias: np.ndarray          # biais de couverture par bin
    gc: np.ndarray            # piste type GC, corrélée à A — sert à orienter E1

    @property
    def span(self) -> int:
        return self.n_bins * self.resolution


def _blocks(n_bins: int, lo: int, hi: int, rng: np.random.Generator) -> np.ndarray:
    """Débuts de blocs découpant [0, n_bins), longueurs tirées dans [lo, hi].

    `n_bins` n'est jamais un début de bloc : la fin du chromosome n'est pas une
    frontière. Planter une frontière fantôme là fabrique un vrai positif que
    personne ne peut appeler — il n'y a pas de bin derrière — et écrase le
    rappel mesuré sans que rien ne soit cassé dans le caller.
    """
    edges = [0]
    while edges[-1] + lo + lo <= n_bins:
        nxt = edges[-1] + int(rng.integers(lo, hi + 1))
        if nxt + lo > n_bins:      # le dernier bloc absorbe le reste
            break
        edges.append(nxt)
    return np.asarray(edges, dtype=np.int64)


def plant(
    n_bins: int = 2_000,
    resolution: int = 10_000,
    chrom: str = "chr1",
    seed: int = 0,
    *,
    alpha: float = 1.1,
    comp_gain: float = 0.35,
    tad_gain: float = 0.9,
    loop_gain: float = 5.0,
    bias_sd: float = 0.35,
    total_contacts: int = 6_000_000,
) -> tuple[pd.DataFrame, pd.DataFrame, Truth]:
    """Fabrique bins, pixels et la vérité terrain correspondante."""
    rng = np.random.default_rng(seed)

    # --- TADs d'abord : blocs de 200 à 800 kb --------------------------------
    tad_edges = _blocks(n_bins, 20, 80, rng)
    domain = np.empty(n_bins, dtype=np.int32)
    for k, lo in enumerate(tad_edges):
        hi = tad_edges[k + 1] if k + 1 < len(tad_edges) else n_bins
        domain[lo:hi] = k

    # --- compartiments ensuite, par groupes de TADs entiers -------------------
    #
    # Les frontières de compartiment sont un SOUS-ENSEMBLE des frontières de TAD,
    # et ce n'est pas une commodité : dans les données réelles, une transition
    # A/B tombe sur une frontière de TAD. Si on tire les deux indépendamment, un
    # changement de compartiment au milieu d'un TAD crée une vraie insulation que
    # le caller détecte correctement — mais qui n'est pas dans notre vérité. On
    # mesure alors comme une erreur du caller ce qui est une erreur du modèle.
    compartment = np.empty(n_bins, dtype=np.int8)
    sign, k = 1, 0
    while k < len(tad_edges):
        group = int(rng.integers(2, 6))               # 2 à 5 TADs par compartiment
        lo = tad_edges[k]
        hi = tad_edges[k + group] if k + group < len(tad_edges) else n_bins
        compartment[lo:hi] = sign
        sign, k = -sign, k + group

    # --- boucles : aux coins des TADs assez longs, là où CTCF forme une ancre -
    loops = []
    for k, lo in enumerate(tad_edges):
        hi = (tad_edges[k + 1] if k + 1 < len(tad_edges) else n_bins) - 1
        if hi - lo >= 25:
            loops.append((lo, hi))
    loops_arr = np.asarray(loops, dtype=np.int64)

    # --- biais de couverture : ce qu'ICE devra retirer -----------------------
    bias = np.exp(rng.normal(0.0, bias_sd, n_bins))

    # --- piste type GC : les compartiments A sont riches en GC ---------------
    gc = 0.42 + 0.05 * compartment + rng.normal(0, 0.012, n_bins)

    # --- matrice attendue ----------------------------------------------------
    idx = np.arange(n_bins)
    dist = np.abs(idx[:, None] - idx[None, :])
    expected = (dist + 1.0) ** (-alpha)

    same_comp = compartment[:, None] == compartment[None, :]
    expected *= np.where(same_comp, 1.0 + comp_gain, 1.0 - comp_gain)

    same_tad = domain[:, None] == domain[None, :]
    expected *= np.where(same_tad, 1.0 + tad_gain, 1.0)

    for a, b in loops_arr:                      # bosse gaussienne au coin
        lo_a, hi_a = max(0, a - 2), min(n_bins, a + 3)
        lo_b, hi_b = max(0, b - 2), min(n_bins, b + 3)
        da = (np.arange(lo_a, hi_a) - a)[:, None]
        db = (np.arange(lo_b, hi_b) - b)[None, :]
        bump = loop_gain * np.exp(-(da**2 + db**2) / 2.0)
        expected[lo_a:hi_a, lo_b:hi_b] *= 1.0 + bump
        expected[lo_b:hi_b, lo_a:hi_a] *= 1.0 + bump.T

    expected *= bias[:, None] * bias[None, :]
    expected *= total_contacts / expected.sum()

    # --- échantillonnage de Poisson sur le triangle supérieur ----------------
    i, j = np.triu_indices(n_bins)
    counts = rng.poisson(expected[i, j]).astype(np.int32)
    keep = counts > 0

    bins = pd.DataFrame(
        {
            "chrom": pd.Categorical([chrom] * n_bins, categories=[chrom]),
            "start": idx * resolution,
            "end": (idx + 1) * resolution,
        }
    )
    pixels = pd.DataFrame(
        {"bin1_id": i[keep], "bin2_id": j[keep], "count": counts[keep]}
    ).sort_values(["bin1_id", "bin2_id"], ignore_index=True)

    truth = Truth(
        chrom=chrom,
        resolution=resolution,
        n_bins=n_bins,
        compartment=compartment,
        domain=domain,
        boundaries=np.asarray(tad_edges[1:], dtype=np.int64),  # 0 n'est pas une frontière
        loops=loops_arr,
        bias=bias,
        gc=gc,
        )
    return bins, pixels, truth


def write_cool(path: Path, bins: pd.DataFrame, pixels: pd.DataFrame) -> Path:
    """Écrit un vrai fichier .cool — même format que les données 4DN."""
    import cooler

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cooler.create_cooler(
        str(path),
        bins=bins,
        pixels=pixels,
        dtypes={"count": "int32"},
        assembly="synthetic",
        ordered=True,
        symmetric_upper=True,
    )
    return path


def gc_track(truth: Truth) -> pd.DataFrame:
    """Piste de phasage pour orienter le vecteur propre.

    Le signe de E1 est mathématiquement arbitraire : seule une information
    extérieure décide quel côté est A. En vrai c'est le GC ou la densité génique.
    """
    return pd.DataFrame(
        {
            "chrom": [truth.chrom] * truth.n_bins,
            "start": np.arange(truth.n_bins) * truth.resolution,
            "end": (np.arange(truth.n_bins) + 1) * truth.resolution,
            "GC": truth.gc,
        }
    )
