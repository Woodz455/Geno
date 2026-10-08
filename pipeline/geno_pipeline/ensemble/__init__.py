"""Ensembles : N repliements du même génome, et ce qu'on a le droit d'en dire.

Le principe n° 1 du projet tient en une phrase — *une structure Hi-C unique est
un artefact statistique* — et la semaine 7 est celle où il cesse d'être une
déclaration d'intention. On produit N structures, et on mesure ce qui, dans une
seule d'entre elles, était une propriété du génome et ce qui n'était qu'un tirage.

`generate` produit, `store` range en `.zarr`, `stats` mesure — et refuse de
mesurer ce qui n'a pas de sens sur un ensemble sans repère commun.
"""

from .generate import generate
from .stats import (
    REGIMES,
    Contacts,
    DamID,
    Frame,
    Reproducibility,
    contact_pairs,
    contacts,
    damid,
    frame,
    homolog_map,
    medoid,
    radial_by_quartile,
    read_bedgraph,
    reproducibility,
    slope,
    self_consistency,
)
from .store import Ensemble, create, pending, put, read

__all__ = [
    "REGIMES",
    "Contacts",
    "DamID",
    "Ensemble",
    "Frame",
    "Reproducibility",
    "contact_pairs",
    "contacts",
    "create",
    "damid",
    "frame",
    "generate",
    "homolog_map",
    "medoid",
    "pending",
    "put",
    "radial_by_quartile",
    "read",
    "read_bedgraph",
    "reproducibility",
    "self_consistency",
    "slope",
]
