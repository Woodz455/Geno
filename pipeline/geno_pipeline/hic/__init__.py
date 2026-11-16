"""Conformation : matrices de contacts, features 3D, et leur validation.

`synthetic` plante une structure connue et l'écrit en .cool ; `features` extrait
compartiments, TADs et boucles via cooler/cooltools. Les deux ensemble donnent
la seule chose qu'on ne peut pas obtenir de données réelles : une preuve que les
callers retrouvent ce qui est réellement là.
"""

from .features import Agreement, balance, compartments, loops, match_pairs, match_positions, tad_boundaries
from .synthetic import Truth, gc_track, plant, write_cool

__all__ = [
    "Agreement",
    "Truth",
    "balance",
    "compartments",
    "gc_track",
    "loops",
    "match_pairs",
    "match_positions",
    "plant",
    "tad_boundaries",
    "write_cool",
]
