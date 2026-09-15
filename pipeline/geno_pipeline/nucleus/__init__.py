"""Noyau diploïde complet : 46 chaînes de billes TAD, placées sans se traverser.

La semaine 5 reconstruisait *une* chaîne à partir de contacts. La semaine 6 pose
la question d'échelle : un noyau entier, 46 polymères, deux mètres d'ADN dans dix
micromètres, et rien qui se traverse.

Le découpage est dans `beads`, le caryotype dans `karyotype`, le solveur dans
`pack`, les mesures dans `metrics`, l'assemblage dans `build`.
"""

from .beads import Beads, capacity_at, segment
from .build import Nucleus, build, save
from .karyotype import Karyotype, chrom_sizes, gm12878
from .metrics import Periphery, Territoriality, gyration, periphery, radial, territoriality
from .pack import Quality, System, quality, relax

__all__ = [
    "Beads",
    "Karyotype",
    "Nucleus",
    "Periphery",
    "Quality",
    "System",
    "Territoriality",
    "build",
    "capacity_at",
    "chrom_sizes",
    "gm12878",
    "gyration",
    "periphery",
    "quality",
    "radial",
    "relax",
    "save",
    "segment",
    "territoriality",
]
