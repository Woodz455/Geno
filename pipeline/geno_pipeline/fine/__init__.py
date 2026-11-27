"""Échelle fine : sous le TAD, un polymère et des cohésines qui extrudent.

La semaine 7 a laissé le modèle de noyau avec un défaut chiffré : au-delà de
15 Mb sa P(s) s'aplatit à 0,0253 quand le Hi-C réel continue de décroître. Le
modèle a le *fait* des territoires, pas leur organisation interne. La semaine 8
descend d'un cran pour lui donner cette organisation là où elle naît — à
l'échelle de la boucle CTCF et du domaine, entre le kilobase et le mégabase.

Quatre étages, chacun vérifiable seul :

- `region` — la région, ses monomères, ses sites CTCF orientés et la vérité qu'on
  y plante (témoins compris) ;
- `extrusion` — le modèle 1D de Fudenberg 2016 : où sont les pieds des cohésines ;
- `polymer` — la dynamique de Langevin sous OpenMM, dont les liaisons suivent
  l'extrusion ;
- `observe` — carte de contacts, insulation, points d'angle, P(s), et le raccord
  `R(s)` avec le noyau entier.

Le rayon d'un monomère n'est **pas** un paramètre de cette semaine : il vient de
la loi de la semaine 6, ce qui rend le raccord falsifiable.
"""

from .extrusion import Extrusion, lef_count, lifetime_steps, simulate
from .observe import (
    HIC_REFERENCE,
    REGIMES,
    Dots,
    Junction,
    Map,
    call_boundaries,
    contact_map,
    dot_score,
    dots_by_kind,
    insulation,
    junction,
    ps,
    regime,
    separation_curve,
    stationarity,
)
from .polymer import Field, bond_catalogue, equilibration, start
from .polymer import simulate as fold
from .region import (
    CONVERGENT,
    DIVERGENT,
    KINDS,
    TANDEM,
    Region,
    Sites,
    Truth,
    blob_radius,
    genome_bp,
    monomer_radius,
    plant,
    read_ctcf,
)

__all__ = [
    "CONVERGENT",
    "DIVERGENT",
    "Dots",
    "Extrusion",
    "Field",
    "HIC_REFERENCE",
    "Junction",
    "KINDS",
    "Map",
    "REGIMES",
    "Stationarity",
    "Region",
    "Sites",
    "TANDEM",
    "Truth",
    "blob_radius",
    "bond_catalogue",
    "call_boundaries",
    "contact_map",
    "dot_score",
    "dots_by_kind",
    "equilibration",
    "fold",
    "genome_bp",
    "insulation",
    "junction",
    "lef_count",
    "lifetime_steps",
    "monomer_radius",
    "plant",
    "ps",
    "read_ctcf",
    "regime",
    "separation_curve",
    "stationarity",
    "simulate",
    "start",
]
