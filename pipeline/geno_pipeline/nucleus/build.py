"""Construction d'un noyau diploïde complet, de bout en bout.

Trois étapes, et la première est la moins évidente.

**1. Placer les territoires.** Les 46 copies reçoivent d'abord un domaine sphérique,
placé par le même solveur que les billes — sphères inégales, volume exclu,
confinement, biais radial proportionnel au contenu LAD de la copie. Le volume
d'un territoire suit sa longueur génomique, mais à une densité *interne* plus
forte que la moyenne nucléaire : la chromatine ne remplit pas le noyau, il reste
entre les territoires l'espace interchromatinien.

**2. Faire pousser chaque chaîne dans son territoire.** Ce n'est pas une
commodité numérique, c'est une hypothèse biologique, et il faut la nommer : une
relaxation sous contraintes ne fait **jamais** se croiser deux chaînes. La
topologie du modèle est donc celle de son initialisation. Partir de chaînes
mélangées produirait un noyau enlacé qu'aucun recuit ne démêlerait ; partir de
chaînes séparées suppose des territoires. Ce qui tranche, c'est la mitose : les
chromosomes se décondensent là où ils se sont retrouvés en fin de télophase, ils
n'ont pas à se trier. L'initialisation territoriale est l'hypothèse fidèle.

Conséquence méthodologique, assumée : **la territorialité n'est pas un résultat
de ce modèle, c'est une entrée.** Ce que le modèle peut dire, et que le rapport
mesure, c'est si le recuit la préserve. Un recuit qui la ferait fondre invaliderait
le modèle ; un recuit qui la laisse intacte ne prouve rien sur la biologie.

**3. Recuire le noyau entier.** Volume exclu partout, chaînes tendues, confinement,
et la préférence LAD vers la périphérie. Le polissage final coupe la préférence :
l'absence d'interpénétration est une condition, la périphérie une envie.

Les billes voisines sont des sphères tangentes, et le volume exclu interdit aux
autres de passer entre elles — le tube est plein. Une chaîne ne peut donc pas se
traverser elle-même tant que les contraintes tiennent, et le nombre d'enlacements
est un invariant du recuit, pas un paramètre libre.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .beads import Beads, segment
from .karyotype import Karyotype, gm12878
from .metrics import Periphery, Territoriality, gyration, periphery, territoriality
from .pack import Quality, System, quality, relax, sealed


@dataclass(frozen=True)
class Nucleus:
    """Un noyau diploïde placé, et tout ce qui permet de le juger."""

    beads: Beads
    x: np.ndarray                      # (n, 3) en µm
    seed: int
    karyotype: str
    assembly: str
    provenance: str
    final: Quality
    before: Territoriality
    after: Territoriality
    rim: Periphery
    territory_radius: np.ndarray       # (46,) rayon du domaine attribué, µm
    sealed: bool                       # le tube est-il resté étanche pendant l'inflation
    trace: list[tuple[int, float, float]]

    @property
    def n(self) -> int:
        return self.beads.n


def _territories(
    beads: Beads,
    *,
    phi_local: float,
    slack: float,
    seed: int,
    steps: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Un domaine sphérique par copie : centre et rayon.

    Les centres sont placés en faisant s'exclure des sphères **plus petites** que
    les domaines (facteur `slack`), si bien que les domaines se recouvrent d'autant.
    Ce n'est pas un relâchement de confort : des territoires strictement disjoints
    n'existent pas — ils s'interdigitent en surface, et c'est ce recouvrement qui
    fait l'entremêlement mesuré entre chromosomes voisins. Les rendre disjoints
    rendrait aussi le pavage infaisable, 46 sphères inégales ne remplissant pas une
    sphère au-delà de l'empilement aléatoire.
    """
    n_copies = len(beads.labels)
    counts = np.bincount(beads.copy_id, minlength=n_copies)
    mean_r = np.array(
        [beads.radius[beads.copy_id == c].mean() for c in range(n_copies)]
    )
    radius = mean_r * np.cbrt(counts / phi_local)

    # Même logique de rang que pour les billes : la copie la plus riche en LAD
    # vise la périphérie, la plus pauvre le centre, et le profil de densité visé
    # reste celui d'une sphère uniforme.
    lad = np.array([beads.lad[beads.copy_id == c].mean() for c in range(n_copies)])
    rank = np.empty(n_copies)
    rank[np.argsort(lad)] = np.arange(n_copies)
    outward = np.cbrt((rank + 0.5) / n_copies)

    rng = np.random.default_rng(seed)
    x = rng.normal(0.0, 1.0, (n_copies, 3))
    x *= (beads.nuclear_radius * 0.55 * rng.random(n_copies) ** (1 / 3)
          / np.linalg.norm(x, axis=1))[:, None]

    sys = System(x=x, r=radius * slack, R=beads.nuclear_radius, outward=outward)
    relax(sys, steps=steps, t0=0.25, t1=0.002, outward_gain=0.05,
          polish=600, tol=0.02, seed=seed + 1)
    return sys.x, radius


def _cell(p: np.ndarray, step: float) -> tuple[int, int, int]:
    """Case de la grille de hachage. Division *plancher* : `int()` tronque vers zéro
    et collerait les deux moitiés de chaque axe dans la même case autour de l'origine."""
    return int(p[0] // step), int(p[1] // step), int(p[2] // step)


_NEIGHBOURHOOD = tuple(
    (a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)
)


def _grow(
    x: np.ndarray,
    radius: np.ndarray,
    grid: dict[tuple[int, int, int], list[int]],
    cell: float,
    order: np.ndarray,
    centre: np.ndarray,
    territory: float,
    nuclear: float,
    rng: np.random.Generator,
    persistence: float,
    trials: int = 10,
) -> None:
    """Fait pousser une chaîne auto-évitante dans son territoire, en place dans `x`.

    L'auto-évitement n'est pas un raffinement. Une marche persistante ordinaire de
    pas `2r` repose régulièrement une bille sur une précédente : mesuré, le
    chevauchement initial atteignait 99 % — deux billes confondues. Le recuit part
    alors d'une configuration que le volume exclu doit défaire avant de pouvoir
    faire quoi que ce soit d'utile, et il n'y arrive jamais complètement.

    La grille est **commune aux 46 copies** et les copies poussent de la plus
    longue à la plus courte. Une grille par copie ne servirait à rien : les
    territoires se recouvrent, donc les collisions qui comptent sont justement
    celles entre chaînes différentes. Mesuré avec une grille par copie, le
    chevauchement initial restait à 99 %.

    À chaque pas on tire jusqu'à `trials` directions et on garde la première qui
    laisse la place ; si aucune ne convient, la moins mauvaise. Chaque essai ne
    regarde que 27 cellules.
    """

    def place(i: int, p: np.ndarray) -> None:
        x[i] = p
        grid.setdefault(_cell(p, cell), []).append(i)

    def room(p: np.ndarray, r: float, skip: int) -> float:
        """Distance à la bille placée la plus proche, en fraction du contact.

        ≥ 1 = la place est libre. `skip` est le prédécesseur immédiat, qui est à
        exactement un pas : le compter ferait échouer le test à l'epsilon flottant
        près, et **aucune** direction ne serait jamais acceptée.
        """
        kx, ky, kz = _cell(p, cell)
        worst = np.inf
        for a, b, c in _NEIGHBOURHOOD:
            for idx in grid.get((kx + a, ky + b, kz + c), ()):
                if idx == skip:
                    continue
                delta = p - x[idx]
                worst = min(worst, float(delta @ delta) / (r + radius[idx]) ** 2)
                if worst < 1.0:
                    return worst
        return worst

    def confine(p: np.ndarray, previous: np.ndarray, d: np.ndarray, r: float):
        """Réflexion spéculaire sur la paroi du territoire, puis sur celle du noyau."""
        for anchor, limit in ((centre, territory), (np.zeros(3), nuclear - r)):
            off = p - anchor
            far = float(np.linalg.norm(off))
            if far <= limit:
                continue
            normal = off / far
            d = d - 2.0 * float(d @ normal) * normal
            d /= np.linalg.norm(d)
            p = previous + float(np.linalg.norm(p - previous)) * d
            off = p - anchor
            far = float(np.linalg.norm(off))
            if far > limit:
                p = anchor + off / far * limit * 0.97
        return p, d

    step = 2.0 * float(radius[order].mean())
    first = int(order[0])
    place(first, centre + rng.normal(0.0, 1.0, 3) * territory * 0.25)
    d = rng.normal(0.0, 1.0, 3)
    d /= np.linalg.norm(d)

    for i in order[1:]:
        i = int(i)
        previous = x[i - 1]
        best_p, best_room, best_d = None, -np.inf, d
        for t in range(trials):
            # Plus les essais échouent, plus on s'autorise à tourner court.
            keep = persistence * (1.0 - t / trials)
            nd = keep * d + (1.0 - keep) * rng.normal(0.0, 1.0, 3)
            nd /= np.linalg.norm(nd)
            candidate, nd = confine(previous + step * nd, previous, nd, float(radius[i]))

            free = room(candidate, float(radius[i]), i - 1)
            if free > best_room:
                best_p, best_room, best_d = candidate, free, nd
            if free >= 1.0:
                break

        d = best_d
        place(i, best_p)


def build(
    karyotype: Karyotype | None = None,
    *,
    bp_per_bead: int = 750_000,
    nuclear_radius: float = 5.0,
    phi: float = 0.30,
    phi_local: float = 0.55,
    territory_slack: float = 0.85,
    persistence: float = 0.72,
    bond_stretch: float = 1.15,
    outward_gain: float = 0.05,
    inflate_from: float = 0.65,
    polish_t: float = 0.02,
    bond_gain: float = 0.5,
    steps: int = 600,
    polish: int = 4_000,
    territory_steps: int = 700,
    tol: float = 0.01,
    seed: int = 0,
    lad_seed: int | None = None,
    trace_every: int = 0,
) -> Nucleus:
    # `seed` tire la conformation ; `lad_seed` tire le génome — le découpage en
    # billes et la piste LAD. Les séparer n'a d'intérêt que pour un ensemble, et
    # là c'est indispensable : faire varier la seule graine ferait varier la piste
    # LAD d'une structure à l'autre, et on mesurerait la variabilité de **deux**
    # génomes différents au lieu de la variabilité de repliement d'un seul. Par
    # défaut les deux coïncident, ce qui garde `build(seed=k)` reproductible.
    karyotype = karyotype or gm12878()
    beads = segment(
        karyotype,
        bp_per_bead=bp_per_bead,
        nuclear_radius=nuclear_radius,
        phi=phi,
        seed=seed if lad_seed is None else lad_seed,
    )

    centres, t_radius = _territories(
        beads,
        phi_local=phi_local,
        slack=territory_slack,
        seed=seed,
        steps=territory_steps,
    )

    rng = np.random.default_rng(seed + 977)
    x = np.empty((beads.n, 3))
    grid: dict[tuple[int, int, int], list[int]] = {}
    cell = 2.0 * float(beads.radius.max())

    counts = np.bincount(beads.copy_id, minlength=len(beads.labels))
    # La plus longue copie passe en premier : elle a le moins de latitude, et la
    # place qui reste s'accommode plus facilement d'un petit chromosome.
    for c in np.argsort(-counts):
        order = np.flatnonzero(beads.copy_id == c)
        _grow(
            x,
            beads.radius,
            grid,
            cell,
            order,
            centres[c],
            max(t_radius[c] - beads.radius[order].mean(), cell),
            nuclear_radius,
            rng,
            persistence,
        )

    before = territoriality(x, beads.copy_id)

    sys = System(
        x=x,
        r=beads.radius,
        R=nuclear_radius,
        bonds=beads.bonds,
        bond_len=beads.bond_limits(bond_stretch),
        outward=beads.radial_targets(seed),
    )
    # `radial_targets` prend la graine de *conformation*, pas celle du génome :
    # départager deux billes de même contenu LAD est une indétermination du
    # placement, pas une propriété de la séquence. L'ensemble doit donc en
    # explorer les deux issues.
    final, trace = relax(
        sys,
        steps=steps,
        polish=polish,
        tol=tol,
        inflate_from=inflate_from,
        polish_t=polish_t,
        bond_gain=bond_gain,
        outward_gain=outward_gain,
        seed=seed + 31,
        trace_every=trace_every,
    )

    capacity, _ = beads.lamina_capacity()
    return Nucleus(
        beads=beads,
        x=sys.x,
        seed=seed,
        karyotype=karyotype.name,
        assembly=karyotype.assembly,
        provenance=karyotype.provenance,
        final=final,
        before=before,
        after=territoriality(sys.x, beads.copy_id),
        rim=periphery(sys.x, beads.radius, beads.lad, nuclear_radius, capacity),
        territory_radius=t_radius,
        sealed=sealed(inflate_from, bond_stretch),
        trace=trace,
    )


def save(nucleus: Nucleus, path: Path) -> Path:
    """Écrit la structure et, à côté, tout ce qu'il faut pour ne pas s'y tromper.

    Le format définitif est le `.g3d` de la semaine 9. Ici un `.npz` suffit, à
    condition qu'il soit inséparable de sa provenance — d'où le `.json` frère qui
    dit d'où viennent les longueurs, d'où vient la piste LAD, et avec quels
    paramètres la configuration a été obtenue.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    b = nucleus.beads

    np.savez_compressed(
        path,
        x=nucleus.x.astype(np.float32),
        radius=b.radius.astype(np.float32),
        copy_id=b.copy_id,
        start=b.start,
        end=b.end,
        lad=b.lad.astype(np.float32),
        acrocentric=b.acrocentric,
    )

    sidecar = path.with_suffix(".json")
    sidecar.write_text(
        json.dumps(
            {
                "karyotype": nucleus.karyotype,
                "assembly": nucleus.assembly,
                "chrom_sizes_provenance": nucleus.provenance,
                "lad_source": b.lad_source,
                "evidence": "simulated",
                "labels": list(b.labels),
                "n_beads": b.n,
                "bp_per_bead": int(np.median(b.end - b.start)),
                "nuclear_radius_um": b.nuclear_radius,
                "phi": b.phi,
                "mean_bead_radius_nm": round(float(b.radius.mean()) * 1000, 1),
                "seed": nucleus.seed,
                "max_overlap": nucleus.final.max_overlap,
                "outside": nucleus.final.outside,
                "bond_stretch": nucleus.final.bond_stretch,
                "territory_index": nucleus.after.index,
                "gyration_um": [round(float(v), 3) for v in gyration(nucleus.x, b.copy_id)],
                "warning": (
                    "Structure simulée. Les positions ne sont contraintes par aucune "
                    "donnée de conformation mesurée ; la piste LAD est synthétique."
                ),
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    return sidecar
