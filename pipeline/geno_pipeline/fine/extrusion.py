"""Extrusion de boucles, en une dimension : où sont les pieds des cohésines.

Le modèle est celui de Fudenberg 2016 et il tient en cinq règles :

1. un complexe (LEF) se charge au hasard sur deux monomères voisins ;
2. ses deux jambes s'écartent, une vers la gauche, une vers la droite ;
3. elles ne se dépassent pas entre complexes — un pied occupé bloque ;
4. un motif CTCF **qui leur fait face** les arrête, avec sa probabilité d'occupation ;
5. le complexe se décroche au bout d'une durée de vie exponentielle.

Rien d'autre. Aucune notion de TAD, de compartiment ni de boucle n'est écrite
ici : les domaines sont censés *sortir* de la règle 4, et c'est exactement ce que
la validation vérifie.

**Pourquoi un étage 1D séparé du 3D.** Parce qu'il est vérifiable seul et qu'il
coûte une milliseconde là où le 3D coûte des minutes. La règle de convergence est
une affirmation sur les *pieds*, pas sur les distances : si la carte d'ancrage 1D
ne piquait pas aux paires convergentes, aucune dynamique moléculaire ne le
rattraperait, et il vaut mieux l'apprendre tout de suite.

**Une unité de temps à ne pas surinterpréter.** Un pas d'extrusion = une jambe
avance d'un monomère. La durée de vie s'en déduit par la processivité :
`processivité = 2 · vitesse · durée_de_vie`, en paires de bases. Convertir ces
pas en secondes demanderait une vitesse mesurée de cohésine, qui n'entre nulle
part dans ce qu'on rapporte — on ne rapporte que des moyennes d'état
stationnaire.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .region import Region, Sites


@dataclass(frozen=True)
class Extrusion:
    """Les instantanés de pieds, et de quoi les lire.

    `legs` garde **un emplacement par complexe**, `-1` quand il est décroché,
    plutôt que la liste compacte des seuls complexes chargés. Ce n'est pas un
    détail de rangement : l'étage 3D tient une liaison OpenMM par emplacement et
    la déplace d'un monomère à la fois. Si l'identité des emplacements changeait
    d'un instantané à l'autre, une liaison sauterait d'une paire à une autre sans
    rapport, et la dynamique encaisserait un choc là où l'extrusion réelle
    n'avance que d'un cran.
    """

    legs: tuple[np.ndarray, ...]    # un tableau (n_lef, 2) par instantané, -1 = décroché
    n: int                          # monomères de la région
    n_lef: int
    lifetime: float                 # pas d'extrusion
    separation: int                 # pb par complexe
    processivity: int               # pb

    @property
    def snapshots(self) -> int:
        return len(self.legs)

    @property
    def occupancy(self) -> float:
        """Fraction des complexes effectivement chargés en moyenne."""
        return float(np.mean([(f[:, 0] >= 0).sum() for f in self.legs])) / self.n_lef

    @property
    def bonds(self) -> tuple[np.ndarray, ...]:
        """Les seules paires réellement tenues, instantané par instantané."""
        return tuple(f[f[:, 0] >= 0] for f in self.legs)

    def foot_density(self) -> np.ndarray:
        """Fréquence à laquelle chaque monomère porte un pied de complexe."""
        d = np.zeros(self.n)
        for b in self.bonds:
            if len(b):
                np.add.at(d, b.ravel(), 1.0)
        return d / max(self.snapshots, 1)

    def anchor_frequency(self, pairs: np.ndarray, *, tol: int = 2) -> np.ndarray:
        """Fréquence d'ancrage d'un complexe sur chacune des paires données.

        `tol` est une tolérance en monomères, et elle n'est pas cosmétique : la
        semaine 3 a mesuré qu'un appel de frontière ne tombe pas sur le bon casier
        — 38 % de rappel à tolérance nulle, 100 % à ±1. Comparer des pieds à des
        positions plantées sans tolérance mesurerait l'arrondi.
        """
        pairs = np.atleast_2d(np.asarray(pairs))
        hits = np.zeros(len(pairs))
        for b in self.bonds:
            if not len(b):
                continue
            lo, hi = b.min(axis=1), b.max(axis=1)
            for p, (i, j) in enumerate(pairs):
                if np.any((np.abs(lo - i) <= tol) & (np.abs(hi - j) <= tol)):
                    hits[p] += 1.0
        return hits / max(self.snapshots, 1)

    def loop_sizes(self, bp_per_bead: int) -> np.ndarray:
        """Taille des boucles tenues, en paires de bases, tous instantanés confondus."""
        sizes = [b[:, 1] - b[:, 0] for b in self.bonds if len(b)]
        if not sizes:
            return np.zeros(0)
        return np.concatenate(sizes).astype(np.float64) * bp_per_bead


def lef_count(region: Region, separation: int) -> int:
    return max(1, int(round(region.span / separation)))


def lifetime_steps(processivity: int, bp_per_bead: int, velocity: float = 1.0) -> float:
    """`processivité = 2 · vitesse · durée_de_vie`, convertie en pas d'extrusion."""
    return max(1.0, processivity / (2.0 * velocity * bp_per_bead))


def simulate(
    region: Region,
    sites: Sites,
    *,
    separation: int = 200_000,
    processivity: int = 200_000,
    velocity: float = 1.0,
    release: float = 0.003,
    burnin: int = 2_000,
    snapshots: int = 200,
    stride: int = 4,
    seed: int = 0,
) -> Extrusion:
    """Fait tourner le modèle 1D jusqu'à l'état stationnaire, puis échantillonne.

    Le rodage tourne **sans** dynamique 3D — il ne coûte rien et il évite de
    payer des dizaines de milliers de pas de Langevin pendant que la
    configuration 1D se remplit encore. `burnin` est en pas d'extrusion ; le
    défaut vaut plusieurs dizaines de durées de vie.

    `release` est la probabilité, par pas, qu'une jambe arrêtée reparte sans que
    le complexe se décroche. Elle borne le temps de séjour à un site par
    `1/release` pas, ce qui empêche une barrière quasi parfaite de figer
    définitivement un complexe et de vider le reste de la région.

    **Les bords de la région sont un artefact, et il faut le dire.** Une jambe
    arrivée au monomère 0 ou `n-1` ne peut plus avancer : elle s'y comporte comme
    devant une barrière parfaite, alors que le chromosome, lui, continue. Les
    premiers et derniers monomères voient donc une extrusion plus longue qu'ils
    ne devraient, et aucune conclusion ne se tire des domaines qui touchent un
    bord. C'est la raison pour laquelle la région choisie déborde largement le
    locus qu'on regarde.
    """
    rng = np.random.default_rng(seed)
    n = region.n
    stop_left, stop_right = sites.barriers(n)

    n_lef = lef_count(region, separation)
    lifetime = lifetime_steps(processivity, region.bp_per_bead, velocity)
    p_unload = 1.0 / lifetime

    left = np.full(n_lef, -1, dtype=np.int64)
    right = np.full(n_lef, -1, dtype=np.int64)
    held_l = np.zeros(n_lef, dtype=bool)     # jambe gauche arrêtée par un CTCF
    held_r = np.zeros(n_lef, dtype=bool)
    busy = np.zeros(n, dtype=bool)

    def unload(k: int) -> None:
        if left[k] >= 0:
            busy[left[k]] = False
            busy[right[k]] = False
        left[k] = right[k] = -1
        held_l[k] = held_r[k] = False

    def load(k: int) -> None:
        for _ in range(8):                   # quelques essais, puis on laisse tomber
            i = int(rng.integers(0, n - 1))
            if not busy[i] and not busy[i + 1]:
                left[k], right[k] = i, i + 1
                busy[i] = busy[i + 1] = True
                held_l[k] = rng.random() < stop_left[i]
                held_r[k] = rng.random() < stop_right[i + 1]
                return

    def advance() -> None:
        """Un pas : décrochage, chargement, relâchement, translocation — dans cet ordre.

        L'ordre compte. Relâcher avant de bouger donne à une jambe libérée sa
        chance au même pas ; l'inverse lui ferait perdre un pas et rallongerait
        artificiellement le temps de séjour aux barrières.
        """
        # 1. décrochage
        rolls = rng.random(n_lef)
        for k in range(n_lef):
            if left[k] >= 0 and rolls[k] < p_unload:
                unload(k)
        # 2. chargement des complexes libres
        for k in range(n_lef):
            if left[k] < 0:
                load(k)
        # 3. relâchement des jambes retenues
        free_l = rng.random(n_lef) < release
        free_r = rng.random(n_lef) < release
        np.logical_and(held_l, ~free_l, out=held_l)
        np.logical_and(held_r, ~free_r, out=held_r)
        # 4. translocation
        moves_l = rng.random(n_lef) < velocity
        moves_r = rng.random(n_lef) < velocity
        for k in range(n_lef):
            if left[k] < 0:
                continue
            if moves_l[k] and not held_l[k]:
                dest = left[k] - 1
                if dest >= 0 and not busy[dest]:
                    busy[left[k]] = False
                    busy[dest] = True
                    left[k] = dest
                    if rng.random() < stop_left[dest]:
                        held_l[k] = True
            if moves_r[k] and not held_r[k]:
                dest = right[k] + 1
                if dest < n and not busy[dest]:
                    busy[right[k]] = False
                    busy[dest] = True
                    right[k] = dest
                    if rng.random() < stop_right[dest]:
                        held_r[k] = True

    for _ in range(burnin):
        advance()

    frames: list[np.ndarray] = []
    for _ in range(snapshots):
        for _ in range(max(1, stride)):
            advance()
        frames.append(np.c_[left, right].astype(np.int64))

    return Extrusion(
        legs=tuple(frames),
        n=n,
        n_lef=n_lef,
        lifetime=lifetime,
        separation=separation,
        processivity=processivity,
    )
