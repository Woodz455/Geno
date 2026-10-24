"""Polymère 3D connu, et la matrice de contacts qui en découle.

La semaine 3 plantait des *enrichissements de contacts* : compartiments, TADs,
boucles. Elle ne plantait aucune position dans l'espace. Or la semaine 5 doit
répondre à une question géométrique — quelle conversion contact → distance
restitue la vraie structure — et pour ça il faut une vraie structure.

D'où l'inversion : on fabrique d'abord une conformation 3D, puis on en **dérive**
la matrice de contacts par un modèle direct explicite. La matrice devient une
observation d'une vérité qu'on connaît, et la reconstruction devient falsifiable.

Modèle direct, un seul paramètre :

    f(i, j) ∝ d(i, j)^(-gamma)

`gamma` est planté. La reconstruction, elle, balaie `alpha` dans
`d ∝ f^(-alpha)` sans le connaître. Si la méthode est correcte, l'optimum doit
tomber sur **alpha = 1 / gamma**. C'est une prédiction chiffrée qui peut échouer,
pas une impression.

Organisation spatiale : les compartiments A occupent le centre du noyau, les B la
périphérie. Ce n'est pas décoratif — c'est l'observation de Cremer & Cremer sur
les territoires, et c'est ce qui rend la structure plantée vérifiable contre un
fait indépendant.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Conformation:
    """Une conformation 3D et ce qui la caractérise."""

    coords: np.ndarray        # (n, 3), unités arbitraires
    compartment: np.ndarray   # +1 = A, -1 = B, par bille
    radius: float             # rayon de confinement

    @property
    def n(self) -> int:
        return len(self.coords)

    def radial(self) -> np.ndarray:
        """Position radiale normalisée : 0 au centre, 1 à la périphérie."""
        centre = self.coords.mean(axis=0)
        return np.linalg.norm(self.coords - centre, axis=1) / self.radius

    def distances(self) -> np.ndarray:
        d = np.linalg.norm(self.coords[:, None, :] - self.coords[None, :, :], axis=2)
        np.fill_diagonal(d, 0.0)
        return d


def _blocks(n: int, lo: int, hi: int, rng: np.random.Generator) -> list[tuple[int, int]]:
    edges = [0]
    while edges[-1] + 2 * lo <= n:
        nxt = edges[-1] + int(rng.integers(lo, hi + 1))
        if nxt + lo > n:
            break
        edges.append(nxt)
    return [(edges[k], edges[k + 1] if k + 1 < len(edges) else n) for k in range(len(edges))]


def chain(
    n: int = 400,
    seed: int = 0,
    *,
    radius: float = 10.0,
    step: float = 1.0,
    persistence: float = 0.65,
    drift: float = 0.12,
    relax_iters: int = 120,
    excluded: float = 0.55,
    radial_pull: float = 0.06,
) -> Conformation:
    """Marche persistante confinée, avec ségrégation A/B et volume exclu.

    La persistance donne une fibre plutôt qu'un nuage ; la dérive attire chaque
    bloc vers le rayon typique de son compartiment ; la relaxation empêche les
    billes de se superposer, ce qui rendrait les distances courtes irréalistes.
    """
    rng = np.random.default_rng(seed)

    # Les blocs suivent la longueur de la chaîne : sur une chaîne courte, des
    # blocs de taille fixe n'en laissent que deux ou trois, et la ségrégation A/B
    # n'a pas la place de s'établir.
    compartment = np.empty(n, dtype=np.int8)
    lo_len, hi_len = max(6, n // 14), max(14, n // 6)
    for k, (lo, hi) in enumerate(_blocks(n, lo_len, hi_len, rng)):
        compartment[lo:hi] = 1 if k % 2 == 0 else -1

    # A au centre, B en périphérie — territoires de Cremer & Cremer.
    target_r = np.where(compartment == 1, 0.35 * radius, 0.80 * radius)

    coords = np.zeros((n, 3))
    coords[0] = rng.normal(0, 1, 3)
    direction = rng.normal(0, 1, 3)
    direction /= np.linalg.norm(direction)

    for i in range(1, n):
        direction = persistence * direction + (1 - persistence) * rng.normal(0, 1, 3)
        direction /= np.linalg.norm(direction)

        r = np.linalg.norm(coords[i - 1])
        if r > 1e-6:
            radial_dir = coords[i - 1] / r
            pull = (target_r[i] - r) / radius
            direction = direction + drift * pull * radial_dir * 3.0
            direction /= np.linalg.norm(direction)

        coords[i] = coords[i - 1] + step * direction
        over = np.linalg.norm(coords[i])
        if over > radius:                       # confinement : on réfléchit vers l'intérieur
            coords[i] *= radius / over * 0.98

    # Recentrage avant relaxation. La dérive vise un rayon depuis l'origine,
    # mais toute mesure faite sur une reconstruction part forcément du centroïde
    # — le centre du noyau n'y est pas connu. Si les deux ne coïncident pas, la
    # ségrégation A/B plantée devient invisible à la mesure qui la cherche, et le
    # test échoue sans que rien ne soit faux.
    coords -= coords.mean(axis=0)

    # Relaxation : volume exclu, chaîne tenue, et ségrégation A/B.
    #
    # La dérive appliquée pendant la marche ne suffit pas : elle lutte contre la
    # diffusion du hasard, et le résultat dépend de la longueur de chaîne et de la
    # graine. Mesuré avant correction : écart radial A/B de +0,157 à 180 billes
    # mais +0,009 à 350. Une propriété qui ne tient qu'à certaines tailles n'est
    # pas une propriété. Ici c'est une force de rappel, donc une propriété
    # convergée.
    for _ in range(relax_iters):
        r = np.linalg.norm(coords, axis=1)
        safe = np.maximum(r, 1e-9)
        # Ferme `radial_pull` de l'écart au rayon cible à chaque itération.
        coords += radial_pull * (target_r - r)[:, None] * (coords / safe[:, None])

        delta = coords[:, None, :] - coords[None, :, :]
        dist = np.linalg.norm(delta, axis=2)
        np.fill_diagonal(dist, np.inf)
        too_close = dist < excluded
        if not too_close.any():
            break
        push = np.where(too_close[:, :, None], delta / np.maximum(dist, 1e-9)[:, :, None], 0.0)
        coords += 0.12 * push.sum(axis=1)

        # Les liaisons de chaîne se rappellent au bon souvenir des billes.
        bond = coords[1:] - coords[:-1]
        length = np.linalg.norm(bond, axis=1, keepdims=True)
        correction = 0.35 * bond * (1 - step / np.maximum(length, 1e-9))
        coords[:-1] += correction
        coords[1:] -= correction

        coords -= coords.mean(axis=0)
        out = np.linalg.norm(coords, axis=1)
        over = out > radius
        if over.any():
            coords[over] *= (radius / out[over])[:, None] * 0.99

    return Conformation(coords=coords, compartment=compartment, radius=radius)


def contacts(
    conf: Conformation,
    *,
    gamma: float = 3.0,
    total: int = 4_000_000,
    seed: int = 0,
    bias_sd: float = 0.0,
) -> np.ndarray:
    """Matrice de comptages échantillonnée depuis la conformation.

    `gamma` est l'exposant du modèle direct, et **c'est lui que la reconstruction
    devra retrouver** sous la forme alpha = 1/gamma. `bias_sd` ajoute un biais de
    couverture par bille, à retirer par équilibrage avant toute reconstruction.
    """
    rng = np.random.default_rng(seed)
    d = conf.distances()
    n = conf.n

    with np.errstate(divide="ignore"):
        expected = np.where(d > 0, d ** (-gamma), 0.0)
    np.fill_diagonal(expected, 0.0)

    if bias_sd > 0:
        bias = np.exp(rng.normal(0, bias_sd, n))
        expected *= bias[:, None] * bias[None, :]

    expected *= total / expected.sum()
    counts = rng.poisson(expected)
    counts = np.triu(counts) + np.triu(counts, 1).T   # symétrique
    np.fill_diagonal(counts, 0)
    return counts.astype(np.float64)
