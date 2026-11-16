"""La région fine : 2–4 Mb découpés en monomères de 1–5 kb, et ses sites CTCF.

Trois choses y sont fixées, et une seule d'entre elles est libre.

**Le rayon d'un monomère n'est pas un paramètre de la semaine 8.** Il vient de la
loi de la semaine 6 — volume de chromatine proportionnel aux paires de bases,
constante fixée par la fraction volumique `phi` mesurée par ChromEMT — et c'est
littéralement la même fonction qui le calcule (`nucleus.beads.capacity_at`).
C'est ce qui rend le raccord de fin de semaine falsifiable : si on avait le droit
de choisir la taille d'un monomère fin, on pourrait toujours faire coïncider les
deux modèles, et la coïncidence ne dirait rien.

**L'orientation des sites CTCF est le seul ingrédient du modèle qui produise des
boucles.** La règle de convergence (Rao 2014, Sanborn 2015) dit qu'une boucle
s'ancre entre un motif *forward* à gauche et un motif *reverse* à droite, les
deux pointant vers l'intérieur de la boucle. Traduite en barrières :

- un motif **forward (+)** arrête une jambe qui se déplace **vers la gauche** ;
- un motif **reverse (−)** arrête une jambe qui se déplace **vers la droite**.

Un cohésine chargé entre un `+` en `i` et un `−` en `j > i` voit donc ses deux
jambes se bloquer en `i` et `j`. Une paire divergente (`−` à gauche, `+` à
droite) n'arrête rien : les jambes lui passent au travers. C'est une prédiction
binaire, et le plan de validation est construit pour qu'elle puisse échouer.

**Ce qui est planté, et pourquoi.** Faute de réseau, aucune piste CTCF publiée
n'est lisible ici (`read_ctcf` l'attend, elle est écrite et testée). `plant`
fabrique donc une région de vérité connue — comme la semaine 3 pour les callers
Hi-C, et pour la même raison : sur des données réelles, un désaccord entre le
modèle et les boucles publiées ne dit pas lequel des deux a tort. La différence
avec une simple illustration tient au **témoin** : une fraction des domaines
reçoit des ancres divergentes ou en tandem, à la même distance génomique que les
convergents. Si le modèle produisait des points aux trois, il ne mesurerait rien.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..nucleus.beads import capacity_at
from ..nucleus.karyotype import gm12878

CONVERGENT, DIVERGENT, TANDEM = 0, 1, 2
KINDS = ("convergent", "divergent", "tandem")


def genome_bp() -> int:
    """Total diploïde de la lignée de référence, en paires de bases."""
    return int(sum(c.length for c in gm12878().copies))


def monomer_radius(
    bp_per_bead: int,
    *,
    nuclear_radius: float = 5.0,
    phi: float = 0.30,
    total_bp: int | None = None,
) -> float:
    """Rayon d'un monomère en nm, par la loi de la semaine 6 et pas une autre.

    À 750 kb elle rend les 166,8 nm du tableau de la semaine 6 ; à 2 kb, 23,2 nm.
    Le rayon suit `L^(1/3)`, donc le rapport entre deux résolutions est exact et
    un test le verrouille.
    """
    _, r_nm, _ = capacity_at(
        total_bp if total_bp is not None else genome_bp(),
        bp_per_bead,
        nuclear_radius=nuclear_radius,
        phi=phi,
    )
    return r_nm


def blob_radius(n: int, phi: float = 0.30) -> float:
    """Rayon, **en diamètres de monomère**, de la sphère que `n` monomères remplissent à `phi`.

    `R = r·(n/phi)^(1/3)`, soit `R/sigma = 0,5·(n/phi)^(1/3)`. C'est le volume que
    le modèle de noyau alloue à cette longueur de chromatine, ni plus ni moins :
    confiner le segment fin là-dedans, c'est lui imposer la densité du noyau sans
    lui imposer sa conformation.
    """
    return 0.5 * (n / phi) ** (1.0 / 3.0)


@dataclass(frozen=True)
class Region:
    """Un intervalle génomique et son découpage en monomères."""

    chrom: str
    start: int
    end: int
    bp_per_bead: int = 2_000

    @property
    def span(self) -> int:
        return self.end - self.start

    @property
    def n(self) -> int:
        return max(1, int(round(self.span / self.bp_per_bead)))

    def bead(self, pos: int) -> int:
        """Index du monomère qui contient une coordonnée génomique."""
        return int(np.clip((pos - self.start) // self.bp_per_bead, 0, self.n - 1))

    def bp(self, bead: int | np.ndarray) -> np.ndarray:
        """Milieu génomique d'un monomère."""
        return self.start + (np.asarray(bead) + 0.5) * self.bp_per_bead

    def radius_nm(self, phi: float = 0.30) -> float:
        return monomer_radius(self.bp_per_bead, phi=phi)

    def __str__(self) -> str:
        kb = self.bp_per_bead / 1000
        return (
            f"{self.chrom}:{self.start:,}-{self.end:,}  "
            f"{self.n:,} monomères de {kb:g} kb  "
            f"(rayon {self.radius_nm():.1f} nm)"
        )


@dataclass(frozen=True)
class Sites:
    """Des sites CTCF orientés, positionnés en index de monomères."""

    pos: np.ndarray          # (k,) index de monomère
    strand: np.ndarray       # (k,) +1 forward, -1 reverse
    occupancy: np.ndarray    # (k,) probabilité d'arrêt à la rencontre, dans [0, 1]
    source: str              # "planted:<seed>" ou le chemin du fichier lu

    @property
    def k(self) -> int:
        return len(self.pos)

    def barriers(self, n: int) -> tuple[np.ndarray, np.ndarray]:
        """Probabilités d'arrêt par monomère, séparées par sens de marche.

        Renvoie `(stop_leftward, stop_rightward)`. Les sites `+` remplissent le
        premier, les `−` le second, conformément à la règle de convergence. Deux
        sites de même sens sur le même monomère — possible à 2 kb, les motifs CTCF
        se groupent — ne s'additionnent pas : on garde le plus occupé, parce
        qu'une probabilité d'arrêt ne dépasse pas 1.
        """
        left = np.zeros(n)
        right = np.zeros(n)
        for p, s, o in zip(self.pos, self.strand, self.occupancy):
            target = left if s > 0 else right
            target[p] = max(target[p], float(o))
        return left, right

    def without_ctcf(self) -> "Sites":
        """Le même jeu de sites, tous inoccupés — le témoin négatif du modèle.

        Ce n'est pas « retirer les sites » : les positions restent, seule
        l'occupation tombe à zéro. Ce qui change entre les deux exécutions est
        donc exactement l'arrêt de l'extrusion, et rien d'autre.
        """
        return Sites(self.pos, self.strand, np.zeros(self.k), self.source + "+knockout")


@dataclass(frozen=True)
class Truth:
    """Ce qu'on a planté, donc ce qu'on a le droit de vérifier.

    Les brins ne sont pas stockés : la classe d'un domaine les détermine
    entièrement — convergent `+ −`, divergent `− +`, tandem `+ +`. Tout ce qui
    suit se recalcule donc depuis `kind` seul.
    """

    left: np.ndarray         # (d,) ancre gauche de chaque domaine, en monomères
    right: np.ndarray        # (d,) ancre droite
    kind: np.ndarray         # (d,) CONVERGENT | DIVERGENT | TANDEM
    boundary: np.ndarray     # (d-1,) monomère de la frontière entre domaines

    @property
    def strand_left(self) -> np.ndarray:
        return np.where(self.kind == DIVERGENT, -1, +1)

    @property
    def strand_right(self) -> np.ndarray:
        return np.where(self.kind == CONVERGENT, -1, +1)

    @property
    def blocking(self) -> np.ndarray:
        """Combien des deux sens de marche chaque frontière arrête : 1 ou 2.

        **Jamais 0**, et c'est tout l'intérêt de la propriété. Une frontière
        porte deux ancres — la droite du domaine `k` et la gauche du domaine
        `k+1` — et une ancre arrête toujours *l'un* des deux sens : un `+` les
        jambes qui vont vers la gauche, un `−` celles qui vont vers la droite.
        Une frontière « non bloquante » n'existe pas dans ce dispositif.

        C'est une erreur commise et sortie par les résultats eux-mêmes : la
        première version marquait une frontière bloquante si l'ancre de gauche
        était `−` **ou** celle de droite `+`, ce qui rate les deux autres cas et
        annonçait 18 frontières bloquantes sur 19. La dix-neuvième insulait
        évidemment autant, et le rapport imprimait ce « témoin » comme retrouvé
        — un témoin qui contredit son propre énoncé n'en est pas un.

        Ce qui distingue vraiment deux frontières est leur **force** : 2 quand
        les deux ancres sont de sens opposés, car les deux sens de marche sont
        alors arrêtés ; 1 sinon.
        """
        return np.where(self.strand_right[:-1] == self.strand_left[1:], 1, 2)

    @property
    def strong(self) -> np.ndarray:
        """Frontières qui arrêtent les deux sens — les seules qu'on attend nettes."""
        return self.blocking == 2

    def pairs(self, kind: int) -> np.ndarray:
        """Les paires d'ancres d'un type donné, `(m, 2)`."""
        who = self.kind == kind
        return np.c_[self.left[who], self.right[who]]

    def __str__(self) -> str:
        counts = [f"{KINDS[k]} {int((self.kind == k).sum())}" for k in (0, 1, 2)]
        return (
            f"{len(self.left)} domaines ({', '.join(counts)}) · "
            f"{int(self.strong.sum())}/{len(self.boundary)} frontières à deux sens bloqués"
        )


def plant(
    region: Region,
    *,
    seed: int = 0,
    mean_domain: int = 400_000,
    min_domain: int = 150_000,
    control: float = 0.35,
    occupancy: tuple[float, float] = (0.7, 0.95),
) -> tuple[Sites, Truth]:
    """Pave la région de domaines à ancres connues, témoins compris.

    Chaque domaine porte une ancre à chaque bout. Avec la probabilité `1 -
    control` elles sont **convergentes** (`+` à gauche, `−` à droite), ce qui est
    l'arrangement canonique d'une boucle CTCF. Sinon le domaine est un **témoin**,
    divergent (`−`, `+`) ou en tandem (`+`, `+`), tiré à pile ou face.

    Les trois classes partagent la même distribution de longueurs, et c'est tout
    l'intérêt : un score de point plus élevé aux convergents ne pourra pas être
    mis sur le compte de la distance génomique. Sans les témoins, on mesurerait
    « le modèle fait des points », ce qui est acquis dès qu'on pose des barrières.

    L'occupation est tirée dans `occupancy`, jamais 1 : un site CTCF réel n'est
    pas occupé en permanence, et une barrière parfaite rendrait les domaines
    étanches, ce que le Hi-C ne montre pas.
    """
    rng = np.random.default_rng(seed)
    n = region.n
    per_bead = region.bp_per_bead
    lo = max(2, min_domain // per_bead)
    mean = max(lo + 1, mean_domain // per_bead)

    edges = [0]
    while True:
        # Longueurs exponentielles décalées : les TADs ont une distribution à
        # queue lourde, pas une taille caractéristique nette.
        step = lo + int(rng.exponential(mean - lo))
        nxt = edges[-1] + step
        if nxt >= n - lo:
            break
        edges.append(nxt)
    edges.append(n)
    if len(edges) < 3:                       # région trop courte pour un témoin
        edges = [0, n // 2, n]

    left = np.array(edges[:-1], dtype=np.int64)
    right = np.array(edges[1:], dtype=np.int64) - 1

    kind = np.full(len(left), CONVERGENT, dtype=np.int64)
    is_control = rng.random(len(left)) < control
    # Les deux types de témoin **alternent** au lieu d'être tirés à pile ou face.
    # Sur une vingtaine de domaines, un tirage indépendant donne couramment 1 et 5
    # — mesuré — et un témoin à un seul membre ne contredit rien. L'équilibre des
    # effectifs est une propriété du plan d'expérience, pas du hasard.
    kind[is_control] = np.where(np.arange(int(is_control.sum())) % 2 == 0, DIVERGENT, TANDEM)

    strand_left = np.where(kind == DIVERGENT, -1, +1)
    strand_right = np.where(kind == CONVERGENT, -1, +1)

    pos = np.concatenate([left, right])
    strand = np.concatenate([strand_left, strand_right])
    occ = rng.uniform(occupancy[0], occupancy[1], len(pos))

    order = np.argsort(pos, kind="stable")
    sites = Sites(pos[order], strand[order], occ[order], f"planted:{seed}")

    truth = Truth(left=left, right=right, kind=kind, boundary=right[:-1])
    return sites, truth


def read_ctcf(path: str, region: Region, *, default_occupancy: float = 0.85) -> Sites:
    """Lit un BED de motifs CTCF orientés, restreint à la région.

    Colonnes attendues : `chrom start end [name] [score] strand`. Le brin est
    obligatoire — un site CTCF sans orientation ne contraint rien dans ce modèle,
    et en inventer une fabriquerait les boucles qu'on prétend prédire. Une ligne
    sans brin lisible est donc refusée, pas devinée.

    Si une colonne `score` est présente et tient dans [0, 1000], elle sert
    d'occupation (échelle BED standard) ; sinon `default_occupancy` s'applique à
    tous les sites, et l'appelant doit le savoir.
    """
    pos, strand, occ = [], [], []
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            f = line.split()
            if len(f) < 6:
                raise ValueError(f"{path}: 6 colonnes attendues (brin compris), vu {len(f)}")
            if f[0] != region.chrom:
                continue
            start, end = int(f[1]), int(f[2])
            mid = (start + end) // 2
            if not (region.start <= mid < region.end):
                continue
            if f[5] not in ("+", "-"):
                raise ValueError(f"{path}: brin illisible « {f[5]} », aucun défaut appliqué")
            pos.append(region.bead(mid))
            strand.append(+1 if f[5] == "+" else -1)
            try:
                score = float(f[4])
            except ValueError:
                score = -1.0
            occ.append(score / 1000.0 if 0.0 <= score <= 1000.0 else default_occupancy)

    if not pos:
        raise ValueError(f"{path}: aucun site CTCF dans {region.chrom}:{region.start}-{region.end}")

    order = np.argsort(pos, kind="stable")
    return Sites(
        np.asarray(pos)[order],
        np.asarray(strand)[order],
        np.clip(np.asarray(occ)[order], 0.0, 1.0),
        path,
    )
