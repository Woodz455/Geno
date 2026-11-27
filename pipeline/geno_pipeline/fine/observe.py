"""Ce qu'on mesure sur un ensemble de conformations fines, et à quoi on le compare.

Quatre observables, et une seule d'entre elles peut contredire le modèle de la
semaine 6.

**La carte de contacts** est ce qu'une expérience Hi-C verrait. On en tire les
frontières (score d'insulation) et les points d'angle (score de coin). Ces deux-là
se comparent à la vérité *plantée*, ce qui n'est pas une validation biologique
mais une validation de méthode — même logique qu'à la semaine 3, et la même
limite : sur des données réelles, un désaccord ne dit pas qui a tort.

**P(s)** se compare à ce que la semaine 7 a mesuré sur le noyau entier. Elle y
tombait à une pente de −0,10 au-delà de 15 Mb, c'est-à-dire plate, là où le Hi-C
réel continue de décroître. La question de la semaine 8 est de savoir si, en
dessous du TAD, le modèle fin retrouve une pente.

**R(s)** est le raccord. Les deux modèles se recouvrent entre 750 kb — la bille
de la semaine 6 — et l'étendue de la région fine. Sur cette fenêtre ils
prétendent tous les deux à la distance spatiale moyenne entre deux morceaux de
chromatine séparés de `s` paires de bases, et ils ne partagent **aucun
paramètre ajusté** : le rayon d'un monomère vient de la même loi, la fraction
volumique est la même, et rien dans la semaine 8 n'a été réglé sur la semaine 6.
S'ils se croisent, le raccord tient ; sinon on a une contradiction chiffrée, et
c'est un résultat.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from ..ensemble.stats import slope
from .region import KINDS, Region, Truth


# Les fenêtres d'ajustement de P(s), choisies **en regardant la courbe** et
# nommées pour ce qu'elles sont. La semaine 7 a payé l'erreur inverse : une
# fenêtre choisie en mégabases sans vérifier combien de billes cela faisait,
# tombée en plein régime de chaîne, et une pente qui ressemblait au Hi-C réel
# par accident.
REGIMES = (
    ("boucle", 20_000, 150_000),            # sous la taille des domaines
    ("domaine", 150_000, 800_000),          # l'échelle des barrières CTCF
    ("confinement", 1_000_000, 3_000_000),  # la sphère décide : artefact, pas résultat
)

# Le seul régime comparable au Hi-C publié est le deuxième. Lieberman-Aiden 2009
# rapporte `s^-1,08` sur **500 kb–7 Mb** — une fenêtre qu'une région de 4 Mb ne
# peut pas couvrir, ce qu'il faut dire en même temps que le chiffre. Au-delà d'environ 1,5 Mb sur une région de
# 4 Mb confinée dans sa propre sphère, ce qu'on mesure est la sphère — le nommer
# « confinement » plutôt que « territoire » évite de republier la P(s) plate de
# la semaine 7 comme si c'était une prédiction.
HIC_REFERENCE = -1.08


@dataclass(frozen=True)
class Map:
    """Une carte de contacts moyennée sur les conformations."""

    counts: np.ndarray        # (m, m) contacts par conformation
    frames: int
    bin_beads: int            # monomères par casier
    bp_per_bin: int
    cutoff: float

    @property
    def m(self) -> int:
        return len(self.counts)

    def bin_of(self, bead) -> np.ndarray:
        return np.asarray(bead) // self.bin_beads

    def expected(self) -> np.ndarray:
        """Moyenne par diagonale — le « attendu » d'une normalisation O/E."""
        m = self.m
        out = np.zeros(m)
        for d in range(m):
            out[d] = self.counts.diagonal(d).mean()
        return out

    def oe(self, floor: float = 1e-9) -> np.ndarray:
        """Carte observée/attendue. Retire la décroissance avec la distance.

        Sans elle, un « point » se lit d'autant plus fort qu'il est près de la
        diagonale, et comparer des paires d'ancres de longueurs différentes
        mesurerait leur longueur.
        """
        exp = np.maximum(self.expected(), floor)
        d = np.abs(np.subtract.outer(np.arange(self.m), np.arange(self.m)))
        return self.counts / exp[d]


def contact_map(
    coords: np.ndarray,
    *,
    cutoff: float = 1.5,
    bin_beads: int = 1,
    bp_per_bead: int = 2_000,
) -> Map:
    """Compte les paires à moins de `cutoff` diamètres de monomère.

    `cutoff` est en diamètres, c'est-à-dire en `r_i + r_j` : exactement la même
    convention qu'à la semaine 7 sur le noyau entier, où le seuil retenu valait
    1,5. Deux cartes mesurées avec deux définitions du contact ne se comparent
    pas, et la comparaison est tout l'objet de la semaine.
    """
    n = coords.shape[1]
    m = -(-n // bin_beads)
    counts = np.zeros((m, m))
    for frame in coords:
        pairs = cKDTree(np.asarray(frame, dtype=np.float64)).query_pairs(
            cutoff, output_type="ndarray"
        )
        if not len(pairs):
            continue
        a, b = pairs[:, 0] // bin_beads, pairs[:, 1] // bin_beads
        np.add.at(counts, (a, b), 1.0)
        np.add.at(counts, (b, a), 1.0)
    counts /= max(len(coords), 1)
    np.fill_diagonal(counts, counts.diagonal() / 2.0)   # comptée deux fois ci-dessus
    return Map(counts, len(coords), bin_beads, bin_beads * bp_per_bead, cutoff)


def insulation(cmap: Map, window: int) -> np.ndarray:
    """Score d'insulation : contacts traversant chaque casier, normalisés.

    `window` est en casiers. La semaine 3 a mesuré que ce choix décide du
    résultat — sur des TADs de 450 kb, le rappel des frontières passait de 100 %
    à 42 % en portant la fenêtre de 100 kb à 600 kb — donc la fenêtre doit rester
    nettement sous l'échelle des domaines, et l'appelant doit la donner.
    """
    m = cmap.m
    out = np.full(m, np.nan)
    c = cmap.counts
    for b in range(window, m - window):
        out[b] = c[b - window:b, b + 1:b + 1 + window].sum()
    good = np.isfinite(out) & (out > 0)
    if good.any():
        out = out / np.exp(np.log(out[good]).mean())     # moyenne géométrique
    return out


def call_boundaries(ins: np.ndarray, *, prominence: float = 0.1, spread: int = 3) -> np.ndarray:
    """Minima locaux du score d'insulation, au-delà d'une proéminence donnée.

    Un minimum local nu appelle du bruit ; on exige que le score remonte d'au
    moins `prominence` (en unités du score normalisé) de part et d'autre, dans
    une fenêtre de `spread` casiers.
    """
    m = len(ins)
    hits = []
    for b in range(spread, m - spread):
        v = ins[b]
        if not np.isfinite(v):
            continue
        left = ins[b - spread:b]
        right = ins[b + 1:b + 1 + spread]
        if not (np.isfinite(left).any() and np.isfinite(right).any()):
            continue
        if v <= np.nanmin(left) and v <= np.nanmin(right):
            rise = min(np.nanmax(left) - v, np.nanmax(right) - v)
            if rise >= prominence:
                hits.append(b)
    return np.array(hits, dtype=np.int64)


def dot_score(
    oe: np.ndarray,
    pairs: np.ndarray,
    *,
    inner: int = 1,
    outer: int = 6,
) -> np.ndarray:
    """Enrichissement en coin : le pixel contre son voisinage en anneau.

    Dans l'esprit de HiCCUPS — le « donut » — mais sur une carte déjà normalisée
    O/E, donc le fond ne porte plus la décroissance en distance. Reste à retirer
    le fond **local** : une paire d'ancres située à l'intérieur d'un domaine
    dense verrait sinon son score gonflé par le domaine.

    La croix passant par le pixel est exclue de l'anneau. La raison n'est pas
    qu'une bande produirait un faux point — elle n'en produit pas — mais
    l'inverse : une ancre de boucle émet une **bande** le long de sa ligne et de
    sa colonne, l'anneau tomberait en plein dedans, le fond serait surestimé et
    le point réel serait *masqué*. C'est la raison pour laquelle HiCCUPS retire
    la croix, et elle se vérifie par un point planté au croisement d'une bande.

    Ce n'est pas HiCCUPS : HiCCUPS compare le pixel à **quatre** fonds — le
    donut, le coin inférieur gauche, l'horizontale et la verticale — justement
    parce que le coin d'un domaine a un fond systématiquement différent de son
    intérieur. Un fond unique en anneau garde donc un biais. Il est **le même
    pour les trois classes d'ancres**, qui partagent la distribution de
    longueurs, donc il ne peut pas expliquer un écart entre elles : c'est la
    comparaison qui porte la conclusion, pas la valeur absolue du score.
    """
    m = len(oe)
    pairs = np.atleast_2d(np.asarray(pairs))
    out = np.full(len(pairs), np.nan)
    for p, (i, j) in enumerate(pairs):
        i, j = int(i), int(j)
        if not (outer <= i < m - outer and outer <= j < m - outer):
            continue
        peak = oe[i - inner:i + inner + 1, j - inner:j + inner + 1].mean()
        block = oe[i - outer:i + outer + 1, j - outer:j + outer + 1]
        mask = np.ones(block.shape, dtype=bool)
        c = outer
        mask[c - inner:c + inner + 1, :] = False          # la croix, lignes
        mask[:, c - inner:c + inner + 1] = False          # la croix, colonnes
        ring = block[mask]
        ring = ring[np.isfinite(ring) & (ring > 0)]
        if len(ring) < 8:
            continue
        out[p] = peak / np.median(ring)
    return out


@dataclass(frozen=True)
class Dots:
    """Scores de coin par classe d'ancres — le test qui peut échouer."""

    score: dict[str, np.ndarray]

    def summary(self) -> list[tuple[str, int, float, float]]:
        rows = []
        for name, v in self.score.items():
            v = v[np.isfinite(v)]
            if len(v):
                rows.append((name, len(v), float(np.median(v)), float(v.mean())))
            else:
                rows.append((name, 0, float("nan"), float("nan")))
        return rows

    def median(self, name: str) -> float:
        v = self.score.get(name, np.zeros(0))
        v = v[np.isfinite(v)]
        return float(np.median(v)) if len(v) else float("nan")

    def separation(self, against: str = "divergent") -> float:
        """Médiane des convergents divisée par celle d'une classe de référence.

        Un seul nombre, et il vaut 1 si l'orientation ne sert à rien.

        Le défaut est **divergent** et non « tous les témoins », parce que les
        témoins ne se valent pas. Un domaine divergent est un vrai négatif :
        aucune de ses deux ancres n'arrête la jambe qui l'atteint. Un domaine en
        tandem ne l'est qu'à moitié — son ancre gauche `+` arrête bien la jambe
        gauche, et si le domaine suivant commence par un `−`, ce `−` se trouve à
        un monomère de l'ancre droite et arrête la jambe droite. Le « témoin »
        porte alors une vraie boucle, ancrée un cran plus loin, que la tolérance
        du score ramasse. Mesuré : un domaine en tandem sur trois marque 3,27,
        autant que le meilleur convergent.

        Ce n'est pas un défaut du modèle, c'est un défaut du témoin, et le
        nommer vaut mieux que moyenner les deux classes ensemble.
        """
        conv, ref = self.median(KINDS[0]), self.median(against)
        if not np.isfinite(conv) or not np.isfinite(ref) or ref == 0:
            return float("nan")
        return float(conv / ref)


def dots_by_kind(cmap: Map, truth: Truth, *, null: int = 400, seed: int = 0, **kw) -> Dots:
    """Score de coin pour chaque classe d'ancres plantée, plus un fond de paires au hasard.

    Le fond n'est pas décoratif. Les classes témoins sont peu nombreuses — sur une
    vingtaine de domaines, trois divergents et trois tandem — et une médiane sur
    deux valeurs ne dit pas grand-chose. Les paires tirées au hasard, elles, se
    comptent par centaines et **à distances appariées** : pour chaque ancre
    plantée on tire `null / d` paires de même séparation ailleurs dans la région.
    C'est ce qui permet de dire si 2,5 est un grand nombre.
    """
    oe = cmap.oe()
    rng = np.random.default_rng(seed)
    out: dict[str, np.ndarray] = {}
    for k, name in enumerate(KINDS):
        pairs = truth.pairs(k)
        out[name] = (
            dot_score(oe, np.c_[cmap.bin_of(pairs[:, 0]), cmap.bin_of(pairs[:, 1])], **kw)
            if len(pairs)
            else np.zeros(0)
        )

    # Les tirages restent dans la zone **scorable** : `dot_score` rend NaN dès
    # qu'un pixel est à moins de `outer` casiers d'un bord, et un fond de NaN ne
    # sert à rien. Sur une petite carte avec de longues portées, l'intervalle
    # peut être vide — on saute alors cette portée plutôt que de tirer à blanc.
    pad = int(kw.get("outer", 6))
    spans = cmap.bin_of(truth.right) - cmap.bin_of(truth.left)
    draws = []
    per = max(1, null // max(len(spans), 1))
    for span in spans:
        span = int(span)
        hi = cmap.m - pad - span
        if span <= 0 or hi <= pad:
            continue
        i = rng.integers(pad, hi, per)
        draws.append(np.c_[i, i + span])
    out["hasard"] = (
        dot_score(oe, np.concatenate(draws), **kw) if draws else np.zeros(0)
    )
    return Dots(out)


def ps(cmap: Map, *, min_sep: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """P(s) depuis la carte : moyenne par diagonale, normalisée à s minimal.

    Les deux dernières décades de séparation reposent sur de moins en moins de
    paires — à `s = m - 1` il n'en reste qu'une — donc la queue est bruitée par
    construction et une pente ne s'y ajuste pas.
    """
    exp = cmap.expected()
    s = np.arange(cmap.m) * cmap.bp_per_bin
    keep = np.arange(cmap.m) >= min_sep
    p = exp / max(exp[min_sep], 1e-12)
    return s[keep], p[keep]


def separation_curve(
    coords: np.ndarray,
    *,
    bp_per_bead: int,
    unit_nm: float,
    copy_id: np.ndarray | None = None,
    seps: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """`R(s)` : distance spatiale moyenne entre monomères séparés de `s`, en nm.

    Marche pour les deux modèles, et c'est le but. Sur la chaîne fine, `copy_id`
    est nul et toutes les paires comptent. Sur le noyau entier, `copy_id` désigne
    les 46 copies et seules les paires **intra-copie** sont retenues : une
    distance entre chromosomes différents ne correspond à aucune séparation
    génomique.

    Renvoie `(séparations en pb, moyenne en nm, écart-type en nm)`.
    """
    coords = np.asarray(coords, dtype=np.float64)
    n = coords.shape[1]
    if seps is None:
        seps = np.unique(np.geomspace(1, max(n - 1, 2), 28).round().astype(int))
    seps = np.asarray([s for s in seps if 1 <= s < n], dtype=np.int64)

    mean = np.empty(len(seps))
    sd = np.empty(len(seps))
    for t, s in enumerate(seps):
        d = np.linalg.norm(coords[:, s:] - coords[:, :-s], axis=2)
        if copy_id is not None:
            d = d[:, copy_id[s:] == copy_id[:-s]]
        mean[t] = d.mean() * unit_nm
        sd[t] = d.std() * unit_nm
    return seps * bp_per_bead, mean, sd


@dataclass(frozen=True)
class Junction:
    """Le raccord entre l'échelle fine et le noyau entier, en chiffres."""

    bp: np.ndarray            # séparations comparées
    fine_nm: np.ndarray
    coarse_nm: np.ndarray
    fine_pairs: np.ndarray    # paires de monomères par conformation, à cette séparation

    @property
    def thin(self) -> np.ndarray:
        """Les points qui reposent sur peu de paires, donc à lire avec prudence.

        À la séparation `s` sur `n` monomères il ne reste que `n − s` paires : le
        dernier point de recouvrement d'une région de 4 Mb découpée en 2 kb en
        compte 125 par conformation contre 1 625 pour le premier. Ce n'est pas
        une raison de le cacher, c'en est une de dire combien.
        """
        return self.fine_pairs < 0.2 * self.fine_pairs.max()

    @property
    def ratio(self) -> np.ndarray:
        return self.fine_nm / np.maximum(self.coarse_nm, 1e-9)

    @property
    def worst(self) -> float:
        r = self.ratio
        return float(max(r.max(), 1.0 / r.min()))

    def exponents(self) -> tuple[float, float]:
        """Pentes log-log de `R(s)` pour les deux modèles, sur la fenêtre commune.

        C'est le fond du raccord, et c'est plus parlant qu'un rapport : deux
        courbes peuvent se croiser en un point et diverger partout ailleurs. Les
        repères utiles sont 0,5 pour une marche aléatoire idéale, ~0,588 pour une
        marche auto-évitante gonflée, et 1/4 à 1/3 pour ce que le traçage de
        chromatine mesure au-dessus du mégabase (Wang 2016, Bintu 2018, Su 2020) —
        une chromatine réelle est nettement plus compacte qu'une marche.
        """
        if len(self.bp) < 2:
            return float("nan"), float("nan")
        x = np.log(self.bp.astype(float))
        fit = lambda y: float(np.polyfit(x, np.log(np.maximum(y, 1e-9)), 1)[0])
        return fit(self.fine_nm), fit(self.coarse_nm)

    def holds(self, tol: float = 1.25) -> bool:
        """Les deux modèles se raccordent-ils à `tol` près sur toute la fenêtre ?

        La tolérance est un facteur, pas une différence : à ces échelles, les
        distances publiées en FISH s'accordent rarement mieux qu'à 20–30 % entre
        laboratoires, et exiger de deux modèles indépendants qu'ils tombent à
        mieux que ça serait exiger plus que la mesure.
        """
        return bool(np.isfinite(self.ratio).all() and self.worst <= tol)


def junction(
    fine: np.ndarray,
    region: Region,
    coarse: np.ndarray,
    *,
    coarse_copy_id: np.ndarray,
    coarse_bp_per_bead: int,
    sigma_nm: float,
    coarse_unit_nm: float = 1_000.0,
    max_bp: int | None = None,
) -> Junction:
    """Compare `R(s)` des deux modèles sur leur fenêtre commune.

    La fenêtre commune commence à une bille de la semaine 6 — en dessous, le
    modèle grossier n'a rien à dire — et s'arrête **strictement avant** l'étendue
    de la région fine : sur `n` monomères, la séparation `n` n'existe pas, il n'y
    a aucune paire à cette distance. Une région de 3 Mb découpée en 750 kb donne
    donc trois points de comparaison et non quatre, et il vaut mieux le calculer
    que laisser une troncature silencieuse le faire.

    Les séparations du modèle grossier tombent sur des multiples de 750 kb ; on
    demande au modèle fin exactement les mêmes, ce qui évite d'interpoler l'une
    des deux courbes pour la faire passer par l'autre.
    """
    span = region.n * region.bp_per_bead
    limit = min(span, max_bp if max_bp else span)
    k = 0
    while (k + 1) * coarse_bp_per_bead <= limit and round(
        (k + 1) * coarse_bp_per_bead / region.bp_per_bead
    ) < region.n:
        k += 1
    if k == 0:
        raise ValueError(
            f"aucune séparation commune : région de {span:,} pb contre des billes "
            f"de {coarse_bp_per_bead:,} pb"
        )

    seps = np.arange(1, k + 1)
    coarse_bp, coarse_nm, _ = separation_curve(
        coarse,
        bp_per_bead=coarse_bp_per_bead,
        unit_nm=coarse_unit_nm,
        copy_id=coarse_copy_id,
        seps=seps,
    )
    fine_seps = np.round(coarse_bp / region.bp_per_bead).astype(int)
    _, fine_nm, _ = separation_curve(
        fine,
        bp_per_bead=region.bp_per_bead,
        unit_nm=sigma_nm,
        seps=fine_seps,
    )
    return Junction(coarse_bp, fine_nm, coarse_nm, np.asarray(fine.shape[1] - fine_seps))


@dataclass(frozen=True)
class Stationarity:
    """La trajectoire a-t-elle fini de se souvenir de sa conformation initiale ?"""

    bp: np.ndarray
    first_nm: np.ndarray
    second_nm: np.ndarray

    @property
    def drift(self) -> np.ndarray:
        return self.second_nm / np.maximum(self.first_nm, 1e-9) - 1.0

    @property
    def worst(self) -> float:
        return float(np.abs(self.drift).max())

    def holds(self, tol: float = 0.05) -> bool:
        return bool(np.isfinite(self.drift).all() and self.worst <= tol)


def stationarity(
    coords: np.ndarray,
    *,
    per_replicate: int,
    bp_per_bead: int,
    unit_nm: float,
    seps: np.ndarray | None = None,
) -> Stationarity:
    """Compare `R(s)` sur la première et la seconde moitié des instantanés.

    C'est le contrôle d'équilibrage, fait **sur les données de production** plutôt
    que sur une étude séparée. Il vaut mieux, pour deux raisons : il porte sur la
    trajectoire qu'on publie, et il ne coûte rien.

    Le mode lent d'un polymère confiné est `R(s)` aux grandes séparations —
    `Rg` n'y voit rien, puisqu'une conformation initiale comprimée a déjà le bon
    `Rg` par construction (§ `polymer.start`). Si les deux moitiés diffèrent, la
    conformation initiale décide encore et le raccord de fin de semaine mesure
    le point de départ plutôt que le modèle.

    `per_replicate` sert à couper **dans chaque réplicat** : les instantanés sont
    empilés bout à bout, donc couper l'ensemble en deux comparerait des réplicats
    entre eux, ce qui est une autre question.
    """
    coords = np.asarray(coords)
    half = per_replicate // 2
    reps = len(coords) // per_replicate
    first = np.concatenate(
        [coords[r * per_replicate: r * per_replicate + half] for r in range(reps)]
    )
    second = np.concatenate(
        [coords[r * per_replicate + half: (r + 1) * per_replicate] for r in range(reps)]
    )
    bp, a, _ = separation_curve(first, bp_per_bead=bp_per_bead, unit_nm=unit_nm, seps=seps)
    _, b, _ = separation_curve(second, bp_per_bead=bp_per_bead, unit_nm=unit_nm, seps=seps)
    return Stationarity(bp, a, b)


def regime(s: np.ndarray, p: np.ndarray, lo: float, hi: float) -> float:
    """Pente log-log sur une fenêtre choisie **en regardant la courbe**.

    C'est la leçon de la semaine 7 : la première P(s) du noyau avait été ajustée
    sur une fenêtre choisie en mégabases sans vérifier combien de billes cela
    faisait, la fenêtre couvrait le régime de chaîne, et la pente obtenue —
    −0,82 — ressemblait au Hi-C réel par accident.
    """
    return slope(s, p, lo, hi)
