"""Le locus WNT6/IHH/EPHA4/PAX3, en hg38 : ce qui est mesuré, et ce qui est reconstruit.

**Mesuré ou annoté** (sources verrouillées) : les gènes (GENCODE v50), les sites CTCF de
GM12878 (pics IDR d'ENCODE, orientés par le scan MA0139.1 de la semaine 13), les domaines et
boucles de Rao 2014. Dans GM12878, le TAD d'EPHA4 est la boucle 219,74 ↔ 221,57 Mb, ancrée
sur deux grappes CTCF convergentes : `+` à 219,736 Mb, `−` au promoteur d'EPHA4 (221,572 Mb).
C'est cette grappe-là, au bord télomérique du TAD, qui sépare EPHA4 de PAX3.

**Reconstruit depuis le texte de Lupiáñez 2015** — la Table S1, qui donne les points de cassure
des patients, n'a pas pu être lue. Le texte contraint :

- les délétions de brachydactylie font **1,75 à 1,9 Mb**, emportent EPHA4 et la frontière
  EPHA4/PAX3, et s'arrêtent dans la partie **non codante** du TAD de PAX3 ;
- les enhancers de membre forment une grappe de **30 kb** dans une région de 150 kb qui contacte
  le promoteur d'EPHA4 « across 1.66 Mb of the intervening gene desert » (distance mesurée chez
  la souris, appliquée ici à l'humain) ; PAX3 est mal exprimé dans le motif de ces enhancers,
  qui survivent donc à la délétion ;
- le témoin DelBS a la même forme, mais laisse la région de la frontière et sa grappe CTCF.

On ne tire donc pas **un** variant : toute la famille de points de cassure que le texte
autorise, et la conclusion doit tenir sur chacun.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..fine.region import Region, Sites

CHROM = "chr2"
REGION = Region(CHROM, 218_000_000, 223_000_000, 2_000)
# Fenêtre constante (révision, voir `sites`) : une délétion fait entrer la séquence en aval dans
# la fenêtre simulée ; 2 Mb de marge couvrent la plus longue (1,9 Mb).
REACH = Region(CHROM, REGION.start, REGION.end + 2_000_000, REGION.bp_per_bead)
ENHANCER_OFFSET = 1_660_000          # bp, du promoteur d'EPHA4 vers le centromère (texte)
ENHANCER_SPAN = 30_000               # la grappe des trois enhancers actifs (texte)
DEL_LENGTHS = (1_750_000, 1_800_000, 1_850_000, 1_900_000)   # « 1.75–1.9 Mb » (texte)
CEN_STEP = 50_000
BOUNDARY_WINDOW = 20_000             # sites CTCF orientés à ± 20 kb du TSS d'EPHA4
PROMOTER_HALF = 4_000                # ± 4 kb autour d'un TSS : 4 à 5 monomères de 2 kb
GENES = ("WNT6", "WNT10A", "IHH", "NHEJ1", "EPHA4", "PAX3", "SGPP2")
SOURCE = ("Lupiáñez DG et al., Cell 161:1012-1025 (2015), doi:10.1016/j.cell.2015.04.004 — "
          "texte seulement ; Table S1 (points de cassure des patients) non lue")


@dataclass(frozen=True)
class Gene:
    name: str
    start: int
    end: int
    strand: str

    @property
    def tss(self) -> int:
        return self.start if self.strand == "+" else self.end - 1


@dataclass(frozen=True)
class Peak:
    start: int
    end: int
    motif: str | None
    signal: float

    @property
    def mid(self) -> int:
        return (self.start + self.end) // 2


@dataclass(frozen=True)
class Variant:
    """Une délétion `[start, end)` en hg38, et ce qu'elle épargne. `keep` : des positions
    génomiques dont les sites CTCF sont **gardés** même dans l'intervalle (témoin in silico)."""

    name: str
    start: int
    end: int
    kind: str                 # "DEL-B", "DelBS", "DEL-B, frontière épargnée"
    keep: tuple[tuple[int, int], ...] = ()
    note: str = ""

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class Locus:
    region: Region
    genes: dict[str, Gene]
    peaks: list[Peak]
    p90: float                               # signal CTCF du 90e centile, génome entier
    enhancer: tuple[int, int]
    boundary: tuple[int, int]                # la grappe CTCF du bord télomérique du TAD d'EPHA4
    tad: tuple[int, int]                     # la boucle de Rao 2014 qui le borne
    variants: list[Variant] = field(default_factory=list)
    reach: Region | None = None              # la séquence au-delà de la région, pour `window`

    # — sites —
    def occupancy(self, signal: float) -> float:
        """Signal ChIP → probabilité d'arrêt : 0,95 au 90e centile du génome et au-delà,
        proportionnelle en dessous, plancher 0,05. La plage des sites plantés (0,7–0,95) est
        celle des pics forts ; un pic faible doit pouvoir laisser passer une cohésine."""
        return float(np.clip(0.95 * signal / self.p90, 0.05, 0.95))

    def sites(self, variant: Variant | None = None, *, constant: float | None = None,
              window: bool = False) -> tuple[Sites, np.ndarray]:
        """Les sites orientés de l'allèle — sauvage, ou délété et recollé —, et la position
        génomique (hg38, milieu) de chaque monomère de cet allèle.

        `window=False` (le plan pré-enregistré) : l'allèle est la région moins la délétion, une
        chaîne plus courte. `window=True` (la révision) : l'allèle garde les `region.n` premiers
        monomères recollés, la délétion faisant entrer la séquence en aval (`reach`) — même
        longueur, donc même confinement, que le type sauvage. Le confinement est réglé sur la
        quantité de chromatine : une chaîne plus courte dans une sphère plus petite a un
        plancher de contacts lointains plus haut (~0,0020 contre ~0,0012 au-delà de 600 kb),
        qui se confond avec l'effet de la frontière."""
        base = self.reach if window else self.region
        if base is None:
            raise ValueError("window=True demande la séquence au-delà de la région (reach)")
        bp = base.bp_per_bead
        mids = base.bp(np.arange(base.n)).astype(np.int64)
        keep = np.ones(base.n, bool)
        if variant is not None:
            keep &= ~((mids >= variant.start) & (mids < variant.end))
        if window:
            keep[np.flatnonzero(keep)[self.region.n:]] = False
        index = np.full(base.n, -1, np.int64)
        index[keep] = np.arange(int(keep.sum()))
        pos, strand, occ = [], [], []
        for p in self.peaks:
            if p.motif is None:
                continue
            k = int((p.mid - base.start) // bp)
            if not (0 <= k < base.n):
                continue
            j = index[k]
            if j < 0:
                kept = (variant is not None and any(a <= p.mid < b for a, b in variant.keep)
                        and variant.start <= p.mid < variant.end)
                if not kept:
                    continue
                # Site épargné dans la délétion : posé au point de recollement.
                j = int(index[: k].max()) if (index[: k] >= 0).any() else 0
            pos.append(int(j))
            strand.append(+1 if p.motif == "+" else -1)
            occ.append(constant if constant is not None else self.occupancy(p.signal))
        order = np.argsort(pos, kind="stable")
        src = "ENCODE GM12878 CTCF (ctcf.peaks.gm12878), orienté MA0139.1"
        return (Sites(np.asarray(pos, np.int64)[order], np.asarray(strand, np.int64)[order],
                      np.asarray(occ, float)[order], src), mids[keep])

    def allele_region(self, variant: Variant | None, *, window: bool = False) -> Region:
        """La région à simuler : contiguë, de la longueur de l'allèle recollé (celle du type
        sauvage avec `window`)."""
        n = self.region.n - (0 if variant is None or window else
                             int(((self.region.bp(np.arange(self.region.n)) >= variant.start) &
                                  (self.region.bp(np.arange(self.region.n)) < variant.end)).sum()))
        return Region(self.region.chrom, self.region.start, self.region.start + n * self.region.bp_per_bead,
                      self.region.bp_per_bead)

    # — repères en monomères d'un allèle —
    @staticmethod
    def monomers(mids: np.ndarray, a: int, b: int) -> np.ndarray:
        return np.flatnonzero((mids >= a) & (mids < b))

    def promoter(self, gene: str, mids: np.ndarray) -> np.ndarray:
        t = self.genes[gene].tss
        return self.monomers(mids, t - PROMOTER_HALF, t + PROMOTER_HALF)


def build(m=None) -> Locus:
    """Le locus depuis les sources verrouillées, et la famille de variants du texte."""
    from ..fiche import sources as S
    from ..sequence import TwoBit

    m = m or S.manifest()
    genes_iv, _ = S.genes(m)
    genes: dict[str, Gene] = {}
    for s, e, i in zip(genes_iv.starts[CHROM], genes_iv.ends[CHROM], genes_iv.idx[CHROM]):
        g = genes_iv.items[int(i)]
        if g["name"] in GENES and REGION.start <= s < REGION.end:
            genes[g["name"]] = Gene(g["name"], int(s), int(e), g["strand"])
    missing = set(GENES) - set(genes)
    if missing:
        raise ValueError(f"gènes absents de GENCODE dans la région : {sorted(missing)}")

    ct = S.ctcf(m, TwoBit(S.locked(m, "hg38.2bit").path))
    p90 = float(np.percentile([d["signal"] for d in ct.items], 90))
    peaks = [Peak(int(s), int(e), ct.items[int(i)]["motif"], float(ct.items[int(i)]["signal"]))
             for s, e, i in zip(ct.starts[CHROM], ct.ends[CHROM], ct.idx[CHROM])
             if REACH.start <= s < REACH.end]

    epha4 = genes["EPHA4"]
    near = [p for p in peaks if p.motif is not None and abs(p.mid - epha4.tss) <= BOUNDARY_WINDOW]
    if not near:
        raise ValueError("aucun site CTCF orienté au promoteur d'EPHA4 : pas de frontière à supprimer")
    boundary = (min(p.start for p in near), max(p.end for p in near))

    loops, _ = S.loops(m, S.chain_hg19(m))
    tad = None
    for s, e, i in zip(loops.starts[CHROM], loops.ends[CHROM], loops.idx[CHROM]):
        d = loops.items[int(i)]
        # La boucle la plus longue ancrée sur la grappe du promoteur : le TAD entier, pas une
        # boucle intérieure qui partage la même ancre télomérique.
        if (d["partner"][0] == CHROM and s <= boundary[1] and e >= boundary[0] - 10_000 and d["partner"][1] < s
                and (tad is None or d["partner"][1] < tad[0])):
            tad = (int(d["partner"][1]), int(e))
    if tad is None:
        raise ValueError("la boucle de Rao 2014 qui borne le TAD d'EPHA4 est introuvable")

    c = epha4.tss - ENHANCER_OFFSET
    enh = (c - ENHANCER_SPAN // 2, c + ENHANCER_SPAN // 2)
    loc = Locus(REGION, genes, peaks, p90, enh, boundary, tad, reach=REACH)
    loc.variants = family(loc)
    return loc


def family(loc: Locus) -> list[Variant]:
    """Toutes les délétions DEL-B que le texte autorise, et pour chacune ses deux témoins.

    DEL-B : longueur dans DEL_LENGTHS, extrémité centromérique au-delà de la grappe
    d'enhancers (qui survit), extrémité télomérique au-delà de la frontière et en deçà de la
    fin 3' de PAX3 (partie non codante de son TAD). Témoins, même extrémité centromérique :
    DelBS s'arrête juste avant la grappe CTCF de la frontière ; « frontière épargnée » a les
    mêmes points de cassure que DEL-B mais garde les sites CTCF de la grappe — même distance
    génomique, seule la frontière diffère."""
    pax3 = loc.genes["PAX3"]
    lo = loc.enhancer[1]
    hi = pax3.start                               # PAX3 sur le brin − : sa fin 3' est son début
    out: list[Variant] = []
    cen = int(np.ceil(lo / CEN_STEP) * CEN_STEP)
    while cen < hi:
        for L in DEL_LENGTHS:
            tel = cen + L
            if loc.boundary[1] < tel < hi:
                tag = f"{cen / 1e6:.2f}–{tel / 1e6:.2f} Mb"
                out.append(Variant(f"DEL-B {tag}", cen, tel, "DEL-B"))
                out.append(Variant(f"DelBS {cen / 1e6:.2f}–{(loc.boundary[0] - 2_000) / 1e6:.3f} Mb", cen,
                                   loc.boundary[0] - 2_000, "DelBS"))
                out.append(Variant(f"DEL-B {tag}, frontière épargnée", cen, tel, "DEL-B, frontière épargnée",
                                   keep=(loc.boundary,)))
        cen += CEN_STEP
    # Les témoins DelBS d'une même extrémité centromérique sont identiques : un seul.
    seen, uniq = set(), []
    for v in out:
        k = (v.kind, v.start, v.end)
        if k not in seen:
            seen.add(k)
            uniq.append(v)
    return uniq


def central(loc: Locus) -> tuple[Variant, Variant, Variant]:
    """Le DEL-B du milieu de la famille (extrémité centromérique médiane, 1,8 Mb), et ses deux
    témoins : ceux-là seuls passent par OpenMM."""
    dels = [v for v in loc.variants if v.kind == "DEL-B"]
    cens = sorted({v.start for v in dels})
    cen = cens[len(cens) // 2]
    d = min((v for v in dels if v.start == cen), key=lambda v: abs(v.length - 1_800_000))
    s = next(v for v in loc.variants if v.kind == "DelBS" and v.start == cen)
    k = next(v for v in loc.variants if v.kind == "DEL-B, frontière épargnée" and v.start == d.start and v.end == d.end)
    return d, s, k
