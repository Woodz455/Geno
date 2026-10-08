"""Caryotype et longueurs de chromosomes — ce qu'on modélise, et en combien d'exemplaires.

GM12878 est une lignée lymphoblastoïde féminine de caryotype quasi normal,
**46,XX**. Le noyau est donc deux exemplaires de chacun des 22 autosomes plus
deux X. Les deux X ne sont pas équivalents : l'un est actif (Xa), l'autre inactif
(Xi), et leur repliement diffère profondément — l'Xi n'a quasiment pas de TADs et
se scinde en deux super-domaines autour de XIST (Rao 2014, Deng 2015). La semaine 6
ne modélise pas cette différence, mais elle **nomme les deux copies séparément**
pour que la semaine où on la modélisera n'ait pas à réécrire le modèle de données.

Sur la provenance des longueurs : une taille de chromosome est une propriété
déterministe de l'assemblage, pas une mesure. Elle reste néanmoins une donnée, et
la règle du projet est qu'une donnée remonte à un fichier empreint. D'où la
hiérarchie : si `data/core/hg38.chrom.sizes` existe (récupéré et verrouillé par
`make data-core`), c'est lui qui fait foi ; sinon on retombe sur la table interne,
et **tout ce qui en découle est estampillé `builtin`** jusqu'à ce que le réseau
s'ouvre. Ça ne se devine pas dans les sorties : c'est imprimé.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent

# GRCh38 / hg38, assemblage primaire. Transcrit ; à revérifier contre
# hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes dès que
# l'hôte est joignable — d'où l'estampille `builtin`.
HG38: dict[str, int] = {
    "chr1": 248_956_422,
    "chr2": 242_193_529,
    "chr3": 198_295_559,
    "chr4": 190_214_555,
    "chr5": 181_538_259,
    "chr6": 170_805_979,
    "chr7": 159_345_973,
    "chr8": 145_138_636,
    "chr9": 138_394_717,
    "chr10": 133_797_422,
    "chr11": 135_086_622,
    "chr12": 133_275_309,
    "chr13": 114_364_328,
    "chr14": 107_043_718,
    "chr15": 101_991_189,
    "chr16": 90_338_345,
    "chr17": 83_257_441,
    "chr18": 80_373_285,
    "chr19": 58_617_616,
    "chr20": 64_444_167,
    "chr21": 46_709_983,
    "chr22": 50_818_468,
    "chrX": 156_040_895,
    "chrY": 57_227_415,
}

# Bras courts acrocentriques : rDNA et satellites. Dans hg38 ils sont
# essentiellement des blocs de N, et *physiquement* ils sont dans le nucléole,
# pas dans la chromatine libre. On les marque pour que la semaine qui modélisera
# le nucléole sache où regarder — et pour que le rapport dise ce qu'il fait d'eux.
ACROCENTRIC = ("chr13", "chr14", "chr15", "chr21", "chr22")


@dataclass(frozen=True)
class Copy:
    """Un exemplaire physique d'un chromosome dans le noyau."""

    chrom: str
    homolog: int          # 0 ou 1 — les deux copies héritées
    length: int
    label: str            # "chr7:a", "chrX:Xi", …

    @property
    def acrocentric(self) -> bool:
        return self.chrom in ACROCENTRIC


@dataclass(frozen=True)
class Karyotype:
    """Le jeu complet de copies à modéliser, et d'où viennent ses longueurs."""

    name: str
    copies: tuple[Copy, ...]
    assembly: str
    provenance: str       # "builtin" ou le chemin du fichier verrouillé

    @property
    def n_copies(self) -> int:
        return len(self.copies)

    @property
    def total_bp(self) -> int:
        return sum(c.length for c in self.copies)

    def __str__(self) -> str:
        return (
            f"{self.name} · {self.n_copies} copies · {self.total_bp / 1e9:.2f} Gb · "
            f"{self.assembly} [{self.provenance}]"
        )


def chrom_sizes(path: Path | None = None) -> tuple[dict[str, int], str]:
    """Longueurs de chromosomes, du fichier verrouillé s'il existe, sinon de la table.

    Renvoie aussi la provenance, parce qu'un modèle construit sur la table interne
    et un modèle construit sur le fichier officiel ne se valent pas et ne doivent
    pas se ressembler dans les sorties.
    """
    path = path or ROOT / "data" / "core" / "hg38.chrom.sizes"
    if path.exists():
        sizes: dict[str, int] = {}
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            name, _, size = line.partition("\t")
            # Assemblage primaire seulement : pas de scaffolds ni d'haplotypes alt.
            if "_" in name:
                continue
            sizes[name] = int(size)
        if sizes:
            return sizes, str(path)
    return dict(HG38), "builtin"


def gm12878(path: Path | None = None) -> Karyotype:
    """46,XX — deux exemplaires de chr1–22, un Xa et un Xi, pas de Y."""
    sizes, provenance = chrom_sizes(path)
    copies: list[Copy] = []
    for i in range(1, 23):
        chrom = f"chr{i}"
        for h, tag in enumerate("ab"):
            copies.append(Copy(chrom, h, sizes[chrom], f"{chrom}:{tag}"))
    for h, tag in enumerate(("Xa", "Xi")):
        copies.append(Copy("chrX", h, sizes["chrX"], f"chrX:{tag}"))
    return Karyotype("GM12878 46,XX", tuple(copies), "GRCh38", provenance)
