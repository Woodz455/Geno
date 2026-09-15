"""Composition d'une fiche de locus.

C'est l'ancêtre direct de la fiche d'information de la semaine 13 : même question
posée au même magasin, rendue en texte au lieu d'un panneau. Quand le viewer
arrivera, il appellera `report()` et rendra le même dictionnaire.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Sequence

from .intervals import Feature, Store

_REGION = re.compile(r"^(?P<chrom>[\w.]+)(?::(?P<start>[\d,_]+)-(?P<end>[\d,_]+))?$")

DISPLAY_ORDER = ["cytoband", "genes", "exons", "ccre", "chromhmm", "ctcf"]
LABELS = {
    "cytoband": "bande cytogénétique",
    "genes": "gènes",
    "exons": "exons",
    "ccre": "éléments cis-régulateurs",
    "chromhmm": "état chromatinien",
    "ctcf": "sites CTCF",
}


class RegionError(ValueError):
    pass


def parse_region(text: str, chrom_sizes: dict[str, int] | None = None) -> tuple[str, int, int]:
    """« chr7:5,527,000-5,530,600 » → ('chr7', 5526999, 5530600).

    L'entrée suit la convention que tout le monde tape : **1-based, bornes
    incluses**. La sortie est interne : 0-based, demi-ouverte.
    """
    m = _REGION.match(text.strip())
    if not m:
        raise RegionError(f"région illisible : {text!r} — attendu chrom:début-fin")
    chrom = m.group("chrom")
    if m.group("start") is None:
        size = (chrom_sizes or {}).get(chrom)
        if size is None:
            raise RegionError(f"{chrom} sans bornes : taille de chromosome inconnue")
        return chrom, 0, size
    start = int(m.group("start").replace(",", "").replace("_", ""))
    end = int(m.group("end").replace(",", "").replace("_", ""))
    if start < 1:
        raise RegionError("les coordonnées 1-based commencent à 1")
    if end < start:
        raise RegionError(f"fin ({end:,}) avant début ({start:,})")
    return chrom, start - 1, end


@dataclass
class LocusReport:
    chrom: str
    start: int
    end: int
    assembly: str
    source: str
    tracks: dict[str, list[Feature]] = field(default_factory=dict)
    flank: int = 0
    """Marge appliquée de part et d'autre de la région demandée.

    Elle n'est pas cosmétique : une boucle CTCF encadre son gène par
    construction, donc ses deux ancres tombent hors des bornes du gène. Sans
    marge, la piste CTCF d'un locus est vide alors que la boucle existe.
    """

    @property
    def span(self) -> int:
        return self.end - self.start

    def locus(self) -> str:
        return f"{self.chrom}:{self.start + 1:,}-{self.end:,}"

    def is_empty(self) -> bool:
        return not any(self.tracks.values())

    def to_dict(self) -> dict:
        return {
            "region": {
                "chrom": self.chrom,
                "start": self.start,
                "end": self.end,
                "display": self.locus(),
                "span": self.span,
            },
            "assembly": self.assembly,
            "source": self.source,
            "flank": self.flank,
            "tracks": {
                name: [
                    {"chrom": f.chrom, "start": f.start, "end": f.end, **f.attrs}
                    for f in feats
                ]
                for name, feats in self.tracks.items()
            },
        }


def report(
    store: Store, region: str, tracks: Sequence[str] | None = None, flank: int = 0
) -> LocusReport:
    if flank < 0:
        raise RegionError("la marge ne peut pas être négative")
    chrom, start, end = parse_region(region)
    lo, hi = max(0, start - flank), end + flank
    found = dict(store.query(chrom, lo, hi, tracks))
    ordered = {
        name: found[name]
        for name in sorted(found, key=lambda n: (DISPLAY_ORDER.index(n) if n in DISPLAY_ORDER else 99, n))
    }
    return LocusReport(chrom, start, end, store.assembly, store.source, ordered, flank)


def _describe(name: str, f: Feature) -> str:
    a = f.attrs
    strand = a.get("strand", "")
    if name == "cytoband":
        return f"{a.get('name', '?')}  {a.get('stain', '')}"
    if name in ("genes", "exons"):
        bits = [a.get("name", "?")]
        if strand:
            bits.append(strand)
        if a.get("biotype"):
            bits.append(a["biotype"])
        if a.get("exon_number"):
            bits.append(f"exon {a['exon_number']}")
        return "  ".join(bits)
    label = a.get("name", "?")
    if strand:
        label += f"  {strand}"
    if "score" in a:
        label += f"  score {a['score']}"
    return label


def render(rep: LocusReport) -> str:
    """Rendu texte de la fiche."""
    head = f"{rep.locus()}   {rep.span:,} pb   {rep.assembly}"
    if rep.flank:
        head += f"   ± {rep.flank:,} pb"
    if rep.source and rep.source != "?":
        head += f"   [{rep.source}]"
    out = [head, ""]

    if rep.is_empty():
        out.append("  aucune annotation sur cette région")
        return "\n".join(out)

    for name, feats in rep.tracks.items():
        label = LABELS.get(name, name)
        if not feats:
            out.append(f"  {label}  (0)")
            continue
        out.append(f"  {label}  ({len(feats)})")
        for f in feats:
            desc = _describe(name, f)
            out.append(f"    {desc:<44} {f.locus():>26}  {f.length:>9,} pb")
            gid = f.attrs.get("gene_id")
            if gid and name == "genes":
                out.append(f"    {gid}")
        out.append("")
    return "\n".join(out).rstrip() + "\n"
