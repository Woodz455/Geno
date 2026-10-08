"""Parseurs des formats 1D.

Une seule règle de conversion, appliquée ici et nulle part ailleurs :

    GTF/GFF  1-based, bornes incluses   →  start - 1, end
    BED      0-based, demi-ouvert        →  inchangé
    cytoBand 0-based, demi-ouvert        →  inchangé   (table UCSC, pas le format d'affichage)

Une erreur d'un pb sur cette conversion décale un exon d'une base et personne ne
le voit avant de lire la séquence au niveau 6. C'est pour ça qu'elle est isolée
et testée.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Iterator

_GTF_ATTR = re.compile(r'(\S+)\s+"([^"]*)"')


def _open(path: Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _rows(path: Path) -> Iterator[list[str]]:
    with _open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            yield line.rstrip("\n").split("\t")


def parse_gtf(path: Path, feature: str) -> Iterator[tuple[str, int, int, dict]]:
    """Lit un GTF et ne garde qu'un type de feature (gene, transcript, exon…)."""
    for f in _rows(path):
        if len(f) < 9 or f[2] != feature:
            continue
        attrs = dict(_GTF_ATTR.findall(f[8]))
        rec = {
            "name": attrs.get("gene_name") or attrs.get("gene_id", ""),
            "strand": f[6],
            "biotype": attrs.get("gene_type") or attrs.get("gene_biotype", ""),
        }
        for key in ("gene_id", "transcript_id", "exon_number"):
            if key in attrs:
                rec[key] = attrs[key]
        yield f[0], int(f[3]) - 1, int(f[4]), rec


def parse_bed(path: Path) -> Iterator[tuple[str, int, int, dict]]:
    """BED3 à BED6. Les colonnes absentes sont simplement omises."""
    for f in _rows(path):
        if len(f) < 3:
            continue
        rec: dict = {}
        if len(f) > 3 and f[3] != ".":
            rec["name"] = f[3]
        if len(f) > 4 and f[4] != ".":
            rec["score"] = _num(f[4])
        if len(f) > 5 and f[5] in "+-":
            rec["strand"] = f[5]
        yield f[0], int(f[1]), int(f[2]), rec


def parse_cytoband(path: Path) -> Iterator[tuple[str, int, int, dict]]:
    """Table cytoBandIdeo d'UCSC : chrom, start, end, nom de bande, coloration Giemsa."""
    for f in _rows(path):
        if len(f) < 5:
            continue
        chrom = f[0]
        band = f[3]
        yield chrom, int(f[1]), int(f[2]), {
            "name": f"{chrom.removeprefix('chr')}{band}",
            "band": band,
            "stain": f[4],
            "arm": band[0] if band[:1] in "pq" else "",
        }


def _num(value: str) -> float | int | str:
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


PARSERS = {
    "gtf": parse_gtf,
    "bed": parse_bed,
    "cytoband": parse_cytoband,
}


def read_track(spec: dict, base: Path) -> Iterator[tuple[str, int, int, dict]]:
    """Instancie le parseur décrit par une entrée de `tracks.json`."""
    fmt = spec["format"]
    parser = PARSERS.get(fmt)
    if parser is None:
        raise ValueError(f"format inconnu : {fmt} (connus : {', '.join(PARSERS)})")
    path = base / spec["file"]
    if fmt == "gtf":
        return parser(path, spec.get("feature", "gene"))
    return parser(path)
