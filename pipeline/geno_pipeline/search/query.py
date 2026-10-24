"""La grammaire d'une requête — la même que `parseQuery` en TypeScript.

Ce qu'on tape suit la convention de tout le monde : **1-based, bornes incluses**
(« chr17:43,044,295 », « chr7:5,527,000-5,530,600 »). Ce qui en sort est interne : 0-based,
demi-ouvert. Un chromosome s'écrit avec ou sans « chr », en majuscules ou non ; les séparateurs
de milliers (virgule, espace fine, souligné) sont ignorés.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

_RS = re.compile(r"^rs([0-9]{1,15})$", re.I)
_POS = re.compile(r"^(?:chr)?([0-9]{1,2}|x|y|m|mt)\s*[:\s]\s*([0-9,_   ]+?)(?:\s*[-–]\s*([0-9,_   ]+))?$", re.I)
_ENSG = re.compile(r"^(ENSG[0-9]{11})(?:\.[0-9]+)?$", re.I)
_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._@-]{0,39}$")
_DIGITS = re.compile(r"[,_   ]")


@dataclass(frozen=True)
class Query:
    kind: str                        # "rs" | "position" | "interval" | "gene-id" | "gene" | "malformed"
    key: str = ""                    # rs : le numéro ; gène : le nom en majuscules ; id : ENSG…
    chrom: str = ""
    start: int = 0                   # 0-based
    end: int = 0                     # demi-ouvert
    error: str = ""


def chrom_name(c: str) -> str:
    c = c.upper()
    return "chrM" if c in ("M", "MT") else f"chr{c}"


def parse(text: str) -> Query:
    t = text.strip()
    if not t:
        return Query("malformed", error="requête vide")
    m = _RS.match(t)
    if m:
        return Query("rs", key=str(int(m.group(1))))
    m = _POS.match(t)
    if m:
        chrom = chrom_name(m.group(1))
        a = int(_DIGITS.sub("", m.group(2)))
        if m.group(3) is None:
            if a < 1:
                return Query("malformed", error="les coordonnées 1-based commencent à 1")
            return Query("position", chrom=chrom, start=a - 1, end=a)
        b = int(_DIGITS.sub("", m.group(3)))
        if a < 1:
            return Query("malformed", error="les coordonnées 1-based commencent à 1")
        if b < a:
            return Query("malformed", error="fin avant début")
        return Query("interval", chrom=chrom, start=a - 1, end=b)
    m = _ENSG.match(t)
    if m:
        return Query("gene-id", key=m.group(1).upper())
    if _NAME.match(t) and not t.isdigit():
        return Query("gene", key=t.upper())
    return Query("malformed", error="ni un gène, ni une position chrom:pos, ni un rsID")
