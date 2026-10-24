"""Conversion de coordonnées par fichier de chaînes (format UCSC `.over.chain`).

Une chaîne aligne un intervalle de l'assemblage source sur un intervalle de l'assemblage
cible, par blocs sans trou : `taille dt dq`, puis la taille du dernier bloc. Lever une
position, c'est trouver le bloc qui la contient et reporter son décalage. Une chaîne sur le
brin − de la cible compte ses coordonnées depuis la fin du chromosome cible.

Deux usages cette semaine, une seule implémentation :
- hg38 → T2T-CHM13 (hs1), pour la fiche : la position T2T de ce qu'on clique ;
- GRCh37 → GRCh38 (chaîne Ensembl), pour Rao 2014 et le Repli-seq UW, publiés en hg19.

Plusieurs chaînes peuvent couvrir la même position source (duplications) : on garde celle de
plus haut score, comme `liftOver` par défaut. Une position hors de tout bloc ne se lève pas —
elle n'a pas d'équivalent dans l'autre assemblage, et on le dit au lieu d'interpoler.

Les noms Ensembl (`1`, `X`, `MT`) sont ramenés aux noms UCSC (`chr1`, `chrX`, `chrM`).
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

import numpy as np


def ucsc_name(name: str) -> str:
    if name.startswith("chr"):
        return name
    return "chrM" if name == "MT" else f"chr{name}"


def _open_text(path: str | Path):
    with open(path, "rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="ascii")
    return open(path, "r", encoding="ascii")


@dataclass(frozen=True)
class Lifted:
    chrom: str
    pos: int          # 0-based, sur le brin + de la cible
    strand: str       # brin de la chaîne : '-' = la région est retournée dans la cible
    chain: int        # identifiant de la chaîne dans le fichier


@dataclass(frozen=True)
class LiftedInterval:
    chrom: str
    start: int
    end: int
    strand: str
    mapped: float     # fraction des bases de la source couvertes par un bloc de la chaîne retenue
    split: bool       # les deux bouts tombent sur des chaînes ou des chromosomes différents


class Chain:
    """Les blocs d'un fichier de chaînes, indexés par chromosome source."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        blocks: dict[str, list[tuple[int, int, int, int, int, int]]] = {}
        self.targets: list[tuple[str, int, str]] = []   # par chaîne : (chrom cible, taille, brin)
        self.scores: list[int] = []
        cur = None
        with _open_text(path) as f:
            for line in f:
                if line.startswith("chain"):
                    c = line.split()
                    t_name, t_start = ucsc_name(c[2]), int(c[5])
                    q_name, q_size, q_strand, q_start = ucsc_name(c[7]), int(c[8]), c[9], int(c[10])
                    score = int(float(c[1]))
                    k = len(self.targets)
                    self.targets.append((q_name, q_size, q_strand))
                    self.scores.append(score)
                    cur = [t_name, t_start, q_start, k, score]
                    blocks.setdefault(t_name, [])
                    continue
                parts = line.split()
                if not parts or cur is None:
                    continue
                size = int(parts[0])
                t_name, t, q, k, score = cur
                blocks[t_name].append((t, t + size, q, k, score, 0))
                if len(parts) == 3:
                    cur[1] = t + size + int(parts[1])
                    cur[2] = q + size + int(parts[2])
                else:
                    cur = None
        self._idx: dict[str, tuple[np.ndarray, ...]] = {}
        for chrom, rows in blocks.items():
            a = np.array(rows, dtype=np.int64)
            o = np.argsort(a[:, 0], kind="stable")
            a = a[o]
            # Plus grande fin vue jusqu'ici : un bloc long peut recouvrir des blocs qui
            # commencent après lui (chaînes imbriquées). Sert à borner la recherche.
            self._idx[chrom] = (a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4], np.maximum.accumulate(a[:, 1]))

    def chroms(self) -> list[str]:
        return list(self._idx)

    def _hits(self, chrom: str, pos: int) -> np.ndarray:
        if chrom not in self._idx:
            return np.empty(0, dtype=np.int64)
        ts, te, _, _, _, emax = self._idx[chrom]
        hi = int(np.searchsorted(ts, pos, side="right"))
        lo = int(np.searchsorted(emax, pos, side="right"))     # avant lo, tout bloc finit ≤ pos
        if hi <= lo:
            return np.empty(0, dtype=np.int64)
        r = np.arange(lo, hi)
        return r[te[lo:hi] > pos]

    def point(self, chrom: str, pos: int) -> Lifted | None:
        """Position 0-based de la source → cible, ou None si elle n'est dans aucun bloc."""
        r = self._hits(chrom, pos)
        if r.size == 0:
            return None
        ts, _, qs, ks, scores, _ = self._idx[chrom]
        best = r[np.argmax(scores[r])]
        k = int(ks[best])
        q_name, q_size, q_strand = self.targets[k]
        q = int(qs[best] + (pos - ts[best]))
        if q_strand == "-":
            q = q_size - 1 - q
        return Lifted(q_name, q, q_strand, k)

    def points(self, chrom: str, pos: np.ndarray) -> tuple[list[str | None], np.ndarray, np.ndarray]:
        """Version vectorisée de `point` : (chromosomes, positions, brins ±1) ; −1 si non levée."""
        out_c: list[str | None] = []
        out_p = np.full(len(pos), -1, dtype=np.int64)
        out_s = np.zeros(len(pos), dtype=np.int8)
        for i, p in enumerate(pos):
            hit = self.point(chrom, int(p))
            out_c.append(hit.chrom if hit else None)
            if hit:
                out_p[i] = hit.pos
                out_s[i] = 1 if hit.strand == "+" else -1
        return out_c, out_p, out_s

    def interval(self, chrom: str, start: int, end: int) -> LiftedInterval | None:
        """Lève `[start, end)` par ses deux bouts ; la fraction couverte vient des blocs."""
        a = self._first_mapped(chrom, start, end, +1)
        b = self._first_mapped(chrom, start, end, -1)
        if a is None or b is None:
            return None
        lo, hi = sorted((a.pos, b.pos))
        split = a.chrom != b.chrom or a.chain != b.chain
        covered = self._covered(chrom, start, end, a.chain)
        return LiftedInterval(a.chrom, lo, hi + 1, a.strand, covered / max(1, end - start), split)

    def _first_mapped(self, chrom: str, start: int, end: int, step: int) -> Lifted | None:
        """Le premier point levable en partant d'un bout (les trous de chaîne sont courts)."""
        if chrom not in self._idx:
            return None
        ts, te, _, _, _, _ = self._idx[chrom]
        p = start if step > 0 else end - 1
        hit = self.point(chrom, p)
        if hit:
            return hit
        # Saute au bord du bloc suivant (ou précédent) au lieu d'avancer base par base.
        if step > 0:
            i = int(np.searchsorted(ts, p, side="right"))
            cand = ts[i:][ts[i:] < end]
            for c in cand[:64]:
                hit = self.point(chrom, int(c))
                if hit:
                    return hit
        else:
            i = int(np.searchsorted(ts, p, side="right"))
            cand = te[:i][(te[:i] > start) & (te[:i] <= p)]
            for c in sorted(cand, reverse=True)[:64]:
                hit = self.point(chrom, int(c) - 1)
                if hit:
                    return hit
        return None

    def _covered(self, chrom: str, start: int, end: int, chain: int) -> int:
        ts, te, _, ks, _, _ = self._idx[chrom]
        m = (ks == chain) & (ts < end) & (te > start)
        return int(np.sum(np.minimum(te[m], end) - np.maximum(ts[m], start)))

    def blocks(self, chrom: str, start: int, end: int) -> list[tuple[int, int, int, int]]:
        """Blocs (tStart, tEnd, qStart, chaîne) qui touchent `[start, end)` — pour l'embarquer."""
        if chrom not in self._idx:
            return []
        ts, te, qs, ks, _, _ = self._idx[chrom]
        m = (ts < end) & (te > start)
        return [(int(a), int(b), int(c), int(d)) for a, b, c, d in zip(ts[m], te[m], qs[m], ks[m])]
