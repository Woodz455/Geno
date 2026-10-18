"""Séquence de référence et annotation ponctuelle — le niveau 6 et sa cible.

Deux lecteurs, aucun dépendance nouvelle :

- **`.2bit`** (format UCSC) : en-tête, index des séquences, puis pour chacune ses blocs
  N, ses blocs masqués et l'ADN empaqueté à 2 bits par base (T=0, C=1, A=2, G=3, bits de
  poids fort d'abord). On ne lit que les octets de l'intervalle demandé : extraire 4 Mb
  d'un fichier de 800 Mo coûte 1 Mo de lecture. Le masquage des répétitions (minuscules)
  est ignoré : la double hélice n'en a pas l'usage, et une base est une base.
- **GTF GENCODE**, filtré sur un transcrit : ses CDS, dans l'ordre du transcrit, donnent la
  position génomique de n'importe quelle base codante. C'est ce qui rend la base visée par
  la plongée vérifiable : « c.818 de TP53 » est une affirmation sur l'annotation, pas une
  coordonnée recopiée.
"""

from __future__ import annotations

import gzip
import re
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

SIGNATURE = 0x1A412743
_BASES = "TCAG"
_LETTERS = np.frombuffer(_BASES.encode("ascii"), dtype=np.uint8)


class TwoBit:
    """Lecteur `.2bit` à accès direct. Version 0 (offsets 32 bits) et 1 (64 bits)."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self._f = open(self.path, "rb")
        head = self._f.read(16)
        sig = struct.unpack("<I", head[:4])[0]
        if sig == SIGNATURE:
            self._e = "<"
        elif struct.unpack(">I", head[:4])[0] == SIGNATURE:
            self._e = ">"
        else:
            raise ValueError(f"{path} n'est pas un fichier .2bit")
        version, count, _ = struct.unpack(self._e + "III", head[4:16])
        if version not in (0, 1):
            raise ValueError(f"version .2bit inconnue : {version}")
        off_fmt = "Q" if version == 1 else "I"
        self.offsets: dict[str, int] = {}
        for _ in range(count):
            (n,) = self._f.read(1)
            name = self._f.read(n).decode("ascii")
            (off,) = struct.unpack(self._e + off_fmt, self._f.read(struct.calcsize(off_fmt)))
            self.offsets[name] = off
        self._records: dict[str, tuple[int, list[tuple[int, int]], int]] = {}

    def close(self) -> None:
        self._f.close()

    def __enter__(self) -> TwoBit:
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def _record(self, chrom: str) -> tuple[int, list[tuple[int, int]], int]:
        """(longueur, blocs N [(début, fin)], offset de l'ADN empaqueté)."""
        if chrom not in self._records:
            if chrom not in self.offsets:
                raise KeyError(f"{chrom} absent de {self.path.name}")
            f, e = self._f, self._e
            f.seek(self.offsets[chrom])
            size, n_n = struct.unpack(e + "II", f.read(8))
            starts = struct.unpack(e + f"{n_n}I", f.read(4 * n_n))
            sizes = struct.unpack(e + f"{n_n}I", f.read(4 * n_n))
            (n_mask,) = struct.unpack(e + "I", f.read(4))
            f.seek(8 * n_mask + 4, 1)          # blocs masqués, puis le champ réservé
            self._records[chrom] = (size, [(s, s + l) for s, l in zip(starts, sizes)], f.tell())
        return self._records[chrom]

    def length(self, chrom: str) -> int:
        return self._record(chrom)[0]

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Bases de `[start, end)` (0-based, demi-ouvert), en majuscules ; N dans les trous."""
        size, nblocks, dna = self._record(chrom)
        if not 0 <= start <= end <= size:
            raise ValueError(f"{chrom}:{start}-{end} hors de [0, {size}]")
        first, last = start // 4, (end + 3) // 4
        self._f.seek(dna + first)
        raw = self._f.read(last - first)
        b = np.frombuffer(raw, dtype=np.uint8)
        codes = np.stack([b >> 6, (b >> 4) & 3, (b >> 2) & 3, b & 3], axis=1).reshape(-1)
        chars = _LETTERS[codes[start - 4 * first: end - 4 * first]]
        for a, z in nblocks:
            lo, hi = max(a, start), min(z, end)
            if hi > lo:
                chars[lo - start: hi - start] = ord("N")
        return chars.tobytes().decode("ascii")


def write_twobit(path: str | Path, seqs: dict[str, str]) -> None:
    """Écrit un `.2bit` version 0 — pour les tests, et pour rien d'autre."""
    names = list(seqs)
    index_size = sum(1 + len(n) + 4 for n in names)
    offset = 16 + index_size
    body = b""
    index = b""
    for name in names:
        s = seqs[name].upper()
        runs = [(m.start(), m.end() - m.start()) for m in re.finditer("N+", s)]
        rec = struct.pack("<II", len(s), len(runs))
        rec += struct.pack(f"<{len(runs)}I", *[a for a, _ in runs])
        rec += struct.pack(f"<{len(runs)}I", *[l for _, l in runs])
        rec += struct.pack("<II", 0, 0)
        packed = bytearray((len(s) + 3) // 4)
        for k, ch in enumerate(s):
            v = _BASES.index(ch) if ch in _BASES else 0
            packed[k // 4] |= v << (6 - 2 * (k % 4))
        rec += bytes(packed)
        index += bytes([len(name)]) + name.encode("ascii") + struct.pack("<I", offset + len(body))
        body += rec
    Path(path).write_bytes(struct.pack("<IIII", SIGNATURE, 0, len(names), 0) + index + body)


# ---------------------------------------------------------------------------
# GENCODE
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Transcript:
    gene: str
    transcript: str
    chrom: str
    strand: str
    tags: tuple[str, ...]
    exons: tuple[tuple[int, int], ...]     # 1-based, fermés, dans l'ordre génomique
    cds: tuple[tuple[int, int], ...]

    @property
    def cds_length(self) -> int:
        return sum(b - a + 1 for a, b in self.cds)

    def coding_position(self, c: int) -> int:
        """Coordonnée génomique 1-based de c.`c` (1 = le A de l'ATG), selon le brin."""
        if not 1 <= c <= self.cds_length:
            raise ValueError(f"c.{c} hors du CDS ({self.cds_length} pb)")
        k = c
        blocks = self.cds if self.strand == "+" else tuple(reversed(self.cds))
        for a, b in blocks:
            n = b - a + 1
            if k <= n:
                return a + k - 1 if self.strand == "+" else b - k + 1
            k -= n
        raise AssertionError("inatteignable")

    def codon(self, number: int) -> tuple[int, int, int]:
        """Les trois positions génomiques du codon `number`, dans le sens du transcrit."""
        return tuple(self.coding_position(3 * (number - 1) + i) for i in (1, 2, 3))  # type: ignore[return-value]


_ATTR = re.compile(r'(\S+) "([^"]*)"')


def read_transcript(gtf: str | Path, transcript: str) -> Transcript:
    """Un transcrit d'un GTF GENCODE (gzippé ou non), par identifiant sans version."""
    opener = gzip.open if open(gtf, "rb").read(2) == b"\x1f\x8b" else open
    exons: list[tuple[int, int]] = []
    cds: list[tuple[int, int]] = []
    meta: dict[str, object] = {}
    key = f'transcript_id "{transcript}'
    with opener(gtf, "rt", encoding="utf-8") as f:
        for line in f:
            if key not in line or line.startswith("#"):
                continue
            col = line.rstrip("\n").split("\t")
            attrs = _ATTR.findall(col[8])
            tid = next(v for k, v in attrs if k == "transcript_id")
            if tid.split(".")[0] != transcript:
                continue
            if col[2] == "transcript":
                meta = dict(
                    gene=next(v for k, v in attrs if k == "gene_name"),
                    transcript=tid,
                    chrom=col[0],
                    strand=col[6],
                    tags=tuple(v for k, v in attrs if k == "tag"),
                )
            elif col[2] == "exon":
                exons.append((int(col[3]), int(col[4])))
            elif col[2] == "CDS":
                cds.append((int(col[3]), int(col[4])))
    if not meta:
        raise KeyError(f"{transcript} absent de {gtf}")
    return Transcript(exons=tuple(sorted(exons)), cds=tuple(sorted(cds)), **meta)  # type: ignore[arg-type]


_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(s: str) -> str:
    return s.translate(_COMPLEMENT)[::-1]


def main(argv: list[str] | None = None) -> int:
    """`python -m geno_pipeline.sequence <fichier.2bit> <chrom> <début> <fin>` : les bases de
    `[début, fin)`, 0-based. Pour les outils qui comparent ce que le navigateur dessine à ce
    que le génome dit, par un autre chemin que le `.g3d`."""
    import sys

    a = sys.argv[1:] if argv is None else argv
    if len(a) != 4:
        print(main.__doc__, file=sys.stderr)
        return 2
    with TwoBit(a[0]) as tb:
        print(tb.fetch(a[1], int(a[2]), int(a[3])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
