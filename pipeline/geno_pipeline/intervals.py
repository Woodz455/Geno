"""Magasin d'intervalles génomiques.

Toutes les coordonnées internes sont **0-based, demi-ouvertes** — la convention BED.
Les parseurs convertissent à l'entrée ; rien dans ce module ne manipule du 1-based.

Index : tableau trié par début, plus un maximum de fin par bloc.

    Un intervalle [s, e) chevauche la requête [qs, qe) si et seulement si
    s < qe  ET  e > qs.

    La première condition se résout par recherche dichotomique sur `starts`.
    La seconde est le piège : les intervalles sont triés par début, pas par fin,
    donc un candidat peut se trouver arbitrairement loin en arrière. Un balayage
    naïf vers l'arrière est O(n) dès qu'un gène long traîne au début du
    chromosome — et il y en a : DMD fait 2,2 Mb, CNTNAP2 2,3 Mb.

    D'où le maximum par bloc : `blockmax[b]` est la plus grande fin du bloc b.
    Un bloc dont le maximum est ≤ qs ne peut contenir aucun chevauchement et se
    saute d'un seul test. Le filtrage des blocs candidats est une opération
    numpy vectorisée sur quelques milliers d'éléments, pas une boucle Python.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import numpy as np

BLOCK = 512
"""Taille de bloc de l'index. 512 tient dans quelques lignes de cache et donne
quelques milliers de blocs pour un chromosome dense, ce que numpy filtre en
microsecondes."""

COALESCE = 65_536
"""Trou maximal, en octets, toléré entre deux enregistrements d'attributs avant
de couper la lecture en deux. Lire 64 ko d'un coup coûte moins cher qu'un second
aller-retour ; au-delà, on lirait pour rien."""


@dataclass(frozen=True, slots=True)
class Feature:
    """Un intervalle et ses attributs."""

    chrom: str
    start: int
    end: int
    attrs: dict

    @property
    def length(self) -> int:
        return self.end - self.start

    def locus(self) -> str:
        return f"{self.chrom}:{self.start + 1:,}-{self.end:,}"


def _blockmax(ends: np.ndarray) -> np.ndarray:
    """Maximum de fin par bloc de BLOCK entrées."""
    n = ends.size
    if n == 0:
        return np.empty(0, dtype=np.int64)
    nb = (n + BLOCK - 1) // BLOCK
    pad = nb * BLOCK - n
    padded = np.concatenate([ends, np.full(pad, np.iinfo(np.int64).min, dtype=np.int64)])
    return padded.reshape(nb, BLOCK).max(axis=1)


def _overlaps(
    starts: np.ndarray, ends: np.ndarray, blockmax: np.ndarray, qs: int, qe: int
) -> np.ndarray:
    """Indices locaux des intervalles chevauchant [qs, qe), en ordre croissant."""
    empty = np.empty(0, dtype=np.int64)
    if starts.size == 0 or qe <= qs:
        return empty

    # Tous les candidats ont start < qe.
    hi = int(np.searchsorted(starts, qe, side="left"))
    if hi == 0:
        return empty

    # Parmi eux, seuls les blocs dont une fin dépasse qs peuvent contenir un hit.
    bh = (hi - 1) // BLOCK
    cand = np.flatnonzero(blockmax[: bh + 1] > qs)
    if cand.size == 0:
        return empty

    idx = (cand[:, None] * BLOCK + np.arange(BLOCK, dtype=np.int64)[None, :]).ravel()
    idx = idx[idx < hi]
    return idx[ends[idx] > qs]


class TrackWriter:
    """Écrit une piste sur disque : tableaux numpy memmappables + attributs JSONL."""

    def __init__(self, name: str, out_dir: Path):
        self.name = name
        self.dir = out_dir / name
        self.dir.mkdir(parents=True, exist_ok=True)

    def write(self, rows: Iterable[tuple[str, int, int, dict]]) -> dict:
        by_chrom: dict[str, list[tuple[int, int, dict]]] = {}
        for chrom, start, end, attrs in rows:
            if end <= start:
                raise ValueError(
                    f"{self.name}: intervalle vide ou inversé {chrom}:{start}-{end}"
                )
            by_chrom.setdefault(chrom, []).append((start, end, attrs))

        starts: list[int] = []
        ends: list[int] = []
        blockmaxes: list[np.ndarray] = []
        meta: dict[str, dict] = {}
        attr_offsets: list[int] = [0]

        blob = self.dir / "attrs.jsonl"
        pos = 0
        with blob.open("wb") as fh:
            for chrom in sorted(by_chrom, key=_chrom_key):
                items = sorted(by_chrom[chrom], key=lambda r: (r[0], r[1]))
                lo = len(starts)
                blo = sum(b.size for b in blockmaxes)

                ce = np.fromiter((e for _, e, _ in items), dtype=np.int64, count=len(items))
                starts.extend(s for s, _, _ in items)
                ends.extend(int(e) for e in ce)
                blockmaxes.append(_blockmax(ce))

                for _, _, attrs in items:
                    raw = (json.dumps(attrs, separators=(",", ":"), ensure_ascii=False) + "\n").encode()
                    fh.write(raw)
                    pos += len(raw)
                    attr_offsets.append(pos)

                meta[chrom] = {
                    "lo": lo,
                    "hi": len(starts),
                    "blo": blo,
                    "bhi": blo + blockmaxes[-1].size,
                }

        np.save(self.dir / "starts.npy", np.asarray(starts, dtype=np.int64))
        np.save(self.dir / "ends.npy", np.asarray(ends, dtype=np.int64))
        np.save(
            self.dir / "blockmax.npy",
            np.concatenate(blockmaxes) if blockmaxes else np.empty(0, dtype=np.int64),
        )
        np.save(self.dir / "attrs.off.npy", np.asarray(attr_offsets, dtype=np.int64))

        info = {"track": self.name, "n": len(starts), "chroms": meta}
        (self.dir / "meta.json").write_text(json.dumps(info, indent=2), encoding="utf-8")
        return info


class Track:
    """Une piste chargée. Les tableaux sont memmappés ; les attributs sont lus à la demande."""

    def __init__(self, path: Path):
        self.path = path
        self.meta = json.loads((path / "meta.json").read_text(encoding="utf-8"))
        self.name: str = self.meta["track"]
        self.starts = np.load(path / "starts.npy", mmap_mode="r")
        self.ends = np.load(path / "ends.npy", mmap_mode="r")
        self.blockmax = np.load(path / "blockmax.npy", mmap_mode="r")
        self.offsets = np.load(path / "attrs.off.npy", mmap_mode="r")
        self._blob = (path / "attrs.jsonl").open("rb")

    def __len__(self) -> int:
        return int(self.meta["n"])

    @property
    def chroms(self) -> list[str]:
        return list(self.meta["chroms"])

    def close(self) -> None:
        self._blob.close()

    def _attrs(self, i: int) -> dict:
        a, b = int(self.offsets[i]), int(self.offsets[i + 1])
        self._blob.seek(a)
        return json.loads(self._blob.read(b - a))

    def _attrs_many(self, rows: np.ndarray) -> list[dict]:
        """Lit les attributs de plusieurs intervalles en coalesçant les lectures.

        Les indices sortent de l'index déjà triés, donc leurs enregistrements
        sont contigus ou presque dans le blob. Un `seek` + `read` par hit coûte
        plus cher que la recherche elle-même dès qu'une requête ramène une
        centaine de features. On regroupe donc les indices séparés par moins de
        COALESCE octets en une seule lecture.

        La borne est ce qui empêche le remède d'être pire que le mal : deux hits
        aux extrémités d'une piste d'un million d'entrées ne doivent pas
        déclencher la lecture du blob entier.
        """
        if rows.size == 0:
            return []

        starts = self.offsets[rows].astype(np.int64)
        stops = self.offsets[rows + 1].astype(np.int64)
        # Nouvelle lecture dès que le trou depuis la fin précédente dépasse COALESCE.
        cuts = np.flatnonzero(starts[1:] - stops[:-1] > COALESCE) + 1
        groups = np.split(np.arange(rows.size), cuts)

        out: list[dict] = []
        for g in groups:
            a, b = int(starts[g[0]]), int(stops[g[-1]])
            self._blob.seek(a)
            buf = self._blob.read(b - a)
            for k in g:
                out.append(json.loads(buf[int(starts[k]) - a : int(stops[k]) - a]))
        return out

    def query(self, chrom: str, start: int, end: int) -> list[Feature]:
        c = self.meta["chroms"].get(chrom)
        if c is None:
            return []
        lo, hi, blo, bhi = c["lo"], c["hi"], c["blo"], c["bhi"]
        idx = _overlaps(
            self.starts[lo:hi], self.ends[lo:hi], self.blockmax[blo:bhi], start, end
        )
        rows = idx + lo
        return [
            Feature(chrom, int(self.starts[r]), int(self.ends[r]), attrs)
            for r, attrs in zip(rows.tolist(), self._attrs_many(rows))
        ]

    def count(self, chrom: str, start: int, end: int) -> int:
        """Compte sans lire les attributs — le chemin chaud du benchmark."""
        c = self.meta["chroms"].get(chrom)
        if c is None:
            return 0
        lo, hi, blo, bhi = c["lo"], c["hi"], c["blo"], c["bhi"]
        return int(
            _overlaps(
                self.starts[lo:hi], self.ends[lo:hi], self.blockmax[blo:bhi], start, end
            ).size
        )


class Store:
    """L'ensemble des pistes, plus la provenance."""

    def __init__(self, path: Path):
        self.path = Path(path)
        index = self.path / "index.json"
        if not index.exists():
            raise FileNotFoundError(
                f"aucun magasin en {self.path} — lancer `geno build` d'abord"
            )
        self.index = json.loads(index.read_text(encoding="utf-8"))
        self.tracks: dict[str, Track] = {
            name: Track(self.path / name) for name in self.index["tracks"]
        }

    @property
    def assembly(self) -> str:
        return self.index.get("assembly", "?")

    @property
    def source(self) -> str:
        return self.index.get("source", "?")

    def close(self) -> None:
        for t in self.tracks.values():
            t.close()

    def query(self, chrom: str, start: int, end: int, tracks: Sequence[str] | None = None):
        names = list(tracks) if tracks else list(self.tracks)
        for name in names:
            track = self.tracks.get(name)
            if track is None:
                raise KeyError(f"piste inconnue : {name} (connues : {', '.join(self.tracks)})")
            yield name, track.query(chrom, start, end)


def _chrom_key(chrom: str) -> tuple[int, str]:
    """Ordre naturel : chr1 … chr22, chrX, chrY, chrM, puis le reste."""
    bare = chrom[3:] if chrom.startswith("chr") else chrom
    if bare.isdigit():
        return (int(bare), "")
    return ({"X": 100, "Y": 101, "M": 102, "MT": 102}.get(bare, 200), bare)


def write_store(
    out: Path, tracks: dict[str, Iterator[tuple[str, int, int, dict]]], meta: dict
) -> dict:
    """Écrit un magasin complet et son index."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    counts = {name: TrackWriter(name, out).write(rows)["n"] for name, rows in tracks.items()}
    index = {**meta, "tracks": counts}
    (out / "index.json").write_text(json.dumps(index, indent=2, ensure_ascii=False), encoding="utf-8")
    return index
