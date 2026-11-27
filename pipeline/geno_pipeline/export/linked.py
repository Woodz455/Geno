"""Ce que montrent les vues liées (semaine 14) : la carte de contacts du modèle, et les
pistes linéaires d'igv.js.

**La carte est celle du modèle.** HiGlass montre les contacts que les 200 repliements du
noyau (semaine 7) et les simulations fines (semaine 8) produisent, pas du Hi-C mesuré — le
même choix que pour la 3D, qui reste simulée tant que le passage aux données réelles n'est
pas fait. Le fichier le dit dans chaque jeu de tuiles (`evidence: simule`).

**Deux jeux de tuiles**, que HiGlass superpose :

- `noyau` : 24 Mb → 750 kb, tout le génome, contacts cis et trans ; carte haploïde (les deux
  homologues tombent dans les mêmes bins, comme en séquençage — `ensemble.compartments`) ;
- `fin` : 512 kb → 2 kb, seulement dans les régions fines (TP53, ACTB), cis ; NaN ailleurs,
  pour que HiGlass laisse voir le noyau dessous.

**La grille de HiGlass est uniforme** en coordonnées absolues (chromosomes mis bout à bout),
et les billes du modèle ne le sont pas : 747 à 753 kb selon le chromosome, et un
chromosome ne commence pas sur un multiple de la résolution. Chaque bille est donc
**répartie** sur les pixels qu'elle recouvre, au prorata de sa longueur : le contact d'une
paire de billes se partage entre les pixels de leurs deux intervalles, et le total se
conserve à chaque résolution — c'est ce que `check` vérifie.

**Valeur** : contacts par repliement, au sens de la semaine 7 (centres à moins de 1,5 ×
(r_i + r_j)) — la somme, sur un pixel, des fréquences de contact des paires qu'il couvre.
Les deux jeux ont la même unité, pas la même physique : le raccord S8 échoue, et la
frontière entre les deux cartes le montrera.

Format `geno-contacts/1`, du même moule que le `.g3d` et le fichier de fiches : préambule
de 64 octets, en-tête JSON deflate, puis une tuile par bloc — 256 × 256 float16, octets
regroupés par rang (`shuffle`), deflate, empreinte. Une tuile `z.x.y` : colonnes = bins x,
lignes = bins y, comme HiGlass la dessine ; seules les tuiles x ≤ y sont écrites, HiGlass
produit les autres par miroir.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import sparse

from .g3d import _deflate, _inflate, _shuffle, _unshuffle, digest

MAGIC = b"\x89GCT\r\n\x1a\n"
VERSION = 1
PREAMBLE = 64
TILE = 256
NUCLEUS_RES = (24_000_000, 12_000_000, 6_000_000, 3_000_000, 1_500_000, 750_000)
FINE_RES = (512_000, 256_000, 128_000, 64_000, 32_000, 16_000, 8_000, 4_000, 2_000)
CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X"])


# ---------------------------------------------------------------------------
# Géométrie : coordonnées absolues et répartition sur la grille
# ---------------------------------------------------------------------------

def genome(sizes: dict[str, int], chroms: tuple[str, ...] = CHROMS) -> tuple[dict[str, int], int]:
    """Décalage absolu de chaque chromosome, dans l'ordre de `chroms`, et la longueur totale."""
    off, at = {}, 0
    for c in chroms:
        off[c] = at
        at += int(sizes[c])
    return off, at


def spread(starts: np.ndarray, ends: np.ndarray, res: int, n_px: int) -> sparse.csr_matrix:
    """(n_bins × n_px) : la part de chaque bin qui tombe dans chaque pixel. Chaque ligne
    somme à 1 — un contact réparti ne se perd ni ne se crée."""
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    rows, cols, vals = [], [], []
    first = starts // res
    last = (ends - 1) // res
    for k in range(int((last - first).max()) + 1 if len(starts) else 0):
        p = first + k
        ok = p <= last
        lo = np.maximum(starts[ok], p[ok] * res)
        hi = np.minimum(ends[ok], (p[ok] + 1) * res)
        rows.append(np.flatnonzero(ok))
        cols.append(p[ok])
        vals.append((hi - lo) / (ends[ok] - starts[ok]))
    if not rows:
        return sparse.csr_matrix((0, n_px))
    return sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(len(starts), n_px))


def to_pixels(m: sparse.spmatrix, s: sparse.csr_matrix) -> sparse.csr_matrix:
    """La matrice de bins (symétrique) sur la grille de pixels : Sᵀ M S."""
    return (s.T @ sparse.csr_matrix(m) @ s).tocsr()


def symmetric(i: np.ndarray, j: np.ndarray, v: np.ndarray, n: int) -> sparse.csr_matrix:
    """Pixels (i ≤ j) → matrice symétrique, la diagonale une seule fois."""
    up = sparse.coo_matrix((v, (i, j)), shape=(n, n)).tocsr()
    return (up + up.T - sparse.diags(up.diagonal())).tocsr()


# ---------------------------------------------------------------------------
# Les deux sources
# ---------------------------------------------------------------------------

@dataclass
class Source:
    """Une carte en bins du modèle, avant la grille : intervalles absolus et contacts."""

    starts: np.ndarray          # absolus, bp
    ends: np.ndarray
    matrix: sparse.csr_matrix   # contacts par repliement, symétrique
    structures: int
    cutoff: float
    note: str


def nucleus(ensemble_path: Path, offsets: dict[str, int], cutoff: float = 1.5) -> Source:
    """La carte haploïde de l'ensemble, avec les appelants de la semaine 7."""
    from ..ensemble.compartments import count_contacts, haploid_bins
    from ..ensemble.store import read

    e = read(ensemble_path)
    b = e.beads
    bin_of, bins, _ = haploid_bins(tuple(b.labels), b.copy_id, b.start, b.end)
    px = count_contacts(e.coords, b.radius, bin_of, len(bins), cutoff)
    m = symmetric(px["bin1_id"].to_numpy(), px["bin2_id"].to_numpy(),
                  px["count"].to_numpy(np.float64) / e.n_structures, len(bins))
    off = np.array([offsets[c] for c in bins["chrom"]], dtype=np.int64)
    return Source(off + bins["start"].to_numpy(np.int64), off + bins["end"].to_numpy(np.int64), m,
                  e.n_structures, cutoff, f"{e.n_structures} repliements du noyau (semaine 7), "
                  f"seuil {cutoff:g} × (r_i + r_j), homologues confondus")


def fine(npz: Path, offsets: dict[str, int]) -> tuple[Source, dict]:
    """La carte d'une région fine, monomère par monomère (convention de la semaine 8)."""
    from ..fine.observe import contact_map
    from ..fine.run import load

    coords, _, _, meta = load(str(npz))
    cutoff = float(meta.get("cutoff", 1.5))
    cm = contact_map(coords, cutoff=cutoff, bp_per_bead=int(meta["bp_per_bead"]))
    n = cm.m
    s0 = offsets[meta["chrom"]] + int(meta["start"])
    starts = s0 + np.arange(n, dtype=np.int64) * int(meta["bp_per_bead"])
    region = {"name": Path(npz).stem, "chrom": meta["chrom"], "start": int(meta["start"]),
              "end": int(meta["end"]), "abs_start": int(starts[0]), "abs_end": int(starts[-1] + meta["bp_per_bead"])}
    return Source(starts, starts + int(meta["bp_per_bead"]), sparse.csr_matrix(cm.counts), len(coords), cutoff,
                  f"{len(coords)} conformations de la région {meta['chrom']}:{meta['start']:,}–{meta['end']:,} "
                  f"(semaine 8, CTCF {meta.get('ctcf_source', '?')})"), region


# ---------------------------------------------------------------------------
# Écriture
# ---------------------------------------------------------------------------

def encode_tile(dense: np.ndarray) -> bytes:
    """float16 : la moitié du float32 après deflate (21 Mo au lieu de 43 pour le fichier
    entier), 11 bits de précision — une couleur de carte n'en lit pas plus, et `check`
    compare donc à 10⁻³ près, pas au bit."""
    return _deflate(_shuffle(np.ascontiguousarray(dense, dtype="<f2").tobytes(), 2))


def decode_tile(data: bytes) -> np.ndarray:
    return np.frombuffer(_unshuffle(_inflate(data), 2), "<f2").astype(np.float32).reshape(TILE, TILE)


def tiles(px: sparse.csr_matrix, n_px: int, masks: list[np.ndarray] | None = None):
    """Tuiles (x ≤ y) non vides : `dense[i, j] = px[y·256 + i, x·256 + j]`. `masks`, si
    donné (un masque (n_px,) par région) : un pixel n'a de valeur que si ses **deux**
    positions sont dans la même région — NaN sinon, y compris entre deux régions fines,
    que rien n'a simulées ensemble."""
    nt = -(-n_px // TILE)
    coo = px.tocoo()
    keep = coo.row >= coo.col            # ligne = y ≥ colonne = x : la moitié que HiGlass demande
    ty, tx = coo.row[keep] // TILE, coo.col[keep] // TILE
    present = set(zip(tx.tolist(), ty.tolist()))
    for x, y in sorted(present):
        if x > y or x >= nt or y >= nt:
            continue
        r0, c0 = y * TILE, x * TILE
        block = px[r0:min(r0 + TILE, n_px), c0:min(c0 + TILE, n_px)].toarray()
        dense = np.zeros((TILE, TILE), dtype=np.float32)
        dense[:block.shape[0], :block.shape[1]] = block
        if masks is not None:
            ok = np.zeros((TILE, TILE), bool)
            for m in masks:
                mr = np.zeros(TILE, bool)
                mc = np.zeros(TILE, bool)
                mr[:block.shape[0]] = m[r0:r0 + block.shape[0]]
                mc[:block.shape[1]] = m[c0:c0 + block.shape[1]]
                ok |= np.outer(mr, mc)
            dense[~ok] = np.nan
        yield x, y, dense


def tileset(uid: str, label: str, sources: list[Source], resolutions: tuple[int, ...], total: int,
            put, regions: list[dict] | None = None) -> dict:
    """Un jeu de tuiles : chaque résolution, ses tuiles non vides, et le total qu'elles portent."""
    levels = []
    for z, res in enumerate(resolutions):
        n_px = -(-total // res)
        px = None
        masks = [] if regions is not None else None
        for src in sources:
            s = spread(src.starts, src.ends, res, n_px)
            part = to_pixels(src.matrix, s)
            px = part if px is None else px + part
            if masks is not None:
                masks.append(np.asarray(s.sum(axis=0)).ravel() > 0)
        specs = {}
        for x, y, dense in tiles(px, n_px, masks):
            specs[f"{z}.{x}.{y}"] = put(encode_tile(dense))
        levels.append({"z": z, "resolution": res, "pixels": n_px, "tiles": specs,
                       "total": float(px.sum()), "max": float(px.max()) if px.nnz else 0.0})
    return {"uid": uid, "label": label, "evidence": "simule", "resolutions": list(resolutions),
            "bins_per_dimension": TILE, "min_pos": [0, 0], "max_pos": [total, total],
            "value": "contacts par repliement (somme des fréquences de contact des paires couvertes)",
            "sources": [{"note": s.note, "structures": s.structures, "cutoff": s.cutoff,
                         "bins": int(len(s.starts)), "total": float(s.matrix.sum())} for s in sources],
            "regions": regions, "levels": levels}


def build(ensemble_path: Path, fine_paths: list[Path], sizes: dict[str, int], out: Path,
          provenance: dict, log=print) -> dict:
    offsets, total = genome(sizes)
    blobs: list[bytes] = []
    offset = 0

    def put(b: bytes) -> dict:
        nonlocal offset
        spec = {"offset": offset, "length": len(b), "sha256": digest(b)}
        blobs.append(b)
        offset += len(b)
        return spec

    nuc = nucleus(ensemble_path, offsets)
    log(f"  noyau : {len(nuc.starts):,} bins, {nuc.matrix.nnz:,} pixels non nuls")
    sets = [tileset("noyau", "contacts du modèle de noyau — simulés", [nuc], NUCLEUS_RES, total, put)]
    fines, regions = [], []
    for p in fine_paths:
        src, reg = fine(p, offsets)
        fines.append(src)
        regions.append(reg)
        log(f"  région fine {reg['name']} : {len(src.starts):,} monomères")
    if fines:
        sets.append(tileset("fin", "contacts des régions fines — simulés", fines, FINE_RES, total, put, regions))
    header = {
        "format": "geno-contacts/1", "version": VERSION, "assembly": "GRCh38", "cell_type": "GM12878",
        "chromsizes": [[c, int(sizes[c])] for c in CHROMS], "tilesets": sets,
        "tile": {"shape": [TILE, TILE], "dtype": "<f2", "filters": ["shuffle", "deflate"],
                 "layout": "dense[i·256 + j] = M[y·256 + i][x·256 + j] ; tuiles x ≤ y"},
        "provenance": provenance,
    }
    hbytes = _deflate(json.dumps(header, ensure_ascii=False, separators=(",", ":")).encode())
    pre = struct.pack("<8sII32s", MAGIC, VERSION, len(hbytes), hashlib.sha256(hbytes).digest()).ljust(PREAMBLE, b"\0")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "wb") as f:
        f.write(pre)
        f.write(hbytes)
        for b in blobs:
            f.write(b)
    return {"path": str(out), "bytes": PREAMBLE + len(hbytes) + offset, "header": len(hbytes),
            "tiles": {s["uid"]: sum(len(lv["tiles"]) for lv in s["levels"]) for s in sets}}


# ---------------------------------------------------------------------------
# Lecture (tests, `check`, lecture croisée)
# ---------------------------------------------------------------------------

class Contacts:
    def __init__(self, path: str | Path):
        self.blob = Path(path).read_bytes()
        magic, version, hlen, hsha = struct.unpack_from("<8sII32s", self.blob)
        if magic != MAGIC:
            raise ValueError("ce n'est pas un fichier de contacts Geno")
        if version != VERSION:
            raise ValueError(f"contacts version {version}, ce lecteur lit la {VERSION}")
        hb = self.blob[PREAMBLE:PREAMBLE + hlen]
        if hashlib.sha256(hb).digest() != hsha:
            raise ValueError("en-tête : empreinte divergente")
        self.header = json.loads(_inflate(hb))
        self.base = PREAMBLE + hlen

    def tileset(self, uid: str) -> dict:
        return next(t for t in self.header["tilesets"] if t["uid"] == uid)

    def tile(self, uid: str, z: int, x: int, y: int) -> np.ndarray | None:
        spec = self.tileset(uid)["levels"][z]["tiles"].get(f"{z}.{x}.{y}")
        if spec is None:
            return None
        d = self.blob[self.base + spec["offset"]: self.base + spec["offset"] + spec["length"]]
        if digest(d) != spec["sha256"]:
            raise ValueError("tuile : empreinte divergente")
        return decode_tile(d)

    def pixel(self, uid: str, z: int, row: int, col: int) -> float:
        """La valeur d'un pixel, miroir compris (0 si la tuile n'existe pas)."""
        if col > row:
            row, col = col, row
        t = self.tile(uid, z, col // TILE, row // TILE)
        return 0.0 if t is None else float(t[row % TILE, col % TILE])


# ---------------------------------------------------------------------------
# Pistes d'igv.js
# ---------------------------------------------------------------------------

def _gz(path: Path, lines) -> dict:
    data = gzip.compress("".join(lines).encode(), mtime=0)
    path.write_bytes(data)
    return {"path": path.name, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def igv_tracks(out_dir: Path, *, genes, ctcf, loops, e1_rows, cytoband: Path, sizes: dict[str, int]) -> dict:
    """Les fichiers que lit igv.js, chacun petit assez pour être lu en entier (pas d'index)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    key = lambda c: (CHROMS.index(c) if c in CHROMS else 99, c)  # noqa: E731
    files = {}
    g = sorted(((c, s, e, d) for c, s, e, d in genes if c in sizes), key=lambda r: (key(r[0]), r[1]))
    files["genes"] = _gz(out_dir / "genes.bed.gz",
                         (f"{c}\t{s}\t{e}\t{d['name']}\t0\t{d['strand']}\n" for c, s, e, d in g))
    k = sorted(((c, s, e, d) for c, s, e, d in ctcf if c in sizes), key=lambda r: (key(r[0]), r[1]))
    files["ctcf"] = _gz(out_dir / "ctcf.bed.gz",
                        (f"{c}\t{s}\t{e}\tCTCF\t{min(1000, int(d['signal'] * 10))}\t{d['motif'] or '.'}\n"
                         for c, s, e, d in k))
    seen, lp = set(), []
    for c, s, e, d in loops:
        pc, ps, pe = d["partner"]
        a, b = (c, s, e), (pc, ps, pe)
        pair = (a, b) if (key(a[0]), a[1]) <= (key(b[0]), b[1]) else (b, a)
        if pair not in seen:
            seen.add(pair)
            lp.append((pair, d["observed"]))
    lp.sort(key=lambda r: (key(r[0][0][0]), r[0][0][1]))
    files["loops"] = _gz(out_dir / "loops.bedpe.gz",
                         (f"{a[0]}\t{a[1]}\t{a[2]}\t{b[0]}\t{b[1]}\t{b[2]}\tboucle\t{n}\n" for (a, b), n in lp))
    files["compartments"] = _gz(out_dir / "compartments.bedgraph.gz",
                                (f"{c}\t{s}\t{e}\t{v:.4f}\n" for c, s, e, v in e1_rows if np.isfinite(v)))
    data = cytoband.read_bytes()
    (out_dir / "cytoband.txt.gz").write_bytes(data)
    files["cytoband"] = {"path": "cytoband.txt.gz", "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    sz = "".join(f"{c}\t{sizes[c]}\n" for c in [*CHROMS, "chrY", "chrM"] if c in sizes).encode()
    (out_dir / "chrom.sizes").write_bytes(sz)
    files["chromsizes"] = {"path": "chrom.sizes", "bytes": len(sz), "sha256": hashlib.sha256(sz).hexdigest()}
    return files
