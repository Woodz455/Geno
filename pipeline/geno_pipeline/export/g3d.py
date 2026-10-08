"""Le format `.g3d` : un fichier qu'on lit par morceaux, dans l'ordre où on en a besoin.

Un `.g3d` porte une structure (le médoïde d'un ensemble) à plusieurs résolutions, et
il est disposé pour la requête HTTP Range. Trois décisions le façonnent, et chacune
vient d'une mesure faite avant lui.

**Tout préfixe est utile.** Le préambule dit où finit l'en-tête et où finissent les
données du premier rendu ; l'en-tête suit, puis — contiguës — les colonnes dont le
premier rendu a besoin, et seulement elles. Un client qui lit les 64 premiers Kio en
une requête a, si le fichier est bien construit, de quoi dessiner. C'est la seule
disposition compatible avec un premier rendu en 1,5 s sous un réseau où chaque requête
coûte une demi-seconde de latence.

**Des colonnes, pas des enregistrements.** La semaine 2 a mesuré que le coût d'une
requête suit le nombre de features décodées, dominé par le parsing JSON par
enregistrement. Ici chaque grandeur est un tableau typé, compressé à part, et un client
ne décompresse que les colonnes qu'il lit.

**Des positions quantifiées sur 16 bits.** La semaine 7 a mesuré que la position
radiale d'une bille varie de 276 nm d'un tirage à l'autre. Stocker cette position en
float32 — sept chiffres significatifs — serait afficher une précision que le modèle n'a
pas. Un pas de quantification de 0,15 nm sur un noyau de 10 µm reste 3 600 fois sous
l'incertitude mesurée, pour moitié moins d'octets.

La compression est **par bloc**, jamais par `Content-Encoding` HTTP : une requête Range
sur une réponse compressée à la volée désigne des octets du flux compressé, ce qui ne
sert à rien. PMTiles et COPC font le même choix pour la même raison.
"""

from __future__ import annotations

import hashlib
import json
import struct
import zlib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

MAGIC = b"\x89G3D\r\n\x1a\n"     # à la PNG : détecte un transfert en mode texte
VERSION = 1
PREAMBLE = 64
U16 = 65535

# Chaînes de filtres essayées pour chaque colonne ; la plus courte gagne, et le
# choix est écrit dans l'en-tête. Le lecteur n'a pas à deviner.
CHAINS = (
    ("deflate",),
    ("shuffle", "deflate"),
    ("delta", "deflate"),
    ("delta", "shuffle", "deflate"),
)


# --------------------------------------------------------------------------
# Filtres
# --------------------------------------------------------------------------


def _deflate(raw: bytes) -> bytes:
    c = zlib.compressobj(level=9, wbits=-15)          # deflate brut, sans en-tête zlib
    return c.compress(raw) + c.flush()


def _inflate(data: bytes) -> bytes:
    return zlib.decompress(data, wbits=-15)


def _delta(a: np.ndarray, planes: int) -> np.ndarray:
    """Différences successives **par plan**, en arithmétique modulaire.

    Modulaire, parce qu'un entier non signé qui « descend » doit rester un entier non
    signé : en u16, 3 − 5 vaut 65534, et la somme cumulée au décodage retombe sur 3.
    Aucune information n'est perdue, et aucun cas particulier n'est nécessaire.
    """
    p = a.reshape(planes, -1)
    out = p.copy()
    out[:, 1:] = p[:, 1:] - p[:, :-1]
    return out.reshape(-1)


def _undelta(a: np.ndarray, planes: int) -> np.ndarray:
    p = a.reshape(planes, -1)
    return np.cumsum(p, axis=1, dtype=a.dtype).reshape(-1)


def _shuffle(raw: bytes, itemsize: int) -> bytes:
    """Regroupe les octets de même rang : tous les octets de poids faible, puis les forts.

    Après un delta, les octets de poids fort sont presque tous nuls ; rangés ensemble,
    ils forment une longue plage que deflate réduit à presque rien.
    """
    if itemsize == 1:
        return raw
    return np.frombuffer(raw, np.uint8).reshape(-1, itemsize).T.tobytes()


def _unshuffle(raw: bytes, itemsize: int) -> bytes:
    if itemsize == 1:
        return raw
    return np.frombuffer(raw, np.uint8).reshape(itemsize, -1).T.tobytes()


def encode(arr: np.ndarray, planes: int, chain: tuple[str, ...]) -> bytes:
    a = np.ascontiguousarray(arr).reshape(-1)
    if "delta" in chain:
        a = _delta(a, planes)
    raw = a.tobytes()
    if "shuffle" in chain:
        raw = _shuffle(raw, a.dtype.itemsize)
    return _deflate(raw)


def decode(data: bytes, dtype: str, planes: int, chain: list[str]) -> np.ndarray:
    raw = _inflate(data)
    dt = np.dtype(dtype)
    if "shuffle" in chain:
        raw = _unshuffle(raw, dt.itemsize)
    a = np.frombuffer(raw, dt).copy()
    if "delta" in chain:
        a = _undelta(a, planes)
    return a


def best(arr: np.ndarray, planes: int) -> tuple[tuple[str, ...], bytes, dict]:
    """Essaie chaque chaîne de filtres, garde la plus courte, et rend les tailles de toutes."""
    sizes = {}
    winner = None
    for chain in CHAINS:
        blob = encode(arr, planes, chain)
        sizes["+".join(chain)] = len(blob)
        if winner is None or len(blob) < len(winner[1]):
            winner = (chain, blob)
    assert winner is not None
    return winner[0], winner[1], sizes


def digest(data: bytes) -> str:
    """sha256 tronqué à 16 octets, en hexadécimal.

    Un hash par colonne, pas par fichier : un client qui ne lit que des morceaux ne
    peut rien faire d'une empreinte du fichier entier. 128 bits suffisent largement à
    détecter une corruption ; ce n'est pas une signature.
    """
    return hashlib.sha256(data).hexdigest()[:32]


# --------------------------------------------------------------------------
# Niveaux
# --------------------------------------------------------------------------


@dataclass
class Level:
    """Une résolution du modèle : des billes, leurs positions, et d'où elles viennent."""

    name: str
    bp_per_bead: int
    positions_nm: np.ndarray          # (n, 3) dans le repère du niveau
    copy: np.ndarray                  # (n,) index dans `copies`
    start: np.ndarray                 # (n,) pb, 0-based
    end: np.ndarray                   # (n,) pb, exclu
    variability_nm: np.ndarray        # (n,) écart-type de profondeur sur l'ensemble
    copies: list[str]                 # étiquette par copie, ex. « chr7:a »
    radius_nm: list[float]            # rayon d'une bille, par copie
    evidence: str = "simulated"
    variability: dict = field(default_factory=lambda: {"kind": "radial_sd"})
    transform: list[float] = field(
        default_factory=lambda: [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    )
    fit: dict | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def n(self) -> int:
        return len(self.start)

    def check(self) -> None:
        n = self.n
        for name, a in (("positions_nm", self.positions_nm), ("copy", self.copy),
                        ("end", self.end), ("variability_nm", self.variability_nm)):
            if len(a) != n:
                raise ValueError(f"{self.name} : {name} a {len(a)} lignes pour {n} billes")
        if self.positions_nm.shape != (n, 3):
            raise ValueError(f"{self.name} : positions de forme {self.positions_nm.shape}")
        if not np.isfinite(self.positions_nm).all():
            raise ValueError(f"{self.name} : positions non finies")
        if (self.end <= self.start).any():
            raise ValueError(f"{self.name} : une bille a une fin avant son début")
        if self.copy.max() >= len(self.copies) or self.copy.min() < 0:
            raise ValueError(f"{self.name} : index de copie hors de la table")
        if len(self.copies) > 255:
            raise ValueError(f"{self.name} : plus de 255 copies ne tiennent pas en u8")
        if self.start.max() >= 2**32 or self.end.max() >= 2**32:
            raise ValueError(f"{self.name} : coordonnée au-delà de 2³² pb")


def octree(points: np.ndarray, max_beads: int) -> list[np.ndarray]:
    """Découpe adaptative : un nœud se coupe en huit tant qu'il porte trop de billes.

    Rend la liste des feuilles, chacune un tableau d'indices. Une feuille est un chunk :
    l'unité de lecture du client. Le découpage suit la densité, pas une grille fixe —
    un noyau a un cœur dense et une périphérie creuse, et une grille régulière
    fabriquerait des chunks vides à côté de chunks énormes.
    """
    leaves: list[np.ndarray] = []

    def split(idx: np.ndarray, depth: int) -> None:
        if len(idx) <= max_beads or depth >= 12:
            leaves.append(idx)
            return
        p = points[idx]
        mid = (p.min(axis=0) + p.max(axis=0)) / 2.0
        code = ((p[:, 0] > mid[0]).astype(int)
                | ((p[:, 1] > mid[1]).astype(int) << 1)
                | ((p[:, 2] > mid[2]).astype(int) << 2))
        for o in range(8):
            sub = idx[code == o]
            if len(sub):
                split(sub, depth + 1)

    split(np.arange(len(points)), 0)
    return leaves


# --------------------------------------------------------------------------
# Écriture
# --------------------------------------------------------------------------


@dataclass
class _Column:
    name: str
    dtype: str
    planes: int
    chain: tuple[str, ...]
    blob: bytes
    raw: int
    sizes: dict


def _column(name: str, arr: np.ndarray, planes: int = 1) -> _Column:
    chain, blob, sizes = best(arr, planes)
    return _Column(name, arr.dtype.str.lstrip("<|"), planes, chain, blob, arr.nbytes, sizes)


def _quantize(x: np.ndarray, lo: np.ndarray, step: np.ndarray) -> np.ndarray:
    q = np.rint((x - lo) / step)
    return np.clip(q, 0, U16).astype(np.uint16)


def precision(variability_nm: np.ndarray, fraction: float) -> float:
    """Erreur de quantification maximale admise : une fraction de l'incertitude du niveau.

    La règle est fixée sur l'incertitude, pas sur un budget d'octets. Mesuré sur le
    noyau de la semaine 7, l'écart-type de profondeur médian vaut 280 nm ; à 1 %, l'erreur
    maximale admise est 2,8 nm et les positions passent de 40,6 à 26,4 ko — 35 % de moins
    qu'en quantification pleine sur 16 bits, dont le pas de 0,14 nm ne décrivait rien que
    le modèle sache.
    """
    med = float(np.median(variability_nm))
    if not np.isfinite(med) or med <= 0:
        raise ValueError("incertitude nulle ou absente : impossible d'en déduire une précision")
    return fraction * med


def _prepare(level: Level, max_beads: int, fraction: float):
    """Ordonne les billes en chunks et construit toutes les colonnes d'un niveau."""
    level.check()
    pos = np.asarray(level.positions_nm, dtype=np.float64)
    var = np.asarray(level.variability_nm, dtype=np.float64)
    err = precision(var, fraction)

    lo = pos.min(axis=0)
    hi = pos.max(axis=0)
    # Le pas vaut deux fois l'erreur admise — sauf si la boîte est trop grande pour
    # 16 bits à ce pas, auquel cas c'est la boîte qui décide.
    step = np.maximum(np.maximum((hi - lo) / U16, 2.0 * err), 1e-9)
    var_step = max(2.0 * err, float(var.max()) / U16, 1e-9)

    chunks = []
    order = []
    for leaf in octree(pos, max_beads):
        leaf = leaf[np.lexsort((level.start[leaf], level.copy[leaf]))]
        first = len(order)
        order.extend(leaf.tolist())
        p = pos[leaf]
        centre = (p.min(axis=0) + p.max(axis=0)) / 2.0
        chunks.append(
            {
                "first": first,
                "n": len(leaf),
                "bbox": [*p.min(axis=0).round(2).tolist(), *p.max(axis=0).round(2).tolist()],
                "sphere": [*centre.round(2).tolist(),
                           round(float(np.linalg.norm(p - centre, axis=1).max()), 2)],
                "leaf": leaf,
            }
        )
    order = np.asarray(order, dtype=np.int64)

    for c in chunks:
        leaf = c.pop("leaf")
        q = _quantize(pos[leaf], lo, step)
        c["first_frame"] = [
            _column("copy", level.copy[leaf].astype(np.uint8)),
            # Trois plans x, y, z à la suite : le delta travaille le long de la chaîne
            # sur un axe à la fois, là où les valeurs voisines se ressemblent.
            _column("position", np.ascontiguousarray(q.T), planes=3),
        ]
        c["lazy"] = [
            _column("variability", _quantize(var[leaf], np.zeros(1), np.array([var_step]))),
            _column("start", level.start[leaf].astype(np.uint32)),
            _column("end", level.end[leaf].astype(np.uint32)),
        ]

    # Index génomique : les identifiants de billes dans l'ordre (copie, début). Les
    # billes d'une copie ne se chevauchent pas, donc une dichotomie sur les débuts
    # suffit — pas besoin du maximum de fin par bloc de la semaine 2, qui servait des
    # gènes qui, eux, se chevauchent.
    ident = np.empty(len(order), dtype=np.int64)
    ident[order] = np.arange(len(order))          # bille d'origine → identifiant stocké
    g = np.lexsort((level.start, level.copy))
    copy_sorted = level.copy[g]
    offsets = np.searchsorted(copy_sorted, np.arange(len(level.copies) + 1)).tolist()
    index = [
        _column("ids", ident[g].astype(np.uint32)),
        _column("starts", level.start[g].astype(np.uint32)),
        _column("ends", level.end[g].astype(np.uint32)),
    ]

    meta = {
        "name": level.name,
        "bp_per_bead": int(level.bp_per_bead),
        "n_beads": level.n,
        "evidence": level.evidence,
        "copies": list(level.copies),
        "radius_nm": [round(float(r), 3) for r in level.radius_nm],
        "frame": {"units": "nm", "transform": [float(v) for v in level.transform],
                  "fit": level.fit},
        # L'origine n'est pas arrondie : un arrondi au 1e-4 nm suffisait à faire
        # dépasser l'erreur annoncée de 3·10⁻⁵ nm, ce qu'un test exact a vu.
        "quant": {"origin": lo.tolist(), "step": step.tolist(),
                  "max_error_nm": float(step.max() / 2.0),
                  "rule": f"erreur max = {fraction:.0%} de l'incertitude médiane du niveau"},
        "variability": {**level.variability, "unit": "nm", "step": var_step},
        "index": {"copy_offsets": offsets},
        "notes": list(level.notes),
    }
    return meta, chunks, index, order


def write(
    path: str | Path,
    levels: list[Level],
    *,
    first: str | None = None,
    max_beads: int = 16_384,
    fraction: float = 0.01,
    header: dict | None = None,
) -> dict:
    """Écrit le fichier et rend un compte des octets, colonne par colonne.

    `first` désigne le niveau du premier rendu : ses colonnes `copy` et `position` sont
    placées juste après l'en-tête, contiguës, et le préambule dit où elles finissent.
    """
    if not levels:
        raise ValueError("un .g3d sans niveau ne sert à rien")
    names = [lv.name for lv in levels]
    if len(set(names)) != len(names):
        raise ValueError(f"noms de niveaux en double : {names}")
    first = first or levels[0].name
    if first not in names:
        raise ValueError(f"niveau du premier rendu inconnu : {first}")

    prepared = [_prepare(lv, max_beads, fraction) for lv in levels]

    # Ordre des blocs : d'abord le premier rendu, puis tout le reste niveau par niveau.
    blocks: list[tuple[dict, _Column]] = []
    out_levels = []
    for lv, (meta, chunks, index, _) in zip(levels, prepared):
        meta = dict(meta)
        meta["chunks"] = []
        for c in chunks:
            cm = {k: v for k, v in c.items() if k not in ("first_frame", "lazy")}
            cm["columns"] = {}
            meta["chunks"].append(cm)
        meta["index"]["columns"] = {}
        out_levels.append((lv, meta, chunks, index))

    def placed(lv_name, group):
        for lv, meta, chunks, index in out_levels:
            if (lv.name == first) != (lv_name == "first"):
                continue
            for cm, c in zip(meta["chunks"], chunks):
                for col in c[group]:
                    blocks.append((cm["columns"], col))
            if group == "lazy":
                for col in index:
                    blocks.append((meta["index"]["columns"], col))

    placed("first", "first_frame")
    n_first = len(blocks)
    placed("first", "lazy")
    placed("rest", "first_frame")
    placed("rest", "lazy")

    offset = 0
    report = {"columns": {}}
    for target, col in blocks:
        target[col.name] = {
            "offset": offset,
            "length": len(col.blob),
            "raw": col.raw,
            "dtype": col.dtype,
            "planes": col.planes,
            "filters": list(col.chain),
            "sha256": digest(col.blob),
        }
        offset += len(col.blob)
    first_frame_data = sum(len(c.blob) for _, c in blocks[:n_first])

    hdr = {
        "format": "geno-g3d",
        "version": VERSION,
        **(header or {}),
        "first": first,
        "levels": [meta for _, meta, _, _ in out_levels],
    }
    hbytes = _deflate(json.dumps(hdr, ensure_ascii=False, separators=(",", ":")).encode())
    data_start = PREAMBLE + len(hbytes)
    first_end = data_start + first_frame_data

    pre = struct.pack(
        "<8sIIII32s",
        MAGIC, VERSION, len(hbytes), first_end, 0, hashlib.sha256(hbytes).digest(),
    )
    pre = pre.ljust(PREAMBLE, b"\0")

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(pre)
        fh.write(hbytes)
        for _, col in blocks:
            fh.write(col.blob)

    for lv, meta, chunks, index in out_levels:
        cols = {}
        for c in chunks:
            for col in [*c["first_frame"], *c["lazy"]]:
                agg = cols.setdefault(col.name, {"raw": 0, "stored": 0, "chains": {}})
                agg["raw"] += col.raw
                agg["stored"] += len(col.blob)
                for k, v in col.sizes.items():
                    agg["chains"][k] = agg["chains"].get(k, 0) + v
        for col in index:
            cols["index." + col.name] = {"raw": col.raw, "stored": len(col.blob),
                                         "chains": dict(col.sizes)}
        report["columns"][lv.name] = cols

    report.update(
        {
            "path": str(path),
            "bytes": data_start + offset,
            "header": len(hbytes),
            "first_frame_end": first_end,
            "chunks": {lv.name: len(meta["chunks"]) for lv, meta, _, _ in out_levels},
        }
    )
    return report


# --------------------------------------------------------------------------
# Lecture (pour les tests, `geno g3d inspect`, et pour garder le format honnête)
# --------------------------------------------------------------------------


class CorruptError(ValueError):
    """Une empreinte ne correspond pas : le bloc lu n'est pas celui qui a été écrit."""


def read_header(blob: bytes) -> tuple[dict, int, int]:
    """Lit préambule et en-tête depuis un préfixe. Rend `(en-tête, début des données, fin du premier rendu)`."""
    if len(blob) < PREAMBLE:
        raise ValueError("préfixe plus court que le préambule")
    magic, version, hlen, first_end, _, hsha = struct.unpack_from("<8sIIII32s", blob)
    if magic != MAGIC:
        raise ValueError("ce n'est pas un .g3d (signature)")
    if version != VERSION:
        raise ValueError(f".g3d version {version}, ce lecteur lit la {VERSION}")
    hbytes = blob[PREAMBLE:PREAMBLE + hlen]
    if len(hbytes) < hlen:
        raise ValueError("préfixe trop court pour l'en-tête")
    if hashlib.sha256(hbytes).digest() != hsha:
        raise CorruptError("en-tête : empreinte divergente")
    return json.loads(_inflate(hbytes)), PREAMBLE + hlen, first_end


@dataclass
class Read:
    header: dict
    levels: dict[str, dict]


def read(path: str | Path) -> Read:
    """Relit tout et vérifie chaque empreinte. Rend les colonnes décodées par niveau."""
    blob = Path(path).read_bytes()
    hdr, base, _ = read_header(blob)

    def col(spec: dict) -> np.ndarray:
        data = blob[base + spec["offset"]: base + spec["offset"] + spec["length"]]
        if digest(data) != spec["sha256"]:
            raise CorruptError("colonne : empreinte divergente")
        return decode(data, spec["dtype"], spec["planes"], spec["filters"])

    out = {}
    for lv in hdr["levels"]:
        q = lv["quant"]
        parts = {k: [] for k in ("copy", "position", "variability", "start", "end")}
        for c in lv["chunks"]:
            for name in parts:
                parts[name].append(col(c["columns"][name]))
        pos_q = np.concatenate([p.reshape(3, -1).T for p in parts["position"]])
        positions = np.asarray(q["origin"]) + pos_q * np.asarray(q["step"])
        out[lv["name"]] = {
            "meta": lv,
            "copy": np.concatenate(parts["copy"]),
            "positions_nm": positions,
            "variability_nm": np.concatenate(parts["variability"]) * lv["variability"]["step"],
            "start": np.concatenate(parts["start"]),
            "end": np.concatenate(parts["end"]),
            "index": {k: col(v) for k, v in lv["index"]["columns"].items()},
        }
    return Read(hdr, out)
