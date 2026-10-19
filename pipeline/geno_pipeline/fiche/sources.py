"""Les sources réelles de la fiche, lues une fois et rangées en intervalles par chromosome.

Chaque lecteur part du fichier verrouillé du manifeste (`data/raw/<tier>/<id>`) et rend des
tableaux triés : la composition des fiches (`build.py`) ne fait plus ensuite que des
recouvrements. Les lectures longues (GENCODE, ClinVar, GWAS) sont mises en cache sous
`data/cache/fiche/`, indexé par l'empreinte du manifeste : un fichier changé invalide son
cache, un cache n'invente jamais un fichier.

Trois sources publiées en **hg19** (Rao 2014 et le Repli-seq UW) sont levées en hg38 par la
chaîne Ensembl. Ce qui ne se lève pas est perdu, et compté.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import math
import pickle
import re
import sys
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..liftover import Chain, ucsc_name
from ..parsers import parse_cytoband
from ..sequence import TwoBit

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RAW = DATA / "raw"
CACHE = DATA / "cache" / "fiche"

MAIN = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])


# ---------------------------------------------------------------------------
# Manifeste
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Entry:
    id: str
    tier: str
    status: str
    license: str
    url: str
    sha256: str
    notes: str

    @property
    def path(self) -> Path:
        return RAW / self.tier / self.id


def manifest(path: Path | None = None) -> dict[str, Entry]:
    out = {}
    for line in (path or DATA / "manifest.tsv").read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("id\t"):
            continue
        c = line.split("\t")
        out[c[0]] = Entry(c[0], c[1], c[2], c[3], c[4], c[5], c[7] if len(c) > 7 else "")
    return out


def locked(m: dict[str, Entry], key: str) -> Entry:
    """Une source utilisable : verrouillée, et présente avec la bonne taille."""
    e = m[key]
    if e.status != "ok":
        raise RuntimeError(f"{key} : statut {e.status} — `make data-lock` d'abord")
    if not e.path.exists():
        raise FileNotFoundError(f"{key} : {e.path} absent — `make data-…`")
    return e


def _text(path: Path):
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if gz else open(path, encoding="utf-8", errors="replace")


def _cached(name: str, entries: list[Entry], build, extra: str = ""):
    """`build()` une fois par jeu d'empreintes (et de paramètres) ; ensuite, le pickle."""
    key = "-".join(e.sha256[:12] for e in entries)
    if extra:
        key += "-" + hashlib.sha256(extra.encode()).hexdigest()[:8]
    p = CACHE / f"{name}.{key}.pkl"
    if p.exists():
        with open(p, "rb") as f:
            return pickle.load(f)
    obj = build()
    CACHE.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".part")
    with open(tmp, "wb") as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    tmp.replace(p)
    return obj


# ---------------------------------------------------------------------------
# Intervalles
# ---------------------------------------------------------------------------

@dataclass
class Intervals:
    """Des intervalles par chromosome, triés par début, chacun avec une charge (un index
    dans `items`). La requête de recouvrement borne sa fenêtre par la plus grande longueur
    du chromosome — assez pour des pistes dont aucun élément n'est un gène de 2 Mb."""

    starts: dict[str, np.ndarray] = field(default_factory=dict)
    ends: dict[str, np.ndarray] = field(default_factory=dict)
    idx: dict[str, np.ndarray] = field(default_factory=dict)
    maxlen: dict[str, int] = field(default_factory=dict)
    items: list = field(default_factory=list)

    @classmethod
    def build(cls, rows: list[tuple[str, int, int, object]]) -> Intervals:
        out = cls()
        by: dict[str, list[tuple[int, int, int]]] = {}
        for chrom, s, e, item in rows:
            by.setdefault(chrom, []).append((s, e, len(out.items)))
            out.items.append(item)
        for chrom, r in by.items():
            a = np.array(r, dtype=np.int64)
            a = a[np.lexsort((a[:, 1], a[:, 0]))]
            out.starts[chrom], out.ends[chrom], out.idx[chrom] = a[:, 0], a[:, 1], a[:, 2]
            out.maxlen[chrom] = int((a[:, 1] - a[:, 0]).max())
        return out

    def query(self, chrom: str, s: int, e: int) -> np.ndarray:
        """Positions (dans l'ordre trié du chromosome) des intervalles qui touchent [s, e)."""
        if chrom not in self.starts:
            return np.empty(0, dtype=np.int64)
        st, en = self.starts[chrom], self.ends[chrom]
        lo = int(np.searchsorted(st, s - self.maxlen[chrom], side="left"))
        hi = int(np.searchsorted(st, e, side="left"))
        r = np.arange(lo, hi)
        return r[en[lo:hi] > s]

    def get(self, chrom: str, s: int, e: int) -> list:
        k = self.query(chrom, s, e)
        return [self.items[i] for i in self.idx[chrom][k]] if k.size else []

    def overlap_bp(self, chrom: str, s: int, e: int) -> list[tuple[object, int]]:
        k = self.query(chrom, s, e)
        if not k.size:
            return []
        st, en, ix = self.starts[chrom][k], self.ends[chrom][k], self.idx[chrom][k]
        bp = np.minimum(en, e) - np.maximum(st, s)
        return [(self.items[i], int(b)) for i, b in zip(ix, bp)]


@dataclass
class Points:
    """Des positions ponctuelles par chromosome, triées, avec une charge par point."""

    pos: dict[str, np.ndarray] = field(default_factory=dict)
    idx: dict[str, np.ndarray] = field(default_factory=dict)
    items: list = field(default_factory=list)

    @classmethod
    def build(cls, rows: list[tuple[str, int, object]]) -> Points:
        out = cls()
        by: dict[str, list[tuple[int, int]]] = {}
        for chrom, p, item in rows:
            by.setdefault(chrom, []).append((p, len(out.items)))
            out.items.append(item)
        for chrom, r in by.items():
            a = np.array(r, dtype=np.int64)
            a = a[np.argsort(a[:, 0], kind="stable")]
            out.pos[chrom], out.idx[chrom] = a[:, 0], a[:, 1]
        return out

    def range(self, chrom: str, s: int, e: int) -> np.ndarray:
        """Indices dans `items` des points de [s, e)."""
        if chrom not in self.pos:
            return np.empty(0, dtype=np.int64)
        p = self.pos[chrom]
        return self.idx[chrom][np.searchsorted(p, s): np.searchsorted(p, e)]


# ---------------------------------------------------------------------------
# Coordonnées
# ---------------------------------------------------------------------------

def cytobands(m) -> Intervals:
    e = locked(m, "hg38.cytoband")
    return Intervals.build([(c, s, t, r["name"]) for c, s, t, r in parse_cytoband(e.path) if r["band"]])


def chain_hs1(m) -> Chain:
    return Chain(locked(m, "chain.hg38.to.hs1").path)


def chain_hg19(m) -> Chain:
    return Chain(locked(m, "chain.hg19.to.hg38").path)


# ---------------------------------------------------------------------------
# Gènes : GENCODE, GTEx, gnomAD
# ---------------------------------------------------------------------------

_ATTR = re.compile(r'(\S+) "([^"]*)"')


def genes(m) -> tuple[Intervals, dict]:
    """Gènes GENCODE (basic) : id sans version, nom, biotype, brin, et s'il porte un
    transcrit MANE Select. Rend aussi la release, lue dans l'URL verrouillée."""
    e = locked(m, "gencode.annotation.gtf")

    def build():
        rows, mane = [], set()
        with _text(e.path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.split("\t", 8)
                if c[2] == "gene":
                    a = dict(_ATTR.findall(c[8]))
                    rows.append((c[0], int(c[3]) - 1, int(c[4]), {
                        "id": a["gene_id"].split(".")[0], "name": a.get("gene_name", a["gene_id"]),
                        "type": a.get("gene_type", ""), "strand": c[6],
                        "start": int(c[3]) - 1, "end": int(c[4])}))
                elif c[2] == "transcript" and 'tag "MANE_Select"' in c[8]:
                    a = dict(_ATTR.findall(c[8]))
                    mane.add(a["gene_id"].split(".")[0])
        for _, _, _, g in rows:
            g["mane"] = g["id"] in mane
        return rows

    rows = _cached("genes", [e], build)
    release = re.search(r"release_(\d+)", e.url)
    return Intervals.build(rows), {"release": f"GENCODE v{release.group(1)}" if release else "GENCODE"}


def gtex(m) -> tuple[dict[str, tuple[float, int, float]], list[str]]:
    """TPM médian par gène : en lymphocytes EBV (≈ GM12878), et le tissu où il est maximal."""
    e = locked(m, "gtex.median.tpm")

    def build():
        with _text(e.path) as f:
            f.readline()
            f.readline()
            head = f.readline().rstrip("\n").split("\t")
            tissues = head[2:]
            lcl = tissues.index("Cells_EBV-transformed_lymphocytes")
            out = {}
            for line in f:
                c = line.rstrip("\n").split("\t")
                v = np.array(c[2:], dtype=np.float64)
                k = int(np.argmax(v))
                out[c[0].split(".")[0]] = (float(v[lcl]), k, float(v[k]))
        return out, tissues

    return _cached("gtex", [e], build)


def gnomad(m) -> dict[str, tuple[float, float]]:
    """pLI et LOEUF (borne haute de l'IC à 90 % de o/e LoF), sur le transcrit MANE Select
    du gène quand il existe, sinon le canonique."""
    e = locked(m, "gnomad.constraint")

    def build():
        out: dict[str, tuple[float, float, int]] = {}
        with _text(e.path) as f:
            head = f.readline().rstrip("\n").split("\t")
            col = {k: i for i, k in enumerate(head)}
            for line in f:
                c = line.rstrip("\n").split("\t")
                gid = c[col["gene_id"]]
                if not gid.startswith("ENSG"):
                    continue
                rank = 2 if c[col["mane_select"]] == "true" else 1 if c[col["canonical"]] == "true" else 0
                if rank == 0:
                    continue

                def num(k):
                    v = c[col[k]]
                    return float(v) if v not in ("NA", "") else math.nan

                if gid not in out or rank > out[gid][2]:
                    out[gid] = (num("lof.pLI"), num("lof.oe_ci.upper"), rank)
        return {k: (a, b) for k, (a, b, _) in out.items()}

    return _cached("gnomad", [e], build)


# ---------------------------------------------------------------------------
# Régulation : cCRE, ChromHMM, CTCF orienté
# ---------------------------------------------------------------------------

def ccres(m) -> Intervals:
    """cCRE SCREEN de GM12878. `Low-DNase` : élément du registre **inactif** dans cette
    lignée — compté à part, il n'a pas de classe ici."""
    e = locked(m, "encode.ccres")
    rows = []
    with _text(e.path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            rows.append((c[0], int(c[1]), int(c[2]), (c[3], c[9])))
    return Intervals.build(rows)


def chromhmm(m) -> Intervals:
    e = locked(m, "chromhmm.gm12878")
    rows = []
    with _text(e.path) as f:
        for line in f:
            c = line.split("\t", 5)
            rows.append((c[0], int(c[1]), int(c[2]), c[3]))
    return Intervals.build(rows)


def _pwm(m) -> np.ndarray:
    """Matrice log-odds (19 × 4, ordre A C G T) de MA0139.1, fond uniforme, pseudo-compte 0,8."""
    pfm = json.loads(locked(m, "ctcf.motif.orientation").path.read_text())["pfm"]
    counts = np.array([pfm[b] for b in "ACGT"], dtype=np.float64).T + 0.8
    prob = counts / counts.sum(axis=1, keepdims=True)
    return np.log2(prob / 0.25)


_CODE = np.full(256, 4, dtype=np.int8)
for _i, _b in enumerate("ACGT"):
    _CODE[ord(_b)] = _i
    _CODE[ord(_b.lower())] = _i


def motif_scan(seq: str, pwm: np.ndarray) -> tuple[float, int, str]:
    """Meilleur score sur les deux brins : (score relatif 0–1, décalage, brin)."""
    w = len(pwm)
    if len(seq) < w:
        return 0.0, -1, "."
    code = _CODE[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
    lo, hi = float(pwm.min(axis=1).sum()), float(pwm.max(axis=1).sum())
    idx = np.lib.stride_tricks.sliding_window_view(code, w)
    cols = np.arange(w)
    plus = np.where(idx < 4, pwm[cols, np.minimum(idx, 3)], -10.0).sum(axis=1)
    rc = np.where(code < 4, 3 - code, 4)[::-1]
    idr = np.lib.stride_tricks.sliding_window_view(rc, w)
    minus = np.where(idr < 4, pwm[cols, np.minimum(idr, 3)], -10.0).sum(axis=1)
    kp, km = int(np.argmax(plus)), int(np.argmax(minus))
    if plus[kp] >= minus[km]:
        return (float(plus[kp]) - lo) / (hi - lo), kp, "+"
    return (float(minus[km]) - lo) / (hi - lo), len(seq) - w - km, "-"


MOTIF_MIN = 0.80     # seuil relatif usuel de JASPAR


def ctcf(m, genome: TwoBit) -> Intervals:
    """Pics CTCF IDR de GM12878, chacun avec l'orientation de son meilleur motif MA0139.1
    (scan des deux brins sur la séquence hg38 du pic). Sous le seuil, pas d'orientation —
    un pic sans motif reconnu existe, et c'est dit."""
    e = locked(m, "ctcf.peaks.gm12878")
    pwm = _pwm(m)

    def build():
        rows = []
        with _text(e.path) as f:
            for line in f:
                c = line.rstrip("\n").split("\t")
                chrom, s, t = c[0], int(c[1]), int(c[2])
                try:
                    seq = genome.fetch(chrom, s, t)
                except (KeyError, ValueError):
                    continue
                score, off, strand = motif_scan(seq, pwm)
                ok = score >= MOTIF_MIN
                rows.append((chrom, s, t, {"signal": round(float(c[6]), 2),
                                           "motif": strand if ok else None,
                                           "motif_at": s + off if ok else None,
                                           "motif_score": round(score, 3)}))
        return rows

    return Intervals.build(_cached("ctcf", [e, locked(m, "ctcf.motif.orientation")], build))


# ---------------------------------------------------------------------------
# Structure mesurée (Rao 2014, hg19 → hg38)
# ---------------------------------------------------------------------------

def _hg19_rows(path: Path):
    with _text(path) as f:
        for line in f:
            if line.startswith(("chr1\tx1", "#")):
                continue
            yield line.rstrip("\n").split("\t")


def subcompartments(m, c19: Chain) -> tuple[Intervals, dict]:
    e = locked(m, "rao2014.gm12878.subcompartments")
    rows, lost = [], 0
    for c in _hg19_rows(e.path):
        if c[3] in ("NA", ""):
            continue
        iv = c19.interval(c[0], int(c[1]), int(c[2]))
        if iv is None or iv.split:
            lost += 1
            continue
        rows.append((iv.chrom, iv.start, iv.end, c[3]))
    return Intervals.build(rows), {"lifted": len(rows), "lost": lost}


def domains(m, c19: Chain) -> tuple[Intervals, dict]:
    e = locked(m, "rao2014.gm12878.domains")
    rows, lost = [], 0
    for c in _hg19_rows(e.path):
        chrom = ucsc_name(c[0])
        iv = c19.interval(chrom, int(c[1]), int(c[2]))
        if iv is None or iv.split:
            lost += 1
            continue
        rows.append((iv.chrom, iv.start, iv.end, (iv.start, iv.end)))
    return Intervals.build(rows), {"lifted": len(rows), "lost": lost}


def loops(m, c19: Chain) -> tuple[Intervals, dict]:
    """Boucles HiCCUPS : chaque ancre est un intervalle, qui porte son partenaire et
    l'orientation du motif CTCF que Rao 2014 y a trouvé (`p` / `n` → `+` / `−`)."""
    e = locked(m, "rao2014.gm12878.loops")
    rows, lost, n = [], 0, 0
    o = {"p": "+", "n": "-"}
    for c in _hg19_rows(e.path):
        if len(c) < 30:
            continue
        n += 1
        a = c19.interval(ucsc_name(c[0]), int(c[1]), int(c[2]))
        b = c19.interval(ucsc_name(c[3]), int(c[4]), int(c[5]))
        if a is None or b is None or a.split or b.split:
            lost += 1
            continue
        oa, ob = o.get(c[23]), o.get(c[28])
        rows.append((a.chrom, a.start, a.end, {"anchor": [a.start, a.end], "motif": oa,
                                               "partner": [b.chrom, b.start, b.end], "partner_motif": ob,
                                               "observed": int(float(c[7]))}))
        rows.append((b.chrom, b.start, b.end, {"anchor": [b.start, b.end], "motif": ob,
                                               "partner": [a.chrom, a.start, a.end], "partner_motif": oa,
                                               "observed": int(float(c[7]))}))
    return Intervals.build(rows), {"loops": n, "lost": lost}


# ---------------------------------------------------------------------------
# Repères nucléaires
# ---------------------------------------------------------------------------

# Hansen et al. 2010 (PNAS), le « WaveSignal » du Repli-seq UW : chaque phase pèse selon
# son rang dans le cycle, G1b le plus tôt. Haut = réplication précoce.
REPLI_PHASES = (("repliseq.gm12878", 0.917), ("repliseq.gm12878.S1", 0.750),
                ("repliseq.gm12878.S2", 0.583), ("repliseq.gm12878.S3", 0.417),
                ("repliseq.gm12878.S4", 0.250), ("repliseq.gm12878.G2", 0.0))
REPLI_BIN = 10_000


def replication(m, c19: Chain) -> tuple[dict[str, tuple[np.ndarray, np.ndarray]], dict]:
    """Timing de réplication de GM12878 par bins de 10 kb hg19, levés en hg38 par leur
    centre. Rend, par chromosome hg38 : positions triées et valeurs."""
    import pyBigWig

    entries = [locked(m, k) for k, _ in REPLI_PHASES]

    def build():
        bws = [pyBigWig.open(str(e.path)) for e in entries]
        out: dict[str, list] = {}
        lost = kept = 0
        for chrom, size in bws[0].chroms().items():
            if chrom not in MAIN:
                continue
            nb = size // REPLI_BIN
            if nb == 0:
                continue
            vals = np.zeros(nb)
            ok = np.ones(nb, dtype=bool)
            for bw, (_, w) in zip(bws, REPLI_PHASES):
                v = np.array(bw.stats(chrom, 0, nb * REPLI_BIN, type="mean", nBins=nb), dtype=float)
                ok &= np.isfinite(v)
                vals += w * np.nan_to_num(v)
            centres = np.arange(nb) * REPLI_BIN + REPLI_BIN // 2
            for p, v in zip(centres[ok], vals[ok]):
                hit = c19.point(chrom, int(p))
                if hit is None:
                    lost += 1
                    continue
                kept += 1
                out.setdefault(hit.chrom, []).append((hit.pos, v))
        res = {}
        for c, r in out.items():
            a = np.array(r)
            a = a[np.argsort(a[:, 0])]
            res[c] = (a[:, 0].astype(np.int64), a[:, 1])
        return res, {"bins": kept, "lost": lost}

    return _cached("repli", entries + [locked(m, "chain.hg19.to.hg38")], build)


class BigWig:
    """Moyenne d'un bigWig sur un intervalle, et fraction des bins au-dessus d'un seuil."""

    def __init__(self, path: Path):
        import pyBigWig
        self.bw = pyBigWig.open(str(path))
        self.chroms = self.bw.chroms()

    def summary(self, chrom: str, s: int, e: int, above: float = 0.0) -> tuple[float, float]:
        if chrom not in self.chroms:
            return math.nan, math.nan
        e = min(e, self.chroms[chrom])
        iv = self.bw.intervals(chrom, s, e) if e > s else None
        if not iv:
            return math.nan, math.nan
        a = np.array(iv, dtype=np.float64)
        w = np.minimum(a[:, 1], e) - np.maximum(a[:, 0], s)
        v = a[:, 2]
        ok = np.isfinite(v) & (w > 0)
        if not ok.any():
            return math.nan, math.nan
        return float(np.average(v[ok], weights=w[ok])), float(np.average(v[ok] > above, weights=w[ok]))


def bedgraph(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    out: dict[str, list] = {}
    with _text(path) as f:
        for line in f:
            c = line.split()
            if len(c) < 4 or c[0] == "track":
                continue
            out.setdefault(c[0], []).append((int(c[1]), int(c[2]), float(c[3])))
    res = {}
    for k, r in out.items():
        a = np.array(r)
        a = a[np.argsort(a[:, 0])]
        res[k] = (a[:, 0].astype(np.int64), a[:, 1].astype(np.int64), a[:, 2])
    return res


def bedgraph_sum(track, chrom: str, s: int, e: int) -> float:
    """Somme de valeur × pb sur [s, e) — pour un rapport de deux pistes RPKM."""
    if chrom not in track:
        return math.nan
    st, en, v = track[chrom]
    lo, hi = np.searchsorted(en, s, side="right"), np.searchsorted(st, e)
    if hi <= lo:
        return math.nan
    w = np.minimum(en[lo:hi], e) - np.maximum(st[lo:hi], s)
    return float(np.sum(v[lo:hi] * w))


# ---------------------------------------------------------------------------
# Clinique : ClinVar, GWAS, ontologies
# ---------------------------------------------------------------------------

def obo_names(path: Path, prefix: str) -> dict[str, str]:
    names, cur = {}, None
    with _text(path) as f:
        for line in f:
            if line.startswith("[Term]"):
                cur = None
            elif line.startswith("id: ") and line[4:].startswith(prefix):
                cur = line[4:].strip()
            elif line.startswith("name: ") and cur:
                names[cur] = line[6:].strip()
                cur = None
    return names


def ontologies(m) -> dict[str, str]:
    em, eh = locked(m, "mondo.obo"), locked(m, "hpo.obo")
    return _cached("obo", [em, eh], lambda: {**obo_names(em.path, "MONDO:"), **obo_names(eh.path, "HP:")})


SIG = {"plp": "pathogène / probablement pathogène", "vus": "signification incertaine",
       "blb": "bénin / probablement bénin", "conf": "interprétations contradictoires",
       "other": "autre (facteur de risque, réponse médicamenteuse, non fourni…)"}

STARS = {"practice_guideline": 4, "reviewed_by_expert_panel": 3,
         "criteria_provided,_multiple_submitters,_no_conflicts": 2,
         "criteria_provided,_conflicting_classifications": 1,
         "criteria_provided,_conflicting_interpretations": 1,
         "criteria_provided,_single_submitter": 1}


def _sig(s: str) -> str:
    s = s.lower()
    if s.startswith("conflicting"):
        return "conf"
    if "pathogenic" in s and "conflicting" not in s and "benign" not in s:
        return "plp"
    if "uncertain" in s:
        return "vus"
    if "benign" in s and "pathogenic" not in s:
        return "blb"
    return "other"


def clinvar(m, keep: list[tuple[str, int, int]]) -> tuple[Points, dict]:
    """ClinVar GRCh38 : toutes les positions avec leur classe (pour les comptes), et le
    détail des seuls variants dont la fiche affiche le détail — les P/LP partout, et tous
    ceux des régions `keep` (les régions fines, où l'on clique jusqu'à la base).

    Maladies : l'identifiant MONDO si ClinVar en donne un, nommé par mondo.obo ; sinon HPO,
    nommé par hp.obo ; sinon le libellé ClinVar. Les numéros OMIM sont gardés comme
    identifiants — aucun texte OMIM, qui n'est pas redistribuable.
    """
    e = locked(m, "clinvar.vcf")

    def build():
        rows = []
        with _text(e.path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.rstrip("\n").split("\t")
                chrom = ucsc_name(c[0])
                pos = int(c[1]) - 1
                info = dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
                sig = _sig(info.get("CLNSIG", ""))
                detail = sig == "plp" or any(k == chrom and a <= pos < b for k, a, b in keep)
                if not detail:
                    rows.append((chrom, pos, sys.intern(sig)))      # un compte, rien d'autre
                    continue
                item = {"sig": sig}
                dn = info.get("CLNDN", "").split("|")
                db = info.get("CLNDISDB", "").split("|")
                dis = []
                for k, name in enumerate(dn):
                    ids = db[k].split(",") if k < len(db) else []
                    ids = [i.replace("MONDO:MONDO:", "MONDO:").replace("Human_Phenotype_Ontology:", "")
                           for i in ids if i.startswith(("MONDO:", "Human_Phenotype_Ontology:", "OMIM:", "Orphanet:"))]
                    name = name.replace("_", " ")
                    if name.lower() in ("not provided", "not specified"):
                        continue
                    dis.append([ids, name])
                item.update({
                    "id": c[2], "rs": "rs" + info["RS"] if "RS" in info else None,
                    "ref": c[3][:12] + ("…" if len(c[3]) > 12 else ""),
                    "alt": c[4][:12] + ("…" if len(c[4]) > 12 else ""),
                    "clnsig": info.get("CLNSIG", "").replace("_", " "),
                    "stars": STARS.get(info.get("CLNREVSTAT", ""), 0),
                    "genes": [g.split(":")[0] for g in info.get("GENEINFO", "").split("|") if g],
                    "hgvs": info.get("CLNHGVS"),
                    "kind": info.get("CLNVC", "").replace("_", " "),
                    "diseases": dis,
                })
                rows.append((chrom, pos, item))
        return rows

    rows = _cached("clinvar", [e], build, extra=repr(sorted(keep)))
    names = ontologies(m)
    for _, _, it in rows:
        if isinstance(it, str):
            continue
        for d in it.get("diseases", []):
            ids = d[0]
            for pre in ("MONDO:", "HP:"):
                hit = next((i for i in ids if i.startswith(pre) and i in names), None)
                if hit:
                    d.append(names[hit])
                    break
    date = re.search(r"clinvar_(\d{8})", e.notes)
    return Points.build(rows), {"release": date.group(1) if date else ""}


def gwas(m) -> Points:
    """GWAS Catalog : une ligne par association (SNP × étude). Coordonnées GRCh38 du
    catalogue ; les lignes sans position (haplotypes, SNP fusionnés) ne sont pas placées."""
    e = locked(m, "gwas.catalog")

    def build():
        rows = []
        with zipfile.ZipFile(e.path) as z:
            name = next(n for n in z.namelist() if n.endswith(".tsv"))
            with z.open(name) as raw:
                f = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
                head = f.readline().rstrip("\n").split("\t")
                col = {k: i for i, k in enumerate(head)}
                for line in f:
                    c = line.rstrip("\n").split("\t")
                    chrom, pos = c[col["CHR_ID"]], c[col["CHR_POS"]]
                    if not chrom or not pos.isdigit():
                        continue
                    try:
                        p = float(c[col["P-VALUE"]])
                    except ValueError:
                        continue
                    rows.append((ucsc_name(chrom), int(pos) - 1, {
                        "rs": c[col["SNPS"]], "trait": c[col["DISEASE/TRAIT"]],
                        "mapped": c[col["MAPPED_TRAIT"]], "uri": c[col["MAPPED_TRAIT_URI"]].split(",")[0].strip(),
                        "p": p, "pmid": c[col["PUBMEDID"]], "study": c[col["STUDY ACCESSION"]],
                        "genes": c[col["MAPPED_GENE"]]}))
        return rows

    return Points.build(_cached("gwas", [e], build))
