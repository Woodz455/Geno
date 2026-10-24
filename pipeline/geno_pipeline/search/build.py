"""L'index de recherche, `gm12878.search` (format `geno-search/1`).

Deux tables triées par clé, découpées en blocs, chaque bloc compressé et signé — un bloc par
requête Range, comme les tuiles de `.contacts` :

- **genes** : clé = nom en majuscules **et** identifiant ENSG (deux clés par gène), GENCODE
  (toutes les lignes `gene`, les mêmes que la fiche). Un nom porté par plusieurs loci (régions
  pseudo-autosomiques, paralogues homonymes) rend toutes ses lignes ;
- **rs** : clé = numéro du rsID, ClinVar (champ `RS`, toutes classes confondues) et catalogue
  GWAS (jetons `rs…` de la colonne SNPS, position GRCh38 du catalogue). Une ligne par position
  distincte, avec ses sources (bit 1 : ClinVar, bit 2 : GWAS).

Ce n'est **pas** dbSNP : un rsID sans variant ClinVar ni association GWAS n'y est pas, et la
page le dit (décision de la semaine 16).

Disposition :

    préambule (64 o) : signature \\x89GSR\\r\\n\\x1a\\n, version, longueur et sha256 de l'en-tête
    en-tête JSON, deflate : chromosomes, chromosomes du modèle, répertoire des blocs
    blocs deflate :
      genes  JSON [[clé, nom, ENSG, chrom, début, fin, brin, biotype], …]
      rs     colonnes petit-boutistes : Δrs u32[n] (depuis le premier du bloc), chrom u8[n],
             pos u32[n] (0-based), sources u8[n]
"""

from __future__ import annotations

import hashlib
import json
import re
import struct
from pathlib import Path

import numpy as np

from ..export.g3d import _deflate, _inflate, digest
from ..liftover import ucsc_name

MAGIC = b"\x89GSR\r\n\x1a\n"
VERSION = 1
PREAMBLE = 64
GENE_BLOCK = 1024
RS_BLOCK = 4096
CLINVAR, GWAS = 1, 2
MODEL = tuple(f"chr{c}" for c in [*range(1, 23), "X"])       # 46,XX
_RS = re.compile(r"\brs(\d+)\b")


def gene_rows(m) -> list[list]:
    """[clé, nom, ENSG, chrom, début, fin, brin, biotype], deux clés par gène, triées."""
    from ..fiche import sources as S

    iv, _ = S.genes(m)
    rows = []
    for chrom in iv.starts:
        for i in iv.idx[chrom]:
            g = iv.items[int(i)]
            name = g["name"]
            if not name.isascii():
                raise ValueError(f"nom de gène non ASCII : {name!r} — la comparaison TS ne serait plus la même")
            base = [g["name"], g["id"], chrom, int(g["start"]), int(g["end"]), g["strand"], g["type"]]
            rows.append([name.upper(), *base])
            rows.append([g["id"].upper(), *base])
    rows.sort(key=lambda r: (r[0], r[3], r[4], r[2]))
    return rows


def clinvar_rs(m) -> list[tuple[int, str, int]]:
    """(rs, chrom, pos 0-based) de chaque ligne du VCF ClinVar qui porte un `RS`."""
    from ..fiche import sources as S

    e = S.locked(m, "clinvar.vcf")

    def build():
        out = []
        with S._text(e.path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.split("\t", 8)
                k = c[7].find("RS=")
                if k < 0 or (k > 0 and c[7][k - 1] != ";"):
                    continue
                v = c[7][k + 3:].split(";", 1)[0].strip()
                for rs in v.split("|"):
                    if rs.isdigit():
                        out.append((int(rs), ucsc_name(c[0]), int(c[1]) - 1))
        return out

    return S._cached("search-clinvar-rs", [e], build, extra="v2")


def gwas_rs(m) -> list[tuple[int, str, int]]:
    """(rs, chrom, pos 0-based) des associations GWAS placées dont SNPS est **un** rsID."""
    from ..fiche import sources as S

    pts = S.gwas(m)
    out = []
    for chrom, pos in pts.pos.items():
        for p, i in zip(pos, pts.idx[chrom]):
            toks = _RS.findall(pts.items[int(i)]["rs"])
            if len(toks) == 1:
                out.append((int(toks[0]), chrom, int(p)))
    return out


def rs_table(clin: list, gw: list, chroms: list[str]) -> np.ndarray:
    """Tableau structuré trié (rs, chrom, pos, sources), une ligne par (rs, chrom, pos)."""
    ci = {c: k for k, c in enumerate(chroms)}
    acc: dict[tuple[int, int, int], int] = {}
    for bit, rows in ((CLINVAR, clin), (GWAS, gw)):
        for rs, chrom, pos in rows:
            if chrom not in ci:
                continue
            key = (rs, ci[chrom], pos)
            acc[key] = acc.get(key, 0) | bit
    a = np.array([(*k, v) for k, v in acc.items()],
                 dtype=[("rs", "<u8"), ("chrom", "u1"), ("pos", "<u4"), ("src", "u1")])
    a.sort(order=("rs", "chrom", "pos"))
    return a


def _rs_block(a: np.ndarray) -> bytes:
    first = a["rs"][0]
    return _deflate(b"".join([(a["rs"] - first).astype("<u4").tobytes(), a["chrom"].astype("u1").tobytes(),
                              a["pos"].astype("<u4").tobytes(), a["src"].astype("u1").tobytes()]))


def decode_rs_block(data: bytes, first: int, n: int) -> np.ndarray:
    raw = _inflate(data)
    d = np.frombuffer(raw, "<u4", n, 0).astype(np.int64) + first
    c = np.frombuffer(raw, "u1", n, 4 * n)
    p = np.frombuffer(raw, "<u4", n, 5 * n)
    s = np.frombuffer(raw, "u1", n, 9 * n)
    out = np.empty(n, dtype=[("rs", "<i8"), ("chrom", "u1"), ("pos", "<u4"), ("src", "u1")])
    out["rs"], out["chrom"], out["pos"], out["src"] = d, c, p, s
    return out


def build(m, out: Path, log=print) -> dict:
    from ..fiche import sources as S

    sizes = [(c, int(n)) for c, n in (ln.split("\t")[:2] for ln in
             S.locked(m, "hg38.chrom.sizes").path.read_text().splitlines() if ln.strip())]
    # Les chromosomes primaires seulement : la grammaire ne lit pas les contigs alternatifs, et le
    # modèle n'en a pas.
    sizes = [(c, n) for c, n in sizes if "_" not in c]
    sizes.sort(key=lambda x: (len(x[0]), x[0]))
    chroms = [c for c, _ in sizes]
    blobs: list[bytes] = []
    offset = 0

    def put(b: bytes) -> dict:
        nonlocal offset
        spec = {"offset": offset, "length": len(b), "sha256": digest(b)}
        blobs.append(b)
        offset += len(b)
        return spec

    genes = gene_rows(m)
    log(f"  gènes : {len(genes) // 2:,} lignes GENCODE, {len(genes):,} clés")
    gblocks = []
    for k in range(0, len(genes), GENE_BLOCK):
        part = genes[k:k + GENE_BLOCK]
        b = _deflate(json.dumps(part, ensure_ascii=True, separators=(",", ":")).encode())
        gblocks.append({"first": part[0][0], "n": len(part), **put(b)})

    clin, gw = clinvar_rs(m), gwas_rs(m)
    rs = rs_table(clin, gw, chroms)
    n_rs = int(len(np.unique(rs["rs"])))
    log(f"  rsID : {n_rs:,} distincts, {len(rs):,} positions (ClinVar {len(clin):,} lignes, GWAS {len(gw):,})")
    rblocks = []
    k = 0
    while k < len(rs):
        j = min(k + RS_BLOCK, len(rs))
        # Un rsID ne chevauche pas deux blocs : le bloc s'allonge jusqu'à sa dernière ligne.
        while j < len(rs) and rs["rs"][j] == rs["rs"][j - 1]:
            j += 1
        # Les écarts au premier rsID du bloc tiennent sur 32 bits (le catalogue GWAS a quelques
        # numéros au-delà de 2³²) : sinon, le bloc s'arrête avant.
        over = np.flatnonzero(rs["rs"][k:j] - rs["rs"][k] >= 2**32)
        if over.size:
            j = k + int(over[0])
        part = rs[k:j]
        rblocks.append({"first": int(part["rs"][0]), "last": int(part["rs"][-1]), "n": int(len(part)),
                        **put(_rs_block(part))})
        k = j

    prov = {key: {"sha256": m[key].sha256, "url": m[key].url}
            for key in ("hg38.chrom.sizes", "gencode.annotation.gtf", "clinvar.vcf", "gwas.catalog")}
    header = {
        "format": "geno-search/1", "version": VERSION, "assembly": "GRCh38", "cell_type": "GM12878",
        "chroms": sizes, "model_chroms": list(MODEL),
        "tables": {
            "genes": {"blocks": gblocks, "rows": len(genes),
                      "columns": ["key", "name", "id", "chrom", "start", "end", "strand", "type"]},
            "rs": {"blocks": rblocks, "rows": int(len(rs)), "distinct": n_rs,
                   "sources": {"clinvar": CLINVAR, "gwas": GWAS},
                   "layout": "Δrs u32[n] · chrom u8[n] (index dans chroms) · pos u32[n] (0-based) · sources u8[n]"},
        },
        "scope": "rsID de ClinVar et du catalogue GWAS seulement — pas dbSNP",
        "provenance": prov,
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
            "genes": len(genes) // 2, "rs": n_rs, "rs_rows": int(len(rs)),
            "blocks": {"genes": len(gblocks), "rs": len(rblocks)}}


class Index:
    """Le lecteur Python : tests, `check`, attendus de la lecture croisée."""

    def __init__(self, path: str | Path):
        self.blob = Path(path).read_bytes()
        magic, version, hlen, hsha = struct.unpack_from("<8sII32s", self.blob)
        if magic != MAGIC:
            raise ValueError("ce n'est pas un index de recherche Geno")
        if version != VERSION:
            raise ValueError(f"index version {version}, ce lecteur lit la {VERSION}")
        hb = self.blob[PREAMBLE:PREAMBLE + hlen]
        if hashlib.sha256(hb).digest() != hsha:
            raise ValueError("en-tête : empreinte divergente")
        self.header = json.loads(_inflate(hb))
        self.base = PREAMBLE + hlen
        self.chroms = [c for c, _ in self.header["chroms"]]
        self.sizes = dict(self.header["chroms"])

    def _blob(self, spec: dict) -> bytes:
        d = self.blob[self.base + spec["offset"]: self.base + spec["offset"] + spec["length"]]
        if digest(d) != spec["sha256"]:
            raise ValueError("bloc : empreinte divergente")
        return d

    def genes(self, key: str) -> list[list]:
        blocks = self.header["tables"]["genes"]["blocks"]
        firsts = [b["first"] for b in blocks]
        # Le dernier bloc qui commence avant la clé peut finir par elle ; les suivants qui
        # commencent par elle la continuent.
        k = max(0, int(np.searchsorted(firsts, key, side="left")) - 1)
        out = []
        while k < len(blocks) and blocks[k]["first"] <= key:
            out += [r for r in json.loads(_inflate(self._blob(blocks[k]))) if r[0] == key]
            k += 1
        return out

    def rs(self, number: int) -> list[tuple[str, int, int]]:
        blocks = self.header["tables"]["rs"]["blocks"]
        firsts = [b["first"] for b in blocks]
        k = int(np.searchsorted(firsts, number, side="right")) - 1
        if k < 0 or number > blocks[k]["last"]:
            return []
        a = decode_rs_block(self._blob(blocks[k]), blocks[k]["first"], blocks[k]["n"])
        sel = a[a["rs"] == number]
        return [(self.chroms[int(c)], int(p), int(s)) for c, p, s in zip(sel["chrom"], sel["pos"], sel["src"])]

    def all_rs(self) -> np.ndarray:
        parts = [decode_rs_block(self._blob(b), b["first"], b["n"]) for b in self.header["tables"]["rs"]["blocks"]]
        return np.concatenate(parts)
