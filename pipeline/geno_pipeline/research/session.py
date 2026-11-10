"""Une session : un `.mcool`, des paramètres, une structure, et le manifeste qui la refait.

La reconstruction est celle de la semaine 5, sans rien d'autre : comptages → équilibrage →
`d = f^(−alpha)` → géodésiques (ShRec3D) → MDS classique. Elle rend **une** structure
consensus, à une isométrie près : l'orientation est fixée par convention (ci-dessous) et la
**chiralité est arbitraire** — des distances ne la portent pas. `alpha` est un paramètre de
l'utilisateur, pas une constante : la semaine 5 a mesuré qu'il dépend de la profondeur.

Reproductibilité : tout ce qui décide du résultat est dans le manifeste — l'empreinte du fichier
d'entrée, les paramètres, les versions. `reproduce` refait la session et compare les **octets**
du `.g3d`. Sur la même machine et les mêmes versions, c'est l'identité ; ailleurs, LAPACK peut
différer au dernier bit, et `reproduce` le dit en donnant l'écart au lieu de le cacher.
"""

from __future__ import annotations

import hashlib
import json
import platform
import re
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ..export import g3d as G
from ..fine.region import monomer_radius
from ..hic.reconstruct import balance_dense, classical_mds, complete_shortest_path, procrustes, to_distance

FORMAT = "geno-session/1"
MAX_BEADS = 3_000          # MDS dense et géodésiques en n³ : au-delà, choisir une résolution plus grossière
BALANCES = ("fichier", "ice", "aucun")


class SessionError(ValueError):
    """Une demande que la session refuse, avec la raison à montrer telle quelle."""


@dataclass(frozen=True)
class Params:
    resolution: int
    region: str                     # « chr2 », « chr2:1-50,000,000 », ou « génome »
    alpha: float = 0.35
    balance: str = "fichier"        # poids du fichier, ICE recalculé, ou comptages bruts
    min_coverage: float = 0.1       # un bin sous cette fraction de la couverture médiane est écarté

    def check(self) -> None:
        if self.balance not in BALANCES:
            raise SessionError(f"équilibrage inconnu : {self.balance!r} (attendu : {', '.join(BALANCES)})")
        if not (0.05 <= self.alpha <= 2.0):
            raise SessionError("alpha doit être entre 0,05 et 2")
        if not (0.0 <= self.min_coverage < 1.0):
            raise SessionError("min_coverage doit être dans [0, 1)")


def sha256_file(path: Path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def versions() -> dict:
    import cooler
    import h5py
    import scipy

    from .. import __version__

    return {"geno": __version__, "python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__, "cooler": cooler.__version__, "h5py": h5py.__version__,
            "machine": platform.machine(), "system": platform.system()}


def inspect(path: Path) -> dict:
    """Ce que le fichier contient : résolutions, chromosomes, équilibrage, profondeur."""
    import cooler

    path = Path(path)
    try:
        groups = cooler.fileops.list_coolers(str(path))
    except Exception as e:
        raise SessionError(f"ce n'est pas un fichier cooler lisible : {e}") from None
    if not groups:
        raise SessionError("aucune carte de contacts dans ce fichier")
    res = []
    for g in groups:
        c = cooler.Cooler(f"{path}::{g}")
        if c.binsize is None:
            continue                 # bins de taille variable : hors de ce que la session reconstruit
        res.append({"resolution": int(c.binsize), "group": g, "bins": int(c.info["nbins"]),
                    "pixels": int(c.info["nnz"]), "contacts": int(c.info.get("sum", 0)),
                    "balanced": "weight" in c.bins().columns})
    if not res:
        raise SessionError("aucune résolution à bins réguliers dans ce fichier")
    res.sort(key=lambda r: r["resolution"])
    c = cooler.Cooler(f"{path}::{res[0]['group']}")
    return {"name": path.name, "bytes": path.stat().st_size, "assembly": c.info.get("genome-assembly"),
            "chroms": [[str(n), int(s)] for n, s in zip(c.chromnames, c.chromsizes)],
            "resolutions": res, "max_beads": MAX_BEADS}


def _cooler(path: Path, resolution: int):
    import cooler

    for r in inspect(path)["resolutions"]:
        if r["resolution"] == resolution:
            return cooler.Cooler(f"{path}::{r['group']}")
    raise SessionError(f"résolution {resolution:,} pb absente du fichier")


def _extent(clr, region: str) -> tuple[str, list[tuple[str, int, int]]]:
    """La région demandée → son nom canonique et ses intervalles par chromosome."""
    sizes = dict(zip(clr.chromnames, (int(s) for s in clr.chromsizes)))
    r = region.strip()
    if r.lower() in ("génome", "genome", "*"):
        return "génome", [(c, 0, s) for c, s in sizes.items()]
    if r in sizes:
        return r, [(r, 0, sizes[r])]
    # Le nom du chromosome est celui du fichier, quel qu'il soit (« chr2 », « 2 », « chrP1 ») ;
    # les bornes s'écrivent comme partout : 1-based, incluses.
    m = re.fullmatch(r"(.+?)\s*:\s*([0-9][0-9,_ ]*)\s*[-\u2013]\s*([0-9][0-9,_ ]*)", r)
    if m and m.group(1) in sizes:
        c = m.group(1)
        a, b = (int(re.sub(r"[,_ ]", "", m.group(k))) for k in (2, 3))
        if a < 1 or b < a:
            raise SessionError("bornes illisibles : début ≥ 1 et fin ≥ début")
        if b > sizes[c]:
            raise SessionError(f"{c} finit à {sizes[c]:,}")
        return f"{c}:{a}-{b}", [(c, a - 1, b)]
    raise SessionError(f"région inconnue du fichier : {region!r} — un chromosome ({', '.join(list(sizes)[:4])}…), "
                       "chrom:début-fin, ou « génome »")


def orient(x: np.ndarray) -> np.ndarray:
    """Orientation de convention : centré, axes du MDS (variance décroissante), et sur chaque axe
    le signe qui rend positive la coordonnée de plus grande valeur absolue. Sans cela, le signe
    d'un vecteur propre est celui que LAPACK a rendu ce jour-là."""
    x = x - x.mean(axis=0)
    for k in range(x.shape[1]):
        if x[np.argmax(np.abs(x[:, k])), k] < 0:
            x[:, k] = -x[:, k]
    return x


@dataclass
class Result:
    positions_nm: np.ndarray
    chrom: np.ndarray               # index dans `chroms`, par bille
    start: np.ndarray
    end: np.ndarray
    chroms: list[str]
    stats: dict
    region: str


def reconstruct(path: Path, p: Params, progress=lambda step, frac: None) -> Result:
    p.check()
    clr = _cooler(Path(path), p.resolution)
    name, extent = _extent(clr, p.region)
    progress("lecture", 0.05)
    bins = clr.bins()[:]
    keep = np.zeros(len(bins), bool)
    for c, s, e in extent:
        keep |= ((bins["chrom"].astype(str) == c) & (bins["start"] < e) & (bins["end"] > s)).to_numpy()
    idx = np.flatnonzero(keep)
    if len(idx) > MAX_BEADS:
        raise SessionError(f"{len(idx):,} bins à {p.resolution:,} pb : au-delà de {MAX_BEADS:,} billes, la "
                           "reconstruction dense ne tient pas — choisir une résolution plus grossière ou une région plus courte")
    if len(idx) < 8:
        raise SessionError(f"{len(idx)} bins seulement : trop peu pour une structure")
    if p.balance == "fichier" and "weight" not in bins.columns:
        raise SessionError("le fichier n'a pas de poids d'équilibrage à cette résolution : choisir « ice » ou « aucun »")

    # Les pixels de l'étendue, en dense. `matrix` par blocs : un fetch par paire de chromosomes
    # relirait le fichier ; l'étendue est contiguë en bins pour un chromosome, sinon on indexe.
    lo, hi = int(idx[0]), int(idx[-1]) + 1
    raw = clr.matrix(balance=False, sparse=True)[lo:hi, lo:hi].toarray().astype(np.float64)
    sel = idx - lo
    raw = raw[np.ix_(sel, sel)]
    np.fill_diagonal(raw, 0.0)
    progress("équilibrage", 0.2)

    cov = raw.sum(axis=1)
    med = float(np.median(cov[cov > 0])) if (cov > 0).any() else 0.0
    good = cov > p.min_coverage * med if med > 0 else np.zeros(len(cov), bool)
    if p.balance == "fichier":
        w = bins["weight"].to_numpy()[idx]
        good &= np.isfinite(w)
        m = raw * np.where(np.isfinite(w), w, 0.0)[:, None] * np.where(np.isfinite(w), w, 0.0)[None, :]
    elif p.balance == "ice":
        m = np.zeros_like(raw)
        g = np.flatnonzero(good)
        m[np.ix_(g, g)] = balance_dense(raw[np.ix_(g, g)])
    else:
        m = raw
    g = np.flatnonzero(good)
    if len(g) < 8:
        raise SessionError("presque aucun bin couvert dans cette région")
    m = m[np.ix_(g, g)]
    m = (m + m.T) / 2

    progress("distances", 0.35)
    target = to_distance(m, p.alpha)
    observed = np.isfinite(target) & ~np.eye(len(g), dtype=bool)
    # Des morceaux sans aucun contact entre eux (deux chromosomes d'une carte sans trans) : la
    # complétion les pose à une distance de convention, et leur position relative ne veut rien dire.
    from scipy.sparse.csgraph import connected_components

    pieces = int(connected_components(observed, directed=False)[0])
    geo = complete_shortest_path(target)
    progress("coordonnées", 0.7)
    x = orient(classical_mds(geo))

    # Échelle : unités arbitraires → nm. Deux billes voisines sur la chaîne se touchent, au rayon
    # que la loi de la semaine 6 donne à cette quantité de chromatine. C'est une convention
    # d'affichage, dite comme telle : les contacts ne portent pas d'échelle absolue.
    radius = monomer_radius(p.resolution)
    b = bins.iloc[idx[g]]
    chroms = [c for c, _, _ in extent]
    ci = np.array([chroms.index(str(c)) for c in b["chrom"].astype(str)], np.int64)
    adj = np.flatnonzero((ci[1:] == ci[:-1]) & (b["start"].to_numpy()[1:] == b["end"].to_numpy()[:-1]))
    step = np.linalg.norm(x[adj + 1] - x[adj], axis=1)
    scale = 2.0 * radius / float(np.median(step)) if len(step) and np.median(step) > 0 else 1.0
    x = x * scale

    progress("mesures", 0.9)
    d = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)
    iu = np.triu_indices(len(g), 1)
    ok = observed[iu]
    with np.errstate(divide="ignore", invalid="ignore"):
        lt, ld = np.log(target[iu][ok]), np.log(d[iu][ok])
    fin = np.isfinite(lt) & np.isfinite(ld)
    stats = {
        "bins": int(len(idx)), "beads": int(len(g)), "masked": int(len(idx) - len(g)),
        "contacts": float(raw.sum() / 2), "observed_pairs": float(ok.mean()), "pieces": pieces,
        # Accord entre la structure et ce qu'on lui demandait : distances cibles (contacts observés
        # seulement) contre distances de la structure, en log. Ce n'est pas une preuve de justesse :
        # une structure peut reproduire sa carte et être fausse (une cellule n'est pas une moyenne).
        "r_log_distance": float(np.corrcoef(lt[fin], ld[fin])[0, 1]) if fin.sum() > 2 else float("nan"),
        "radius_nm": float(radius), "scale_nm_per_unit": float(scale),
        "gyration_nm": float(np.sqrt((x**2).sum(axis=1).mean())),
    }
    return Result(x, ci, b["start"].to_numpy().astype(np.int64), b["end"].to_numpy().astype(np.int64),
                  chroms, stats, name)


def write_g3d(res: Result, p: Params, source: dict, out: Path) -> dict:
    radius = res.stats["radius_nm"]
    level = G.Level(
        name="structure", bp_per_bead=p.resolution, positions_nm=res.positions_nm.astype(np.float64),
        copy=res.chrom, start=res.start, end=res.end,
        variability_nm=np.zeros(len(res.start)), copies=res.chroms, radius_nm=[radius] * len(res.chroms),
        evidence="reconstructed", precision_nm=0.01 * radius,
        variability={"kind": "none", "why": "une structure consensus, pas un ensemble"},
        notes=["reconstruite depuis une carte de contacts importée : équilibrage, d = f^(−alpha), "
               "géodésiques (ShRec3D), MDS classique",
               "une structure consensus à une isométrie près : orientation de convention, chiralité arbitraire",
               "échelle de convention : deux billes voisines se touchent (rayon de la loi de la semaine 6)"],
    )
    return G.write(out, [level], header={"assembly": source.get("assembly"), "cell_type": None,
                                         "session": {"input": source["name"], "sha256": source["sha256"],
                                                     "params": asdict(p), "region": res.region}})


def run(mcool: Path, p: Params, out_dir: Path, progress=lambda step, frac: None, input_sha256: str | None = None) -> dict:
    """Reconstruit, écrit `structure.g3d` et `session.json` dans `out_dir`, rend le manifeste."""
    mcool = Path(mcool)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    info = inspect(mcool)
    sha = input_sha256 or sha256_file(mcool)
    res = reconstruct(mcool, p, progress)
    g3d = out_dir / "structure.g3d"
    write_g3d(res, p, {"name": mcool.name, "sha256": sha, "assembly": info["assembly"]}, g3d)
    manifest = {
        "format": FORMAT,
        "input": {"name": mcool.name, "bytes": info["bytes"], "sha256": sha, "assembly": info["assembly"]},
        "params": asdict(p),
        "region": res.region,
        "method": {"chain": ["équilibrage", "d = f^(−alpha)", "géodésiques (ShRec3D, Lesne 2014)", "MDS classique (Torgerson)"],
                   "orientation": "centrée, axes du MDS, signe : la plus grande |coordonnée| de chaque axe est positive",
                   "chirality": "arbitraire — des distances ne la déterminent pas",
                   "scale": "deux billes voisines se touchent ; rayon par la loi de la semaine 6"},
        "versions": versions(),
        "output": {"file": g3d.name, "bytes": g3d.stat().st_size, "sha256": sha256_file(g3d),
                   "positions_sha256": hashlib.sha256(np.ascontiguousarray(res.positions_nm, "<f8").tobytes()).hexdigest()},
        "stats": res.stats,
        # Hors de ce qui se reproduit : la date et la durée de cette exécution-ci.
        "run": {"created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t0)), "seconds": round(time.time() - t0, 2)},
    }
    (out_dir / "session.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=1))
    progress("fait", 1.0)
    return manifest


def reproduce(manifest: dict, mcool: Path, out_dir: Path, progress=lambda step, frac: None) -> dict:
    """Refait la session d'un manifeste sur `mcool`, et compare. Refuse un fichier qui n'est pas
    celui du manifeste : reproduire autre chose n'est pas reproduire."""
    if manifest.get("format") != FORMAT:
        raise SessionError(f"ce n'est pas un manifeste {FORMAT}")
    sha = sha256_file(Path(mcool))
    if sha != manifest["input"]["sha256"]:
        raise SessionError("le fichier fourni n'est pas celui du manifeste (empreinte sha256 différente)")
    try:
        p = Params(**manifest["params"])
    except TypeError as e:
        raise SessionError(f"paramètres du manifeste illisibles : {e}") from None
    new = run(mcool, p, out_dir, progress, input_sha256=sha)
    identical = new["output"]["sha256"] == manifest["output"]["sha256"]
    out = {"identical": identical, "expected_sha256": manifest["output"]["sha256"], "sha256": new["output"]["sha256"],
           "versions_differ": {k: [manifest["versions"].get(k), v] for k, v in new["versions"].items()
                               if manifest["versions"].get(k) != v},
           "manifest": new}
    if not identical:
        out["note"] = ("les octets diffèrent : mêmes paramètres et même fichier, donc une différence de versions ou de "
                       "bibliothèque d'algèbre linéaire — voir versions_differ")
    return out


def compare(a: np.ndarray, b: np.ndarray) -> float:
    """Écart entre deux structures après alignement (réflexion permise), en fraction du rayon de giration."""
    _, rmsd = procrustes(a, b)
    return rmsd / float(np.sqrt(((b - b.mean(axis=0)) ** 2).sum(axis=1).mean()))
