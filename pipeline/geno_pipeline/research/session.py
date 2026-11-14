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
from . import extras as X
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
    ensemble: int = 0               # semaine 18 : réplicats par rééchantillonnage de la carte (0 : aucun)
    seed: int = 0

    def as_dict(self) -> dict:
        """Sans ensemble, le dictionnaire de la semaine 17 tel quel : ses manifestes se rejouent
        encore à l'octet près."""
        d = asdict(self)
        if not self.ensemble:
            d.pop("ensemble")
            d.pop("seed")
        return d

    def check(self) -> None:
        if not (0 <= self.ensemble <= X.ENSEMBLE_MAX):
            raise SessionError(f"ensemble : entre 0 et {X.ENSEMBLE_MAX} réplicats")
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
    ctx: dict | None = None     # les comptages lus, pour l'ensemble et les statistiques (semaine 18)


def _read(path: Path, p: Params):
    """Les comptages bruts de la région, en dense, et de quoi les situer."""
    p.check()
    clr = _cooler(Path(path), p.resolution)
    name, extent = _extent(clr, p.region)
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
    lo, hi = int(idx[0]), int(idx[-1]) + 1
    raw = clr.matrix(balance=False, sparse=True)[lo:hi, lo:hi].toarray().astype(np.float64)
    sel = idx - lo
    raw = raw[np.ix_(sel, sel)]
    np.fill_diagonal(raw, 0.0)
    w = bins["weight"].to_numpy()[idx] if "weight" in bins.columns else None
    cov = raw.sum(axis=1)
    med = float(np.median(cov[cov > 0])) if (cov > 0).any() else 0.0
    good = cov > p.min_coverage * med if med > 0 else np.zeros(len(cov), bool)
    if p.balance == "fichier":
        good &= np.isfinite(w)
    g = np.flatnonzero(good)
    if len(g) < 8:
        raise SessionError("presque aucun bin couvert dans cette région")
    b = bins.iloc[idx[g]]
    chroms = [c for c, _, _ in extent]
    ci = np.array([chroms.index(str(c)) for c in b["chrom"].astype(str)], np.int64)
    return {"raw": raw, "w": w, "g": g, "name": name, "chroms": chroms, "ci": ci, "bins": int(len(idx)),
            "start": b["start"].to_numpy().astype(np.int64), "end": b["end"].to_numpy().astype(np.int64)}


def _solve(raw: np.ndarray, ctx: dict, p: Params):
    """Comptages (toute la région) → coordonnées en nm des billes gardées, orientées et mises à l'échelle."""
    g, w = ctx["g"], ctx["w"]
    if p.balance == "fichier":
        ww = np.where(np.isfinite(w), w, 0.0)
        m = (raw * ww[:, None] * ww[None, :])[np.ix_(g, g)]
    elif p.balance == "ice":
        m = balance_dense(raw[np.ix_(g, g)])
    else:
        m = raw[np.ix_(g, g)]
    m = (m + m.T) / 2
    target = to_distance(m, p.alpha)
    observed = np.isfinite(target) & ~np.eye(len(g), dtype=bool)
    x = orient(classical_mds(complete_shortest_path(target)))
    # Échelle : unités arbitraires → nm. Deux billes voisines sur la chaîne se touchent, au rayon
    # que la loi de la semaine 6 donne à cette quantité de chromatine. C'est une convention
    # d'affichage, dite comme telle : les contacts ne portent pas d'échelle absolue.
    radius = monomer_radius(p.resolution)
    ci = ctx["ci"]
    adj = np.flatnonzero((ci[1:] == ci[:-1]) & (ctx["start"][1:] == ctx["end"][:-1]))
    step = np.linalg.norm(x[adj + 1] - x[adj], axis=1)
    scale = 2.0 * radius / float(np.median(step)) if len(step) and np.median(step) > 0 else 1.0
    return x * scale, target, observed, radius, scale


def reconstruct(path: Path, p: Params, progress=lambda step, frac: None) -> Result:
    progress("lecture", 0.05)
    ctx = _read(path, p)
    raw, g = ctx["raw"], ctx["g"]
    progress("distances", 0.3)
    x, target, observed, radius, scale = _solve(raw, ctx, p)
    # Des morceaux sans aucun contact entre eux (deux chromosomes d'une carte sans trans) : la
    # complétion les pose à une distance de convention, et leur position relative ne veut rien dire.
    from scipy.sparse.csgraph import connected_components

    pieces = int(connected_components(observed, directed=False)[0])
    progress("mesures", 0.6)
    d = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)
    iu = np.triu_indices(len(g), 1)
    ok = observed[iu]
    with np.errstate(divide="ignore", invalid="ignore"):
        lt, ld = np.log(target[iu][ok]), np.log(d[iu][ok])
    fin = np.isfinite(lt) & np.isfinite(ld)
    stats = {
        "bins": ctx["bins"], "beads": int(len(g)), "masked": int(ctx["bins"] - len(g)),
        "contacts": float(raw.sum() / 2), "observed_pairs": float(ok.mean()), "pieces": pieces,
        # Accord entre la structure et ce qu'on lui demandait : distances cibles (contacts observés
        # seulement) contre distances de la structure, en log. Ce n'est pas une preuve de justesse :
        # une structure peut reproduire sa carte et être fausse (une cellule n'est pas une moyenne).
        "r_log_distance": float(np.corrcoef(lt[fin], ld[fin])[0, 1]) if fin.sum() > 2 else float("nan"),
        "radius_nm": float(radius), "scale_nm_per_unit": float(scale),
        "gyration_nm": float(np.sqrt((x**2).sum(axis=1).mean())),
    }
    res = Result(x, ctx["ci"], ctx["start"], ctx["end"], ctx["chroms"], stats, ctx["name"])
    res.ctx = ctx
    return res


def write_g3d(res: Result, p: Params, source: dict, out: Path, variability: np.ndarray | None = None) -> dict:
    radius = res.stats["radius_nm"]
    level = G.Level(
        name="structure", bp_per_bead=p.resolution, positions_nm=res.positions_nm.astype(np.float64),
        copy=res.chrom, start=res.start, end=res.end,
        variability_nm=np.zeros(len(res.start)) if variability is None else variability,
        copies=res.chroms, radius_nm=[radius] * len(res.chroms),
        evidence="reconstructed", precision_nm=0.01 * radius,
        variability={"kind": "none", "why": "une structure consensus, pas un ensemble"} if variability is None else
        {"kind": "bootstrap_rms", "replicates": p.ensemble, "seed": p.seed,
         "what": "dispersion sous rééchantillonnage des comptes de la carte — bruit de comptage, pas variabilité entre cellules"},
        notes=["reconstruite depuis une carte de contacts importée : équilibrage, d = f^(−alpha), "
               "géodésiques (ShRec3D), MDS classique",
               "une structure consensus à une isométrie près : orientation de convention, chiralité arbitraire",
               "échelle de convention : deux billes voisines se touchent (rayon de la loi de la semaine 6)"],
    )
    return G.write(out, [level], header={"assembly": source.get("assembly"), "cell_type": None,
                                         "session": {"input": source["name"], "sha256": source["sha256"],
                                                     "params": p.as_dict(), "region": res.region}})


def is_hic(path: Path) -> bool:
    with open(path, "rb") as f:
        return f.read(3) == b"HIC"


def as_cooler(path: Path, cache: Path) -> Path:
    """Un `.hic` (Juicer) est converti une fois en `.mcool` par hic2cool, rangé sous son empreinte ;
    un fichier cooler est rendu tel quel. L'empreinte d'une session reste celle du fichier **donné**."""
    path = Path(path)
    if not is_hic(path):
        return path
    try:
        from hic2cool import hic2cool_convert
    except ImportError:
        raise SessionError("lire un .hic demande hic2cool (pip install hic2cool)") from None
    out = Path(cache) / f"{sha256_file(path)[:24]}.mcool"
    if not out.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
        tmp = out.with_suffix(".part.mcool")
        try:
            hic2cool_convert(str(path), str(tmp), 0, silent=True)
        except BaseException as e:                         # hic2cool appelle sys.exit sur un fichier qu'il ne lit pas
            tmp.unlink(missing_ok=True)
            raise SessionError(f"ce .hic n'a pas pu être converti : {e}") from None
        tmp.replace(out)
    return out


def run(mcool: Path, p: Params, out_dir: Path, progress=lambda step, frac: None, input_sha256: str | None = None,
        source: Path | None = None) -> dict:
    """Reconstruit, écrit `structure.g3d`, ses exports, ses statistiques, son ensemble s'il est
    demandé, et `session.json` dans `out_dir` ; rend le manifeste. `source` : le fichier donné par
    l'utilisateur quand `mcool` en est la conversion (un `.hic`)."""
    mcool = Path(mcool)
    source = Path(source) if source else mcool
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    info = inspect(mcool)
    sha = input_sha256 or sha256_file(source)
    res = reconstruct(mcool, p, progress)
    ctx = res.ctx
    variability = None
    files: dict[str, Path] = {}
    if p.ensemble:
        if res.stats["beads"] > X.ENSEMBLE_MAX_BEADS:
            raise SessionError(f"ensemble : {res.stats['beads']:,} billes, au-delà de {X.ENSEMBLE_MAX_BEADS:,} — "
                               "l'ensemble ne se calcule qu'en résolution grossière")
        reps, variability = X.ensemble(res.positions_nm, ctx, p, _solve, progress)
        files["ensemble"] = out_dir / "ensemble.f32"
        files["ensemble"].write_bytes(np.ascontiguousarray(reps, "<f4").tobytes())
    g3d = out_dir / "structure.g3d"
    write_g3d(res, p, {"name": source.name, "sha256": sha, "assembly": info["assembly"]}, g3d, variability)

    # Exports et statistiques depuis les positions **du fichier** : ce qu'on exporte est ce qu'on voit.
    lv = G.read(g3d).levels["structure"]
    x = np.asarray(lv["positions_nm"], np.float64)
    var = np.asarray(lv["variability_nm"], np.float64)
    title = f"{source.name} {res.region} {p.resolution} pb"
    texts = {"pdb": ("structure.pdb", X.pdb(x, res.chrom, res.chroms, res.start, var, title)),
             "3dg": ("structure.3dg", X.threedg(x, res.chrom, res.chroms, res.start)),
             "gltf": ("structure.gltf", X.gltf(x, res.chrom, res.chroms, {"input": source.name, "region": res.region,
                                                                             "resolution": p.resolution, "evidence": "reconstructed"})),
             "statistics": ("statistics.json", json.dumps(X.statistics(x, ctx, p.resolution, variability), ensure_ascii=False, indent=1))}
    for key, (name, text) in texts.items():
        files[key] = out_dir / name
        files[key].write_text(text)
    manifest = {
        "format": FORMAT,
        "input": {"name": source.name, "bytes": source.stat().st_size, "sha256": sha, "assembly": info["assembly"]},
        "params": p.as_dict(),
        "region": res.region,
        "method": {"chain": ["équilibrage", "d = f^(−alpha)", "géodésiques (ShRec3D, Lesne 2014)", "MDS classique (Torgerson)"],
                   "orientation": "centrée, axes du MDS, signe : la plus grande |coordonnée| de chaque axe est positive",
                   "chirality": "arbitraire — des distances ne la déterminent pas",
                   "scale": "deux billes voisines se touchent ; rayon par la loi de la semaine 6"},
        "versions": versions(),
        "output": {"file": g3d.name, "bytes": g3d.stat().st_size, "sha256": sha256_file(g3d),
                   "positions_sha256": hashlib.sha256(np.ascontiguousarray(res.positions_nm, "<f8").tobytes()).hexdigest()},
        # Semaine 18 : ce qui accompagne la structure, chaque fichier avec son empreinte.
        "files": {k: {"file": f.name, "bytes": f.stat().st_size, "sha256": sha256_file(f)} for k, f in files.items()},
        "stats": res.stats,
        # Hors de ce qui se reproduit : la date et la durée de cette exécution-ci.
        "run": {"created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t0)), "seconds": round(time.time() - t0, 2)},
    }
    if p.ensemble:
        manifest["method"]["ensemble"] = ("rééchantillonnage de Poisson des comptes, graine dans params ; réplicats alignés sur le "
                                          "consensus (rotation, réflexion) ; variabilité = écart quadratique moyen par bille")
        manifest["ensemble"] = {"replicates": p.ensemble, "beads": res.stats["beads"], "dtype": "<f4", "unit": "nm",
                                "layout": "réplicat, bille, xyz"}
    if source != mcool:
        import hic2cool

        manifest["input"]["converted_by"] = f"hic2cool {hic2cool.__version__}"
        manifest["versions"]["hic2cool"] = hic2cool.__version__
    (out_dir / "session.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=1))
    progress("fait", 1.0)
    return manifest


def reproduce(manifest: dict, mcool: Path, out_dir: Path, progress=lambda step, frac: None, cache: Path | None = None) -> dict:
    """Refait la session d'un manifeste sur `mcool`, et compare. Refuse un fichier qui n'est pas
    celui du manifeste : reproduire autre chose n'est pas reproduire."""
    if manifest.get("format") != FORMAT:
        raise SessionError(f"ce n'est pas un manifeste {FORMAT}")
    source = Path(mcool)
    sha = sha256_file(source)
    if sha != manifest["input"]["sha256"]:
        raise SessionError("le fichier fourni n'est pas celui du manifeste (empreinte sha256 différente)")
    try:
        p = Params(**manifest["params"])
    except TypeError as e:
        raise SessionError(f"paramètres du manifeste illisibles : {e}") from None
    new = run(as_cooler(source, cache or out_dir.parent / "converted"), p, out_dir, progress, input_sha256=sha, source=source)
    # Tout ce que le manifeste d'origine décrit doit revenir à l'octet près ; un manifeste de la
    # semaine 17 ne décrit que la structure.
    differ = [k for k, f in manifest.get("files", {}).items() if new["files"].get(k, {}).get("sha256") != f["sha256"]]
    identical = new["output"]["sha256"] == manifest["output"]["sha256"] and not differ
    out = {"identical": identical, "expected_sha256": manifest["output"]["sha256"], "sha256": new["output"]["sha256"],
           "files_compared": sorted(manifest.get("files", {})), "files_differ": differ,
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
