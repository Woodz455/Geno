"""Semaine 18 : ce qu'une session ajoute à sa structure — un ensemble, des statistiques, des exports.

**Ensemble par rééchantillonnage.** La carte est un tirage : chaque pixel est un compte. On retire
`n` cartes (Poisson autour des comptes observés, graine dans le manifeste), on reconstruit chacune,
on l'aligne sur la structure consensus (rotation et réflexion, sans échelle), et la dispersion de
chaque bille devient sa **variabilité**. Ce que cela mesure : ce que le **bruit de comptage** de la
carte laisse indéterminé. Ce que cela ne mesure pas : la variabilité d'une cellule à l'autre —
une carte de population ne la contient pas sous cette forme (R1).

**Exports.** PDB (un atome CA par bille, une chaîne par chromosome), 3DG (le format texte de
hickit : chrom, position, x, y, z), glTF 2.0 (points et lignes, tampon embarqué). Tous écrits
depuis les positions **du `.g3d`** : ce qu'on exporte est ce qu'on voit.
"""

from __future__ import annotations

import base64
import json

import numpy as np

ENSEMBLE_MAX_BEADS = 1_500       # l'ensemble ne se calcule qu'en résolution grossière (feuille de route)
ENSEMBLE_MAX = 64


def rigid(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """`x` posé sur `y` : translation, rotation, réflexion — pas d'échelle (elle est de convention)."""
    xc, yc = x - x.mean(axis=0), y - y.mean(axis=0)
    u, _, vt = np.linalg.svd(xc.T @ yc)
    return xc @ (u @ vt) + y.mean(axis=0)


def resample(raw: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    up = np.triu(rng.poisson(np.triu(raw, 1)).astype(np.float64), 1)
    return up + up.T


def ensemble(consensus: np.ndarray, ctx: dict, p, solve, progress=lambda step, frac: None) -> tuple[np.ndarray, np.ndarray]:
    """(réplicats alignés (n, billes, 3) en nm, variabilité par bille en nm)."""
    out = np.empty((p.ensemble, len(consensus), 3))
    for k in range(p.ensemble):
        rng = np.random.default_rng([p.seed, k])
        x = solve(resample(ctx["raw"], rng), ctx, p)[0]
        out[k] = rigid(x, consensus)
        progress(f"ensemble {k + 1}/{p.ensemble}", 0.6 + 0.35 * (k + 1) / p.ensemble)
    var = np.sqrt(((out - consensus[None]) ** 2).sum(axis=2).mean(axis=0))
    return out, var


def _curve(sep: np.ndarray, val: np.ndarray, per_decade: int = 8) -> list[list[float]]:
    """Moyenne par classes logarithmiques de séparation : [[séparation en bins, valeur, paires], …]."""
    ok = (sep > 0) & np.isfinite(val)
    sep, val = sep[ok], val[ok]
    if not len(sep):
        return []
    k = np.floor(np.log10(sep) * per_decade).astype(int)
    out = []
    for b in np.unique(k):
        m = k == b
        out.append([float(sep[m].mean()), float(val[m].mean()), int(m.sum())])
    return out


def _slope(curve: list[list[float]], lo: float, hi: float) -> float | None:
    pts = [(s, v) for s, v, _ in curve if lo <= s <= hi and v > 0]
    if len(pts) < 3:
        return None
    x, y = np.log([a for a, _ in pts]), np.log([b for _, b in pts])
    return float(np.polyfit(x, y, 1)[0])


def statistics(x: np.ndarray, ctx: dict, resolution: int, variability: np.ndarray | None) -> dict:
    """La carte et la structure, côte à côte : P(s) des comptes lus, R(s) de la structure, par
    paires du même chromosome ; rayon de giration par chromosome ; variabilité si ensemble."""
    g, ci = ctx["g"], ctx["ci"]
    raw = ctx["raw"][np.ix_(g, g)]
    bin_i = ctx["start"] // resolution
    i, j = np.triu_indices(len(g), 1)
    same = ci[i] == ci[j]
    i, j = i[same], j[same]
    sep = (bin_i[j] - bin_i[i]).astype(float)
    ps = _curve(sep, raw[i, j])
    rs = _curve(sep, np.linalg.norm(x[i] - x[j], axis=1))
    top = float(sep.max()) if len(sep) else 0.0
    lo, hi = 2.0, max(4.0, top / 4)
    out = {
        "resolution": resolution, "fit_range_bins": [lo, hi],
        "contacts_vs_separation": ps, "distance_vs_separation": rs,
        # P(s) ∝ s^slope ; R(s) ∝ s^nu. Lus sur la même plage, hors des deux premiers bins et du
        # dernier quart, où les moyennes portent sur trop peu de paires.
        "ps_slope": _slope(ps, lo, hi), "rs_nu": _slope(rs, lo, hi),
        "chromosomes": [{"name": c, "beads": int((ci == k).sum()),
                         "gyration_nm": float(np.sqrt(((x[ci == k] - x[ci == k].mean(axis=0)) ** 2).sum(axis=1).mean())) if (ci == k).any() else None}
                        for k, c in enumerate(ctx["chroms"])],
    }
    if variability is not None:
        out["variability_nm"] = {"median": float(np.median(variability)), "p90": float(np.percentile(variability, 90)),
                                 "max": float(variability.max())}
    return out


# — Exports —

def _chain_ids(n: int) -> list[str]:
    abc = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    return [abc[k % len(abc)] for k in range(n)]


def pdb(x_nm: np.ndarray, ci: np.ndarray, chroms: list[str], start: np.ndarray, variability: np.ndarray, title: str) -> str:
    """Un CA par bille. Les colonnes de coordonnées du PDB (8 caractères, 3 décimales) tiennent de
    −999,999 à 9999,999 : l'unité est le nm, ou un multiple de 10 nm si la structure n'y tient pas
    — dite dans REMARK 2, jamais devinée."""
    unit = 1.0
    while np.abs(x_nm).max() / unit >= 999:
        unit *= 10.0
    ids = _chain_ids(len(chroms))
    bmax = max(1.0, float(variability.max()) / 999.0) if len(variability) else 1.0
    lines = [f"TITLE     {title[:60]}",
             "REMARK   1 GENO — STRUCTURE RECONSTRUITE (ShRec3D + MDS), CHIRALITE ARBITRAIRE",
             f"REMARK   2 UNIT_NM {unit:g}  (1 unité de coordonnée = {unit:g} nm)",
             f"REMARK   3 BFACTOR_NM {bmax:g}  (facteur B × {bmax:g} = variabilité de la bille, nm)"]
    for k, c in enumerate(chroms):
        lines.append(f"REMARK   4 CHAIN {ids[k]} = {c}")
    serial = 0
    for k in range(len(chroms)):
        idx = np.flatnonzero(ci == k)
        for r, i in enumerate(idx, 1):
            serial += 1
            a, b, c = x_nm[i] / unit
            lines.append(f"ATOM  {serial % 100000:5d}  CA  GLY {ids[k]}{r % 10000:4d}    {a:8.3f}{b:8.3f}{c:8.3f}  1.00{variability[i] / bmax:6.2f}           C")
        if len(idx):
            lines.append(f"TER   {(serial + 1) % 100000:5d}      GLY {ids[k]}{len(idx) % 10000:4d}")
    lines.append("END")
    return "\n".join(lines) + "\n"


def threedg(x_nm: np.ndarray, ci: np.ndarray, chroms: list[str], start: np.ndarray) -> str:
    lines = ["#unit: nm", "#columns: chrom\tposition (début du bin, 0-based)\tx\ty\tz"]
    for i in range(len(x_nm)):
        lines.append(f"{chroms[ci[i]]}\t{int(start[i])}\t{x_nm[i, 0]:.3f}\t{x_nm[i, 1]:.3f}\t{x_nm[i, 2]:.3f}")
    return "\n".join(lines) + "\n"


def gltf(x_nm: np.ndarray, ci: np.ndarray, chroms: list[str], meta: dict) -> str:
    """glTF 2.0 autonome : un nœud par chromosome, une ligne brisée (la chaîne) et ses points.
    Positions en micromètres (float32), tampon embarqué en base64."""
    pos = (x_nm / 1000.0).astype("<f4")
    buf = pos.tobytes()
    views = [{"buffer": 0, "byteOffset": 0, "byteLength": len(buf), "target": 34962}]
    accessors, meshes, nodes, materials = [], [], [], []
    for k, c in enumerate(chroms):
        idx = np.flatnonzero(ci == k)
        if not len(idx):
            continue
        if not np.array_equal(idx, np.arange(idx[0], idx[0] + len(idx))):
            raise ValueError("billes d'un chromosome non contiguës")
        sub = pos[idx]
        accessors.append({"bufferView": 0, "byteOffset": int(idx[0]) * 12, "componentType": 5126, "count": int(len(idx)),
                          "type": "VEC3", "min": [float(v) for v in sub.min(axis=0)], "max": [float(v) for v in sub.max(axis=0)]})
        h = k / max(1, len(chroms)) * 6.283 + 0.6
        materials.append({"name": c, "pbrMetallicRoughness": {"baseColorFactor": [0.6 + 0.4 * float(np.cos(h)), 0.6 + 0.4 * float(np.cos(h - 2.094)),
                                                                                0.6 + 0.4 * float(np.cos(h + 2.094)), 1.0], "metallicFactor": 0.0}})
        a = len(accessors) - 1
        meshes.append({"name": c, "primitives": [{"attributes": {"POSITION": a}, "mode": 3, "material": len(materials) - 1},
                                                 {"attributes": {"POSITION": a}, "mode": 0, "material": len(materials) - 1}]})
        nodes.append({"name": c, "mesh": len(meshes) - 1})
    doc = {"asset": {"version": "2.0", "generator": "Geno", "extras": {"unit": "micromètre", **meta}},
           "scene": 0, "scenes": [{"nodes": list(range(len(nodes)))}], "nodes": nodes, "meshes": meshes, "materials": materials,
           "accessors": accessors, "bufferViews": views,
           "buffers": [{"byteLength": len(buf), "uri": "data:application/octet-stream;base64," + base64.b64encode(buf).decode()}]}
    return json.dumps(doc, ensure_ascii=False, separators=(",", ":"))


def parse_pdb(text: str) -> tuple[np.ndarray, float, list[str]]:
    unit = 1.0
    xyz, chains = [], []
    for l in text.splitlines():
        if l.startswith("REMARK   2 UNIT_NM"):
            unit = float(l.split()[3])
        elif l.startswith("ATOM"):
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            chains.append(l[21])
    return np.array(xyz) * unit, unit, chains


def parse_gltf(text: str) -> list[np.ndarray]:
    d = json.loads(text)
    raw = base64.b64decode(d["buffers"][0]["uri"].split(",", 1)[1])
    assert len(raw) == d["buffers"][0]["byteLength"] == d["bufferViews"][0]["byteLength"]
    out = []
    for m in d["meshes"]:
        a = d["accessors"][m["primitives"][0]["attributes"]["POSITION"]]
        v = np.frombuffer(raw, "<f4", a["count"] * 3, a["byteOffset"]).reshape(-1, 3)
        assert np.allclose(v.min(axis=0), a["min"]) and np.allclose(v.max(axis=0), a["max"])
        out.append(v)
    return out


__all__ = ["ensemble", "statistics", "pdb", "threedg", "gltf", "parse_pdb", "parse_gltf", "rigid"]
