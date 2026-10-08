"""Construit les niveaux d'un `.g3d` à partir des ensembles des semaines 7 et 8.

Un fichier porte **une** structure — le médoïde d'un ensemble — et non l'ensemble :
deux cents structures pèseraient deux cents fois plus, et la semaine 7 a montré qu'il
n'y a pas de repère commun dans lequel les superposer. Ce qui survit de l'ensemble est
la **variabilité de profondeur** par bille, la seule grandeur invariante par rotation
qu'il mesure. Le fichier dit lequel des deux on regarde, et comment il a été obtenu.

Trois niveaux :

- **aperçu** (3 Mb) — quatre billes du noyau fusionnées en une. Il n'existe que pour le
  premier rendu, et n'entre dans le fichier que si la mesure montre qu'il en a besoin.
- **noyau** (750 kb) — le médoïde de la semaine 7.
- **fin** (2 kb, chr7:4–8 Mb) — le médoïde de la semaine 8, posé dans le repère du
  noyau. La pose ne raccorde pas les deux modèles, elle **mesure** à quel point ils ne
  se raccordent pas : la semaine 8 a établi l'échec, le fichier le transporte.
"""

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path

import numpy as np

from ..ensemble import medoid
from ..ensemble import read as read_ensemble
from ..ensemble.stats import kabsch
from ..fine.run import load as load_fine
from .g3d import Level


def _sha(*arrays: np.ndarray) -> str:
    h = hashlib.sha256()
    for a in arrays:
        h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()


def pipeline_version() -> str:
    """Le commit qui a produit le fichier, et s'il était propre."""
    root = Path(__file__).resolve().parents[3]
    try:
        head = subprocess.run(["git", "-C", str(root), "rev-parse", "--short=12", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(root), "status", "--porcelain"],
                               capture_output=True, text=True, check=True).stdout.strip()
        return head + ("+modifié" if dirty else "")
    except (OSError, subprocess.CalledProcessError):
        return "inconnu"


def nucleus(path: str | Path) -> tuple[Level, dict, dict]:
    """Le médoïde de l'ensemble de la semaine 7, et l'écart-type de profondeur par bille.

    Rend aussi les coordonnées de **toutes** les structures, réduites à ce dont l'aperçu
    a besoin : l'aperçu recalcule sa propre variabilité sur l'ensemble au lieu de
    moyenner des écarts-types, ce qui serait faux.
    """
    e = read_ensemble(path)
    b = e.beads
    radial = e.radial()
    mi, _ = medoid(radial)
    scale = b.nuclear_radius * 1000.0
    radius = [float(np.median(b.radius[b.copy_id == c])) * 1000.0
              for c in range(len(b.labels))]
    level = Level(
        name="noyau",
        bp_per_bead=int(e.meta["bp_per_bead"]),
        positions_nm=e.coords[mi] * 1000.0,
        copy=b.copy_id.astype(np.int64),
        start=b.start.astype(np.int64),
        end=b.end.astype(np.int64),
        variability_nm=radial.std(axis=0) * scale,
        copies=list(b.labels),
        radius_nm=radius,
        variability={"kind": "radial_sd", "over": f"{e.n_structures} structures",
                     "reference": "centre du noyau"},
        notes=[
            f"médoïde de {e.n_structures} structures pour la distance entre profils "
            f"radiaux (graine {int(e.seeds[mi])}) — un tirage typique, pas « le » génome",
            "piste LAD synthétique, longueurs de chromosomes intégrées (réseau fermé)",
        ],
    )
    source = {
        "kind": "geno-ensemble/1",
        "path": str(path),
        "structures": e.n_structures,
        "medoid_seed": int(e.seeds[mi]),
        "sha256": _sha(e.coords[mi].astype(np.float32), b.start, b.end, b.copy_id),
    }
    return level, source, {"coords_um": e.coords, "beads": b, "medoid": mi}


def coarsen(level: Level, ens: dict, k: int = 4) -> Level:
    """Fusionne `k` billes consécutives d'une même copie. Le centroïde, pas un choix.

    La variabilité est **recalculée** sur les deux cents structures fusionnées de la
    même façon, pas obtenue en moyennant les écarts-types des billes filles : la
    moyenne d'écarts-types n'est l'écart-type de rien.
    """
    b = ens["beads"]
    copy = b.copy_id
    # Les billes d'une copie sont contiguës et triées : un groupe commence à chaque
    # changement de copie, puis toutes les `k` billes à l'intérieur d'une copie.
    copy_start = np.r_[0, np.flatnonzero(np.diff(copy)) + 1]
    first_of = np.concatenate([
        np.arange(s, e, k) for s, e in zip(copy_start, np.r_[copy_start[1:], len(copy)])
    ])
    last_of = np.r_[first_of[1:], len(copy)] - 1
    counts = (last_of - first_of + 1).astype(np.float64)

    all_merged = np.add.reduceat(ens["coords_um"], first_of, axis=1) / counts[:, None]
    radial = np.linalg.norm(all_merged, axis=2) / b.nuclear_radius
    gcopy = copy[first_of]
    radius = []
    for c in range(len(b.labels)):
        m = gcopy == c
        # r ∝ L^(1/3) : la loi de la semaine 6, appliquée à la longueur fusionnée.
        radius.append(level.radius_nm[c] * float(np.cbrt(np.median(counts[m]))) if m.any()
                      else level.radius_nm[c])
    return Level(
        name="apercu",
        bp_per_bead=level.bp_per_bead * k,
        positions_nm=all_merged[ens["medoid"]] * 1000.0,
        copy=gcopy.astype(np.int64),
        start=b.start[first_of].astype(np.int64),
        end=b.end[last_of].astype(np.int64),
        variability_nm=radial.std(axis=0) * b.nuclear_radius * 1000.0,
        copies=list(level.copies),
        radius_nm=radius,
        variability=dict(level.variability),
        notes=[f"{k} billes du noyau fusionnées par centroïde ; variabilité recalculée "
               f"sur l'ensemble fusionné"],
    )


def pose(target: np.ndarray, local: np.ndarray) -> tuple[np.ndarray, dict]:
    """Transformation rigide qui envoie les ancres `local` sur les ancres `target`.

    Rend une matrice 4×4 **en ligne** (`p_cible = M · [p; 1]`) et le bilan de l'ajustement.
    Le facteur d'échelle optimal (Umeyama) est calculé et publié, **pas appliqué** : il
    sert à dire si l'écart restant est une affaire de taille ou de forme.
    """
    rot, rmsd = kabsch(target, local)
    lc = local - local.mean(axis=0)
    tc = target - target.mean(axis=0)
    scale = float(np.trace(tc.T @ (lc @ rot)) / (lc**2).sum())
    scaled = float(np.sqrt(((scale * (lc @ rot) - tc) ** 2).sum(axis=1).mean()))
    m = np.eye(4)
    m[:3, :3] = rot.T                       # (p − l̄) @ R  ==  Rᵀ (p − l̄) en colonne
    m[:3, 3] = target.mean(axis=0) - rot.T @ local.mean(axis=0)
    return m, {
        "method": "kabsch, réflexion permise, sans échelle",
        "rmsd_nm": round(rmsd, 1),
        "scale_if_allowed": round(scale, 3),
        "rmsd_if_scaled_nm": round(scaled, 1),
    }


def fine(path: str | Path, nuc: Level, *, copy_label: str = "chr7:a") -> tuple[Level, dict]:
    """Le médoïde de la semaine 8, posé dans le repère du noyau — et l'écart que ça laisse.

    **Pose.** Pour chaque bille du noyau couverte à moitié au moins par la région fine,
    on prend le centroïde des monomères fins qui tombent dans son intervalle. Ces
    ancres sont alignées sur les billes du noyau par Kabsch, rotation et translation
    seulement, réflexion permise. Le RMSD restant est publié avec la pose.

    **Pas d'échelle dans la pose**, et c'est voulu. La semaine 8 a montré que les deux
    modèles ne partagent pas la même `R(s)` — un facteur 2,42 au pire. Une similitude
    absorberait cet écart dans un facteur d'échelle et le rendrait invisible ; on calcule
    ce facteur, on l'écrit, et on ne l'applique pas.

    **Pourquoi chr7:a.** Le modèle fin n'a pas d'haplotype : il ne sait rien des deux
    homologues. Le poser sur `a` plutôt que sur `b` est arbitraire, et le fichier le dit.
    """
    coords, sites, _, meta = load_fine(str(path))
    sigma = float(meta["sigma_nm"])
    start0, bp = int(meta["start"]), int(meta["bp_per_bead"])
    n = int(meta["monomers"])

    radial = np.linalg.norm(coords, axis=2).astype(np.float64)   # centre du confinement
    mi, _ = medoid(radial)
    local = coords[mi].astype(np.float64) * sigma
    mstart = start0 + np.arange(n) * bp
    mmid = mstart + bp // 2

    c = nuc.copies.index(copy_label)
    on = np.flatnonzero(nuc.copy == c)
    anchors_fine, anchors_nuc, used = [], [], []
    for j in on:
        s, e = int(nuc.start[j]), int(nuc.end[j])
        inside = (mmid >= s) & (mmid < e)
        if inside.sum() * bp >= 0.5 * (e - s):
            anchors_fine.append(local[inside].mean(axis=0))
            anchors_nuc.append(nuc.positions_nm[j])
            used.append(f"{copy_label}:{s:,}-{e:,}")
    a = np.asarray(anchors_nuc)
    f = np.asarray(anchors_fine)
    if len(a) < 3:
        raise ValueError(f"{len(a)} ancre(s) seulement : une pose 3D en demande au moins trois")

    m, fit = pose(a, f)
    level = Level(
        name="fin",
        bp_per_bead=bp,
        positions_nm=local,
        copy=np.zeros(n, dtype=np.int64),
        start=mstart.astype(np.int64),
        end=(mstart + bp).astype(np.int64),
        variability_nm=radial.std(axis=0) * sigma,
        copies=[copy_label],
        radius_nm=[sigma / 2.0],
        variability={"kind": "radial_sd", "over": f"{len(coords)} conformations",
                     "reference": "centre du confinement de la région"},
        transform=m.reshape(-1).tolist(),
        fit={
            **fit,
            "anchors": used,
            "verdict": "le raccord de la semaine 8 ne tient pas (2,42× sur R(s)) ; "
                       "la pose le montre, elle ne le corrige pas",
        },
        notes=[
            f"médoïde de {len(coords)} conformations, extrusion de boucles sous OpenMM",
            "sites CTCF plantés (vérité connue de la semaine 8), pas mesurés",
            f"posé sur {copy_label} arbitrairement : le modèle fin n'a pas d'haplotype",
        ],
    )
    source = {
        "kind": meta.get("format", "geno-fine/1"),
        "path": str(path),
        "conformations": len(coords),
        "ctcf": sites.source,
        "sha256": _sha(coords[mi]),
    }
    return level, source


def header(sources: list[dict], *, warnings: list[str]) -> dict:
    return {
        "assembly": "GRCh38",
        "cell_type": "GM12878",
        "karyotype": "46,XX",
        "provenance": {"pipeline": pipeline_version(), "sources": sources},
        "warnings": warnings,
    }


WARNINGS = [
    "Toutes les positions sont simulées : aucune donnée de conformation ne les contraint "
    "(réseau fermé, cf. DATA_SOURCES.md § 8).",
    "Une structure n'est pas « le » génome : c'est le médoïde d'un ensemble, et la "
    "variabilité par bille dit de combien un tirage s'en écarte.",
    "Le niveau fin ne se raccorde pas au noyau (VALIDATION.md § S8).",
]
