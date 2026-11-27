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
- **fin-tp53** (2 kb, chr17:5,67–9,67 Mb) — la même méthode autour de TP53 (semaine 12).

Et, à la fin du fichier, la **séquence de référence** des régions fines (semaine 12), lue
dans le `.2bit` verrouillé du manifeste, avec les **cibles** de la plongée situées par
l'annotation GENCODE plutôt que recopiées.
"""

from __future__ import annotations

import hashlib
import re
import subprocess
from pathlib import Path

import numpy as np

from ..ensemble import medoid
from ..ensemble import read as read_ensemble
from ..ensemble.stats import kabsch
from ..fine.run import load as load_fine
from .g3d import Attribute, Level, Links, Sequence


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


def compartments(ens: dict, track_path: str | Path | None = None):
    """L'attribut `compartment` du niveau noyau, et le bilan qui dit ce qu'il vaut.

    Toujours calculé depuis l'ensemble (`ensemble.compartments`). Si une piste E1
    mesurée est fournie — un bedGraph, Rao 2014 par exemple —, c'est **elle** qui est
    affichée, et le bilan dit à quel point le modèle s'en approche. Aucune ne l'a été :
    le réseau de l'environnement est fermé, comme pour le DamID et le CTCF.
    """
    from ..ensemble.compartments import compute
    from ..ensemble.stats import overlap_mean, read_bedgraph

    b = ens["beads"]
    model = compute(ens["coords_um"], b.radius, b.copy_id, tuple(b.labels), b.start, b.end,
                    b.lad, b.nuclear_radius)
    meta = {
        "kind": "E1",
        "meaning": "premier vecteur propre cis ; > 0 : compartiment A, < 0 : B",
        "evidence": "simulated",
        "source": f"carte de contacts de l'ensemble ({model.n_structures} structures, "
                  f"seuil {model.cutoff:g}), appelants de la semaine 3",
        "phasing": "piste LAD : A du côté pauvre en LAD (pas de GC lisible ici)",
        "strength_cross_validated": round(model.strength, 3),
        "split_half_r": round(model.split_half_r, 3),
        "absent": "X inactif (pas de compartiments A/B in vivo) et bins non appelés",
    }
    values = model.e1
    agreement = None
    if track_path:
        track, _ = read_bedgraph(track_path)
        measured = overlap_mean(b.start, b.end, b.copy_id, tuple(b.labels), track)
        measured[np.array([b.labels[c].endswith(":Xi") for c in b.copy_id])] = np.nan
        ok = np.isfinite(measured) & np.isfinite(model.e1)
        agreement = {
            "r": float(np.corrcoef(measured[ok], model.e1[ok])[0, 1]) if ok.sum() > 2 else float("nan"),
            "same_sign": float(((measured[ok] > 0) == (model.e1[ok] > 0)).mean()) if ok.any() else float("nan"),
            "covered": int(ok.sum()),
        }
        values = measured
        meta = {**meta, "evidence": "measured", "source": str(track_path),
                "model_agreement": {k: round(v, 3) if isinstance(v, float) else v
                                    for k, v in agreement.items()}}
    return Attribute(values, meta), model, agreement


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
    # Les attributs suivent : moyenne des billes filles qui en ont une, absente
    # seulement si aucune n'en a. Une bille de 3 Mb à moitié dans l'Xi n'existe pas —
    # les copies ne se mélangent jamais — donc « absente » reste une propriété de copie.
    attributes = {}
    for key, a in level.attributes.items():
        v = np.asarray(a.values, dtype=np.float64)
        ok = np.isfinite(v)
        total = np.add.reduceat(np.where(ok, v, 0.0), first_of)
        n_ok = np.add.reduceat(ok.astype(np.float64), first_of)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.where(n_ok > 0, total / n_ok, np.nan)
        attributes[key] = Attribute(mean, {**a.meta, "coarsened": f"moyenne de {k} billes"})
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
        attributes=attributes,
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


LINK_KINDS = ["cohesin_ctcf", "cohesin_one_side", "cohesin_free", "ctcf_pair"]


def anchor_frequency(legs_: np.ndarray, pairs: np.ndarray, tol: int = 2) -> np.ndarray:
    """Fraction des conformations où une cohésine tient la paire `(i, j)`, à ±`tol` monomères.

    Même définition que `Extrusion.anchor_frequency` (semaine 8), sur le tableau de jambes
    de toutes les conformations au lieu d'une seule trajectoire.
    """
    lo = np.where(legs_[..., 0] >= 0, np.minimum(legs_[..., 0], legs_[..., 1]), -10**9)
    hi = np.where(legs_[..., 0] >= 0, np.maximum(legs_[..., 0], legs_[..., 1]), -10**9)
    out = np.zeros(len(pairs))
    for p, (i, j) in enumerate(pairs):
        hit = (np.abs(lo - i) <= tol) & (np.abs(hi - j) <= tol)
        out[p] = hit.any(axis=1).mean()
    return out


def loops(path: str | Path, frame: int, *, max_span_bp: int = 1_000_000):
    """Sites CTCF par monomère, boucles de la conformation `frame`, paires d'ancres de l'ensemble.

    Les jambes viennent du fichier si elles y sont (exécutions depuis la semaine 11), sinon
    elles sont **régénérées** depuis le 1D déterministe et **vérifiées** contre les
    coordonnées : une paire tenue doit être à ~1 σ, une paire quelconque à même séparation
    non. Sans cet accord, mieux vaut refaire `make fine` que dessiner des boucles fausses.

    Les paires candidates ne viennent **pas** de la vérité plantée : ce sont toutes les
    paires (site `+`, site `−` à sa droite) à moins de `max_span_bp`, ce qu'on ferait sur des
    sites mesurés. La vérité ne sert qu'à juger, après, lesquelles passent le seuil.
    """
    from ..fine.run import check_legs, legs as regenerate

    coords, sites, _, meta = load_fine(str(path))
    z = np.load(str(path))
    if "legs" in z.files:
        L, how = z["legs"].astype(np.int64), "lues dans le fichier"
    else:
        L, how = regenerate(meta, sites, z["seeds"]), "régénérées depuis le 1D"
        held, control = check_legs(coords, L)
        if not (held < 1.6 and control > 2.0 * held):
            raise ValueError(f"jambes régénérées en désaccord avec le 3D ({held:.2f} σ contre {control:.2f})")
    n = int(meta["monomers"])
    bp = int(meta["bp_per_bead"])

    fwd = np.full(n, np.nan)
    rev = np.full(n, np.nan)
    for pos, strand, occ in zip(sites.pos, sites.strand, sites.occupancy):
        target = fwd if strand > 0 else rev
        target[pos] = np.nanmax([target[pos], occ])

    is_f = np.isfinite(fwd)
    is_r = np.isfinite(rev)
    held = L[frame][L[frame][:, 0] >= 0]
    i_c = np.minimum(held[:, 0], held[:, 1])
    j_c = np.maximum(held[:, 0], held[:, 1])
    near = lambda mask, k: mask[max(0, k - 1):k + 2].any()     # pied arrêté sur le site, ±1
    kind_c = np.array([
        0 if near(is_f, a) and near(is_r, b) else 1 if near(is_f, a) or near(is_r, b) else 2
        for a, b in zip(i_c, j_c)
    ], dtype=np.int64)

    plus = np.flatnonzero(is_f)
    minus = np.flatnonzero(is_r)
    pairs = np.array([(a, b) for a in plus for b in minus if 0 < (b - a) * bp <= max_span_bp],
                     dtype=np.int64).reshape(-1, 2)
    freq = anchor_frequency(L, pairs)
    # Les motifs CTCF se groupent (deux sites de même brin sur des monomères voisins, comme
    # en vrai) : un même domaine donnait deux candidates quasi identiques, 1147–1200 et
    # 1148–1200. On garde la plus fréquente de toute grappe dont les deux bouts sont à ±2
    # monomères — la tolérance même de la fréquence d'ancrage.
    keep: list[int] = []
    for k in np.argsort(-freq, kind="stable"):
        if all(abs(pairs[k, 0] - pairs[q, 0]) > 2 or abs(pairs[k, 1] - pairs[q, 1]) > 2 for q in keep):
            keep.append(int(k))
    keep_idx = np.sort(np.asarray(keep, dtype=np.int64))
    pairs, freq = pairs[keep_idx], freq[keep_idx]

    links = Links(
        i=np.r_[i_c, pairs[:, 0]],
        j=np.r_[j_c, pairs[:, 1]],
        kind=np.r_[kind_c, np.full(len(pairs), 3)],
        weight=np.r_[np.full(len(i_c), np.nan), freq],
        kinds=LINK_KINDS,
        meta={
            # Court : l'en-tête fait partie de la première requête (§ S9). Le sens long des
            # sortes est dans ARCHITECTURE.md § 3.
            "weight": f"ctcf_pair : fraction des {len(L)} conformations où une cohésine tient la paire",
            "legs": how,
            "frame": int(frame),
            "candidates": f"paires + / − à moins de {max_span_bp // 1000} kb",
        },
    )
    attrs = {
        "ctcf_forward": Attribute(fwd, {"kind": "occupancy", "evidence": "simulated",
                                        "source": sites.source}),
        "ctcf_reverse": Attribute(rev, {"kind": "occupancy", "evidence": "simulated",
                                        "source": sites.source}),
    }
    return links, attrs


ARC_THRESHOLD = 0.05   # le même seuil que le client (apps/dive/fiber.ts), fixé avant la mesure


def loop_report(path: str | Path, threshold: float = ARC_THRESHOLD) -> dict:
    """Les paires qui passent le seuil, jugées contre la vérité plantée (±2 monomères).

    Les candidates ne doivent rien à la vérité (`loops`) ; c'est ce qui rend le jugement
    possible. Rend, par classe plantée, combien passent le seuil, et combien de paires qui
    passent ne correspondent à aucun domaine planté.
    """
    from ..fine.region import KINDS

    links, _ = loops(path, 0)
    _, _, truth, _ = load_fine(str(path))
    m = links.kind == LINK_KINDS.index("ctcf_pair")
    pi, pj, w = links.i[m].astype(np.int64), links.j[m].astype(np.int64), links.weight[m]
    drawn = w >= threshold
    out = {"threshold": threshold, "candidates": int(m.sum()), "drawn": int(drawn.sum()),
           "max_weight": float(w.max()) if len(w) else 0.0, "by_kind": {}, "unplanted": 0}
    if truth is None:
        return out
    matched = np.zeros(len(pi), dtype=bool)
    for kind_id, name in enumerate(KINDS):
        sel = truth.kind == kind_id
        found = 0
        for l, r in zip(truth.left[sel], truth.right[sel]):
            hit = drawn & (np.abs(pi - l) <= 2) & (np.abs(pj - r) <= 2)
            matched |= hit
            found += int(hit.any())
        out["by_kind"][name] = {"planted": int(sel.sum()), "drawn": found}
    out["unplanted"] = int((drawn & ~matched).sum())
    return out


def fine(path: str | Path, nuc: Level, *, copy_label: str = "chr7:a", name: str = "fin") -> tuple[Level, dict]:
    """Le médoïde de la semaine 8, posé dans le repère du noyau — et l'écart que ça laisse.

    **Pose.** Pour chaque bille du noyau couverte à moitié au moins par la région fine,
    on prend le centroïde des monomères fins qui tombent dans son intervalle. Ces
    ancres sont alignées sur les billes du noyau par Kabsch, rotation et translation
    seulement, réflexion permise. Le RMSD restant est publié avec la pose.

    **Pas d'échelle dans la pose**, et c'est voulu. La semaine 8 a montré que les deux
    modèles ne partagent pas la même `R(s)` — un facteur 2,2 à 2,4 au pire selon l'exécution. Une similitude
    absorberait cet écart dans un facteur d'échelle et le rendrait invisible ; on calcule
    ce facteur, on l'écrit, et on ne l'applique pas.

    **Pourquoi la copie `a`.** Le modèle fin n'a pas d'haplotype : il ne sait rien des deux
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
    links, attrs = loops(path, int(mi))
    level = Level(
        name=name,
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
            "verdict": "le raccord de la semaine 8 ne tient pas (2,2 à 2,4× sur R(s) selon l'exécution, VALIDATION.md § S8 et § S10) ; "
                       "la pose le montre, elle ne le corrige pas",
        },
        attributes=attrs,
        links=links,
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


#: Les régions fines et leur niveau. `actb` garde le nom « fin » des semaines 9 à 11.
FINE = (("actb", "fin"), ("tp53", "fin-tp53"))

#: Les bases que la plongée sait viser, désignées par l'annotation et pas par une coordonnée.
TARGETS = (
    {"name": "TP53 c.818", "transcript": "ENST00000269305", "coding": 818, "codon": 273,
     "why": "codon 273 (CGT, Arg) ; R273H, c.818G>A, est le point chaud le plus fréquent des cancers"},
    {"name": "ACTB c.1", "transcript": "ENST00000646664", "coding": 1, "codon": 1,
     "why": "le A de l'ATG d'ACTB, cible des semaines 10 et 11"},
)


def fine_regions(fine_dir: str | Path, nuc: Level) -> list[tuple[Level, dict, Path]]:
    """Chaque région fine présente, posée sur la copie `a` de son chromosome."""
    out = []
    for stem, name in FINE:
        path = Path(fine_dir) / f"{stem}.npz"
        if not path.exists():
            continue
        chrom = load_fine(str(path))[3]["chrom"]
        level, source = fine(path, nuc, copy_label=f"{chrom}:a", name=name)
        out.append((level, source, path))
    return out


def manifest_entry(manifest: str | Path, ident: str) -> dict:
    for line in Path(manifest).read_text(encoding="utf-8").splitlines():
        col = line.split("\t")
        if col and col[0] == ident and len(col) >= 7:
            return {"id": ident, "status": col[2], "url": col[4], "sha256": col[5], "bytes": int(col[6])}
    raise KeyError(f"{ident} absent du manifeste")


def sequences(levels: list[Level], genome: str | Path, manifest: str | Path) -> list[Sequence]:
    """La séquence de référence sous chaque niveau fin, du `.2bit` verrouillé."""
    from ..sequence import TwoBit

    entry = manifest_entry(manifest, Path(genome).name)
    if entry["status"] != "ok":
        raise ValueError(f"{entry['id']} n'est pas verrouillé dans le manifeste ({entry['status']})")
    out = []
    with TwoBit(genome) as tb:
        for lv in levels:
            chrom = lv.copies[0].split(":")[0]
            start, end = int(lv.start.min()), int(lv.end.max())
            out.append(Sequence("hg38", chrom, start, end, tb.fetch(chrom, start, end),
                                {"file": entry["id"], "url": entry["url"], "sha256": entry["sha256"]}))
    return out


def targets(levels: list[Level], gtf: str | Path, manifest: str | Path,
            genome: str | Path | None = None) -> None:
    """Attache à chaque niveau fin les cibles qu'il couvre, recalculées depuis GENCODE.

    La base est lue dans le génome si on l'a : la cible porte alors la base du brin +,
    et celle du brin codant, de quoi vérifier l'étiquette affichée sans autre source.
    """
    from ..sequence import TwoBit, read_transcript

    tb = TwoBit(genome) if genome and Path(genome).exists() else None
    entry = manifest_entry(manifest, "gencode.annotation.gtf")
    release = re.search(r"gencode\.(v\d+)", entry["url"])
    source = f"GENCODE {release.group(1) if release else '?'} basic, sha256 {entry['sha256'][:12]}…"
    for t in TARGETS:
        tr = read_transcript(gtf, t["transcript"])
        pos = tr.coding_position(t["coding"])
        for lv in levels:
            if lv.copies[0].split(":")[0] == tr.chrom and lv.start.min() < pos <= lv.end.max():
                spec = {"name": t["name"], "gene": tr.gene, "chrom": tr.chrom, "position": pos,
                        "strand": tr.strand, "transcript": tr.transcript, "codon": t["codon"],
                        "coding": f"c.{t['coding']}", "why": t["why"],
                        "tags": [x for x in tr.tags if x in ("MANE_Select", "Ensembl_canonical")],
                        "source": source}
                if tb:
                    plus = tb.fetch(tr.chrom, pos - 1, pos)
                    spec["base_plus"] = plus
                    spec["base_coding"] = plus if tr.strand == "+" else plus.translate(str.maketrans("ACGT", "TGCA"))
                lv.targets.append(spec)
    if tb:
        tb.close()


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
    "encore (cf. DATA_SOURCES.md § 8). La séquence, elle, est la référence hg38.",
    "Une structure n'est pas « le » génome : c'est le médoïde d'un ensemble, et la "
    "variabilité par bille dit de combien un tirage s'en écarte.",
    "Le niveau fin ne se raccorde pas au noyau (VALIDATION.md § S8).",
]
