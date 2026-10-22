"""Ce que vaut le mode variant structural, mesuré.

Critères fixés **avant** de regarder les allèles mutants (seuils dans `CRITERIA`) :

1. **Le type sauvage isole** : l'enhancer touche le promoteur d'EPHA4 au moins 3 fois plus que
   celui de PAX3 (3D).
2. **La frontière décide** : sur le DEL-B central, le contact enhancer → promoteur de PAX3 est au
   moins 1,5 fois celui des deux témoins — DelBS (frontière gardée, comme chez la souris) et
   « frontière épargnée » (mêmes points de cassure, seuls les sites CTCF restent). En 3D, et en
   1D pour **chaque** DEL-B de la famille autorisée par le texte.
3. **Les autres frontières tiennent** : sur le DEL-B, le contact enhancer → IHH reste dans
   [0,67 ; 1,5] fois celui du type sauvage (« all novel interactions in the fused TADs respected
   the adjacent boundaries »).
4. **La prédiction 1D suit le 3D** sur les allèles qu'elle n'a pas vus : Pearson des log du 4C
   virtuel depuis PAX3 ≥ 0,9, et le rapport du critère 2 du même côté de 1,5.

Ces critères sont jugés tels quels sur le plan pré-enregistré (`check`). Le plan révisé
(`revised` : fenêtre constante, loi 1D sans plancher ; raisons dans `REVISION`) les rejuge à
seuils inchangés ; le rapport garde les deux.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..fine.extrusion import simulate as extrude
from . import locus as L
from . import predict as P
from .run import alleles, load
from .run import path as run_path

CRITERIA = {"wt_ratio": 3.0, "boundary": 1.5, "ihh_lo": 0.67, "ihh_hi": 1.5, "r_log": 0.9}
PARAMS = {"separation": 200_000, "processivity": 200_000, "velocity": 1.0, "release": 0.003,
          "burnin": 2_000, "snapshots": 150, "stride": 4, "bpPerBead": 2_000}


def _ex(region, sites, seed: int, snapshots: int = PARAMS["snapshots"]):
    return extrude(region, sites, separation=PARAMS["separation"], processivity=PARAMS["processivity"],
                   velocity=PARAMS["velocity"], release=PARAMS["release"], burnin=PARAMS["burnin"],
                   snapshots=snapshots, stride=PARAMS["stride"], seed=seed)


def landmarks(loc: L.Locus, mids: np.ndarray) -> dict[str, np.ndarray]:
    out = {"enhancer": loc.monomers(mids, *loc.enhancer)}
    for g in ("EPHA4", "PAX3", "IHH", "WNT6"):
        out[g] = loc.promoter(g, mids)
    return out


def contacts(rows: np.ndarray, src: np.ndarray, lm: dict[str, np.ndarray]) -> dict[str, float]:
    out = {}
    for g in ("EPHA4", "PAX3", "IHH", "WNT6"):
        out[g] = P.contact(rows, src, lm["enhancer"], lm[g]) if len(lm[g]) else float("nan")
    return out


def predicted(loc: L.Locus, v: L.Variant | None, law: P.Fit, seeds: range, sources: np.ndarray | None = None):
    sites, mids = loc.sites(v)
    region = loc.allele_region(v)
    lm = landmarks(loc, mids)
    src = np.unique(np.r_[lm["enhancer"], lm["PAX3"]] if sources is None else sources)
    legs = np.concatenate([np.stack(_ex(region, sites, s).legs) for s in seeds])
    return P.rows(region.n, legs, src, law.alpha, law.d0, law.c), src, lm, mids


def measure(loc: L.Locus, runs: dict, fit: P.Fit, log=print) -> tuple[dict, dict]:
    """Contacts 3D et 1D de chaque allèle, accord du 4C, et les rapports des critères 2 et 3."""
    out: dict = {}
    for t, (coords, legs, meta) in runs.items():
        mids = np.asarray(meta["mids_hg38"])
        lm = landmarks(loc, mids)
        src = np.unique(np.r_[lm["enhancer"], lm["PAX3"]])
        obs = P.observed(coords, src)
        pre = P.rows(coords.shape[1], legs, src, fit.alpha, fit.d0, fit.c)
        k = {int(s): i for i, s in enumerate(src)}
        # Le 4C depuis PAX3, en casiers de 20 kb (à 2 kb, la plupart des monomères n'ont aucun
        # contact sur 1 200 instantanés) ; le point de vue lui-même exclu.
        r4c = obs[[k[int(i)] for i in lm["PAX3"]]].mean(axis=0)
        p4c = pre[[k[int(i)] for i in lm["PAX3"]]].mean(axis=0)
        far = np.abs(np.arange(len(r4c)) - lm["PAX3"].mean()) >= 10
        r4c, p4c = (P.binned(np.where(far, x, np.nan)[None])[0] for x in (r4c, p4c))
        ok = np.isfinite(r4c) & (r4c > 0)
        # Le plafond : le 3D contre lui-même, réplicats 1–4 contre 5–8 (le bruit du 3D borne
        # l'accord qu'aucune prédiction ne peut dépasser).
        half = len(coords) // 2
        h = [P.binned(np.where(far, P.observed(c, lm["PAX3"]).mean(axis=0), np.nan)[None])[0]
             for c in (coords[:half], coords[half:])]
        okh = np.isfinite(h[0]) & (h[0] > 0) & (h[1] > 0)
        split = float(np.corrcoef(np.log(h[0][okh]), np.log(h[1][okh]))[0, 1])
        out[t] = {
            "name": (meta.get("variant") or {}).get("name", "type sauvage"),
            "monomers": int(coords.shape[1]), "frames": int(coords.shape[0]),
            "3d": contacts(obs, src, lm), "1d": contacts(pre, src, lm),
            "r_log_4c_pax3": float(np.corrcoef(np.log(p4c[ok]), np.log(r4c[ok]))[0, 1]),
            "r_log_4c_pax3_split_half": split,
        }
        a = out[t]
        log(f"  {t:10s} 3D  enh→EPHA4 {a['3d']['EPHA4']:.4f}  enh→PAX3 {a['3d']['PAX3']:.4f}  enh→IHH {a['3d']['IHH']:.4f}"
            f"   | 1D enh→PAX3 {a['1d']['PAX3']:.4f}  · r(log, 4C PAX3, 20 kb) {a['r_log_4c_pax3']:.3f}"
            f" (3D contre 3D, moitiés : {a['r_log_4c_pax3_split_half']:.3f})")
    verdicts = {}
    for dim in ("3d", "1d"):
        verdicts[f"verdict_{dim}"] = {
            "wt_ratio": out["wt"][dim]["EPHA4"] / out["wt"][dim]["PAX3"],
            "vs_delbs": out["delb"][dim]["PAX3"] / out["delbs"][dim]["PAX3"],
            "vs_keep": out["delb"][dim]["PAX3"] / out["delb-keep"][dim]["PAX3"],
            "vs_wt": out["delb"][dim]["PAX3"] / out["wt"][dim]["PAX3"],
            "ihh": out["delb"][dim]["IHH"] / out["wt"][dim]["IHH"],
        }
    return out, verdicts


REVISION = {
    "window": (
        "Chaque allèle garde 5 Mb (2 500 monomères), la délétion faisant entrer la séquence en aval. "
        "Raison : le confinement est réglé sur la longueur de la chaîne, et une chaîne plus courte a un "
        "plancher de contacts lointains plus haut (3D, séparations de 500 à 800 monomères : 0,0012 pour "
        "le type sauvage à 2 500 monomères, 0,0023 pour DEL-B à 1 600) ; DEL-B (1 600) et DelBS (1 815) "
        "n'étaient donc pas comparés à confinement égal."),
    "law": (
        "Loi 1D sans plancher (c = 0), calée sur les paires à moins de 1 Mb du 3D du type sauvage. "
        "Raison : le plancher c = 0,0008 calé sur toutes les paires reproduit le plateau que la sphère "
        "de confinement impose au-delà d'un mégabase — l'artefact déjà relevé à la semaine 7 — et il "
        "dépasse les contacts lointains eux-mêmes : en 1D, le type sauvage n'isolait plus (1,17)."),
    "thresholds": "Seuils inchangés.",
}


def _law(path: Path) -> P.Fit:
    return P.Fit(**json.loads(path.read_text()))


def family(loc: L.Locus, fit: P.Fit, seeds: int, log=print) -> tuple[list[dict], dict]:
    """Toute la famille, en 1D (graines indépendantes de celles du 3D)."""
    log(f"famille : {sum(v.kind == 'DEL-B' for v in loc.variants)} DEL-B, {seeds} exécutions 1D par allèle …")
    rng = range(91_000, 91_000 + seeds)
    cache: dict = {}

    def c(v):
        key = (v.kind, v.start, v.end)
        if key not in cache:
            rows, src, lm, _ = predicted(loc, v, fit, rng)
            cache[key] = contacts(rows, src, lm)["PAX3"]
        return cache[key]

    fam = []
    for v in loc.variants:
        if v.kind != "DEL-B":
            continue
        s = next(x for x in loc.variants if x.kind == "DelBS" and x.start == v.start)
        k = next(x for x in loc.variants if x.kind == "DEL-B, frontière épargnée" and x.start == v.start and x.end == v.end)
        fam.append({"name": v.name, "start": v.start, "end": v.end, "delb": c(v), "delbs": c(s), "keep": c(k)})
    rb = [f["delb"] / f["delbs"] for f in fam]
    rk = [f["delb"] / f["keep"] for f in fam]
    log(f"  min DEL-B / DelBS = {min(rb):.2f}, min DEL-B / frontière épargnée = {min(rk):.2f} sur {len(fam)}")
    return fam, {"vs_delbs": float(min(rb)), "vs_keep": float(min(rk))}


def verdict(res: dict) -> dict:
    C, A = CRITERIA, res["alleles"]
    v3, v1, fm = res["verdict_3d"], res["verdict_1d"], res["family_min"]
    return {
        "1_wt_isolates": v3["wt_ratio"] >= C["wt_ratio"],
        "2_boundary_3d": v3["vs_delbs"] >= C["boundary"] and v3["vs_keep"] >= C["boundary"],
        "2_boundary_family_1d": fm["vs_delbs"] >= C["boundary"] and fm["vs_keep"] >= C["boundary"],
        "3_other_boundaries": C["ihh_lo"] <= v3["ihh"] <= C["ihh_hi"],
        "4_predictor": all(A[t]["r_log_4c_pax3"] >= C["r_log"] for t in ("delb", "delbs", "delb-keep"))
                       and v1["vs_delbs"] >= C["boundary"] and v1["vs_keep"] >= C["boundary"],
    }


def revised(sv_dir: Path, *, family_seeds: int = 8, log=print) -> dict:
    """Les quatre critères, seuils inchangés, sur le plan à fenêtre constante
    (`epha4-<allèle>.w.npz`) et, si elle est calée (`law-nofloor.json`), la loi sans plancher."""
    loc = L.build()
    nofloor = sv_dir / "law-nofloor.json"
    fit = _law(nofloor if nofloor.exists() else sv_dir / "law.json")
    why = [REVISION["window"]] + ([REVISION["law"]] if nofloor.exists() else []) + [REVISION["thresholds"]]
    runs = {t: load(run_path(sv_dir, t, window=True)) for t in alleles(loc)}
    log(f"plan révisé : fenêtre constante de 5 Mb ; loi α = {fit.alpha}, d0 = {fit.d0}, c = {fit.c:.3g}")
    res: dict = {"why": why, "law": fit.__dict__, "criteria": CRITERIA}
    res["alleles"], verdicts = measure(loc, runs, fit, log)
    res.update(verdicts)
    res["family"], res["family_min"] = family(loc, fit, family_seeds, log)
    res["pass"] = verdict(res)
    return res


def check(sv_dir: Path, *, family_seeds: int = 8, log=print) -> dict:
    loc = L.build()
    tags = alleles(loc)
    runs = {t: load(sv_dir / f"epha4-{t}.npz") for t in tags}
    coords, legs, _ = runs["wt"]
    cache = sv_dir / "law.json"
    if cache.exists():
        fit = _law(cache)
        log("loi relue (law.json, calée sur le 3D du type sauvage) :")
    else:
        log("calage de la loi sur le 3D du type sauvage …")
        fit = P.calibrate(coords, legs)
        cache.write_text(json.dumps(fit.__dict__))
    log(f"  α = {fit.alpha}, d0 = {fit.d0} monomère(s), plancher c = {fit.c:.3g}, r(log, casiers de 20 kb) = {fit.r_log:.3f}")

    res: dict = {"law": fit.__dict__, "criteria": CRITERIA}
    res["alleles"], verdicts = measure(loc, runs, fit, log)
    res.update(verdicts)
    res["family"], res["family_min"] = family(loc, fit, family_seeds, log)
    res["pass"] = verdict(res)
    return res


BATCHES = 16


def dump(out: Path, seeds: int = 40) -> None:
    """Pour la lecture croisée Python → TypeScript (`tools/test-variant.mjs`) : les sites de
    l'allèle sauvage, les paramètres, et ce que l'extrusion Python en tire sur `seeds` graines —
    densité de pieds, et lignes de contact depuis l'enhancer et PAX3 sous une loi fixe."""
    loc = L.build()
    sites, mids = loc.sites(None)
    region = loc.allele_region(None)
    lm = landmarks(loc, mids)
    src = np.unique(np.r_[lm["enhancer"][::5], lm["PAX3"]])
    def batch(first: int):
        exs = [_ex(region, sites, first + s) for s in range(seeds)]
        legs = np.concatenate([np.stack(e.legs) for e in exs])
        return exs, np.mean([e.foot_density() for e in exs], axis=0), P.rows(region.n, legs, src, 1.5, 1.0)

    exs, foot, rows = batch(70_000)
    # Un second lot Python, graines disjointes : le plancher de bruit contre lequel l'écart
    # TypeScript–Python se juge.
    _, foot_b, rows_b = batch(80_000)
    # Les contacts enhancer → promoteur sur `BATCHES` lots indépendants : un lot seul est trop
    # bruité (les instantanés d'une exécution sont corrélés, les boucles longues rares) pour
    # juger un écart de quelques pour cent ; on compare des distributions (Welch).
    enh = np.array([i for i in lm["enhancer"] if i in set(src.tolist())])
    batches: dict[str, list[float]] = {g: [] for g in ("EPHA4", "PAX3", "IHH")}
    for b in range(BATCHES):
        lg = np.concatenate([np.stack(_ex(region, sites, 30_000 + 1_000 * b + s).legs) for s in range(seeds)])
        R = P.rows(region.n, lg, src, 1.5, 1.0)
        for g in batches:
            batches[g].append(P.contact(R, src, enh, lm[g]))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "batches": batches,
        "n": region.n, "params": PARAMS, "seeds": seeds, "law": {"alpha": 1.5, "d0": 1.0},
        "sites": {"pos": sites.pos.tolist(), "strand": sites.strand.tolist(), "occupancy": sites.occupancy.tolist()},
        "sources": src.tolist(), "landmarks": {k: v.tolist() for k, v in lm.items()},
        "foot": foot.tolist(), "rows": rows.tolist(), "foot_b": foot_b.tolist(), "rows_b": rows_b.tolist(),
        "occupancy_lef": float(np.mean([e.occupancy for e in exs])),
    }))


def clinical_beads(fiche: Path, out: Path, levels: tuple[str, ...] = ("apercu", "noyau")) -> None:
    """Pour la lecture croisée de la surimpression (`tools/overlay.mjs`) : par bille, ce que le
    lecteur Python lit dans le fichier de fiches — P/LP ClinVar, associations GWAS, LOEUF
    minimal des gènes de la bille."""
    from ..fiche.build import Fiche

    f = Fiche(fiche)
    res: dict = {}
    for lv in levels:
        n = f.level(lv)["n_beads"]
        cv, gw, lo = [], [], []
        for k in range(1, n + 1):
            c = f.card(lv, k)
            cv.append(c["clinvar"].get("v", {}).get("plp_total"))
            gw.append(c["gwas"].get("v", {}).get("n"))
            vals = [x.get("loeuf") for x in (c["contrainte"].get("v") or {}).values() if x.get("loeuf") is not None]
            lo.append(min(vals) if vals else None)
        res[lv] = {"clinvar": cv, "gwas": gw, "loeuf": lo}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res))
