"""Le rapport de la semaine 21 : les 13 critères, les témoins, les figures.

    python -m geno_pipeline validate          (make validate)

Écrit `data/validation/s21.json` et `docs/figures/s21/*.png`. Les permutations ont une graine
fixe : deux exécutions donnent le même fichier.
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

from ..fiche.sources import locked, manifest
from ..liftover import Chain
from . import lamina, sprite, tracing
from .common import DATA, ROOT, Model, halves, spearman

OUT = DATA / "validation" / "s21.json"
FIGS = ROOT.parent / "docs" / "figures" / "s21"


def _clean(o):
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else round(float(o), 6)
    if isinstance(o, np.integer):
        return int(o)
    return o


def positive_controls(model: Model, chain: Chain, genes) -> dict:
    """Une moitié des repliements contre l'autre, mesurée comme des données (seuil 0,9)."""
    A, B = halves(model)
    radial = spearman(A.radial(), B.radial())
    ca, cb = sprite.load_intra(A, chain, genes), sprite.load_intra(B, chain, genes)
    oa = np.concatenate([sprite.oe_vec(1.0 / c.med[c.i, c.j], c.k) for c in ca])
    ob = np.concatenate([sprite.oe_vec(1.0 / c.med[c.i, c.j], c.k) for c in cb])
    intra = spearman(oa, ob)
    return {"radial": radial, "intra_oe": intra, "ok": bool(radial >= 0.9 and intra >= 0.9)}


def figures(results, extras) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGS.mkdir(parents=True, exist_ok=True)
    out = []
    by = {r.id: r for r in results}

    def save(fig, name):
        p = FIGS / name
        fig.tight_layout(); fig.savefig(p, dpi=110); plt.close(fig); out.append(str(p))

    # P(s) — SPRITE et modèle, normalisées à 1 Mb.
    e = by["C1d"].extra
    fig, ax = plt.subplots(figsize=(5, 4))
    for key, label in (("curve_sprite", "SPRITE GM12878"), ("curve_model", "modèle (contacts)")):
        x, y = np.array(e[key][0][1:]), np.array(e[key][1][1:])
        ax.loglog(x, y / y[0], "o-", ms=3, label=label)
    ax.axvspan(1e6, 1e7, color="0.9"); ax.set_xlabel("séparation (pb)"); ax.set_ylabel("P(s) / P(1 Mb)")
    ax.set_title(f"C1d : pentes {e['sprite']:.2f} (SPRITE) et {e['model']:.2f} (modèle)"); ax.legend()
    save(fig, "ps.png")

    # R(s) — traçage chr21 et modèle.
    e = by["C3a"].extra
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.loglog(e["curve_tracing"][0], e["curve_tracing"][1], "o-", ms=3, label="traçage chr21 (IMR90)")
    ax.loglog(e["curve_model"][0], e["curve_model"][1], "s-", ms=3, label="modèle")
    ax.set_xlabel("séparation (pb)"); ax.set_ylabel("distance médiane (µm)")
    ax.set_title(f"C3a : ν = {e['tracing']:.2f} (traçage), {e['model']:.2f} (modèle)"); ax.legend()
    save(fig, "rs.png")

    # Positions radiales — traçage et modèle, par locus.
    g = extras["genome"]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(g["radial_tracing"], g["radial_model"], s=4, alpha=0.5)
    ax.set_xlabel("position radiale, traçage (IMR90)"); ax.set_ylabel("position radiale, modèle")
    ax.set_title(f"C4a : ρ = {by['C4a'].value:.2f}")
    save(fig, "radial.png")

    # Lamina — DamID et modèle, par bille.
    l = extras["lamina"]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(l["damid"], l["radial"], s=2, alpha=0.3)
    ax.set_xlabel("pA-DamID LMNB1, K562 (log2)"); ax.set_ylabel("position radiale du modèle")
    ax.set_title(f"C5 : ρ = {by['C5'].value:.2f}")
    save(fig, "damid.png")

    # Le tableau de bord : chaque ρ contre son seuil.
    ids = [r.id for r in results if r.id not in ("C1c", "C1d", "C3a", "C3c", "C6")]
    vals = [by[i].value for i in ids]
    colors = {"accord": "#3a8f5c", "partiel": "#d39b2a", "échec": "#c0453a"}
    fig, ax = plt.subplots(figsize=(6, 3.6))
    ax.bar(ids, vals, color=[colors[by[i].verdict] for i in ids])
    ax.axhline(0.3, color="0.4", lw=0.8, ls="--"); ax.axhline(-0.3, color="0.4", lw=0.8, ls="--"); ax.axhline(0, color="0.2", lw=0.6)
    ax.set_ylabel("ρ de Spearman"); ax.set_title("Critères en corrélation (vert accord, orange partiel, rouge échec)")
    save(fig, "summary.png")
    return out


def run() -> int:
    t0 = time.time()
    step = lambda m: print(f"  [{time.time() - t0:5.0f} s] {m}", flush=True)
    model = Model.load()
    step(f"modèle : {model.coords.shape[0]} repliements, {model.coords.shape[1]} billes, piste LAD {model.lad_source}")
    m = manifest()
    chain = Chain(locked(m, "chain.hg19.to.hg38").path)
    genes = sprite.gene_starts()
    step("gènes GENCODE lus")
    chroms = sprite.load_intra(model, chain, genes)
    sizes = []
    import tarfile, io, pandas as pd
    with tarfile.open(locked(m, "sprite.gm12878.intra1mb").path) as tar:
        for ch in sprite.CHROMS:
            sizes.append(sum(1 for _ in io.TextIOWrapper(tar.extractfile(f"human_{ch}_1Mb.txt"))))
    inter = sprite.load_inter(model, sizes)
    step("SPRITE intra et inter lus")
    results = sprite.criteria(model, chroms, inter)
    step("C1, C2")
    c21 = tracing.load_chr21()
    results += tracing.chr21_criteria(model, c21)
    step(f"C3 ({c21.copies} copies de chr21)")
    gen = tracing.load_genome()
    r4, fig_genome = tracing.genome_criteria(model, gen)
    results += r4
    step(f"C4 ({gen.cells} cellules)")
    r5, fig_lamina = lamina.criteria(model)
    results += r5
    step("C5, C6")
    controls = positive_controls(model, chain, genes)
    step("témoin positif")

    lost = {c.name: c.lost for c in chroms}
    payload = _clean({
        "criteria": [asdict(r) for r in results],
        "positive_control": controls,
        "model": {"structures": int(model.coords.shape[0]), "beads": int(model.coords.shape[1]),
                  "lad_source": model.lad_source},
        "sprite": {"bins_lost_in_liftover": lost, "same_bead_pairs": sum(c.same_bead for c in chroms)},
        "valid": controls["ok"],
    })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=1))
    figs = figures(results, {"genome": fig_genome, "lamina": fig_lamina})
    step(f"écrit {OUT} et {len(figs)} figures")

    print("\n  Validation — semaine 21 (critères : docs/S21_CRITERES.md)")
    for r in results:
        ref = f" · N1 {r.reference:+.3f}" if r.reference is not None and r.id == "C1a" else ""
        null = f" · permutations p95 {r.null_p95:+.3f}" if r.null_p95 is not None else ""
        print(f"  {r.id:<4} {r.verdict:<8} {r.value:+.3f}  ({r.threshold}){ref}{null}  — {r.measure} · {r.data}")
    print(f"  témoin positif : radial {controls['radial']:.3f}, intra o/e {controls['intra_oe']:.3f} "
          f"{'✓' if controls['ok'] else '✗ rapport invalide'}")
    counts = {v: sum(r.verdict == v for r in results) for v in ("accord", "partiel", "échec")}
    print(f"  {counts['accord']} accord · {counts['partiel']} partiel · {counts['échec']} échec sur {len(results)}")
    return 0 if controls["ok"] else 1
