"""CLI `geno` — construire le magasin d'intervalles et l'interroger."""

from __future__ import annotations

import argparse
import json
import random
import statistics
import sys
import time
from pathlib import Path

import numpy as np

from . import locus as locus_mod
from .intervals import Store, Track, write_store
from .locus import RegionError, render, report
from .parsers import read_track

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_STORE = ROOT / "data" / "store"
DEFAULT_TRACKS = ROOT / "fixtures" / "tracks.json"
BENCH_STORE = ROOT / "data" / "bench"


def cmd_build(args: argparse.Namespace) -> int:
    spec_path = Path(args.tracks)
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    base = spec_path.parent
    out = Path(args.out)

    tracks = {t["name"]: read_track(t, base) for t in spec["tracks"]}
    meta = {
        "assembly": spec.get("assembly", "?"),
        "source": spec.get("source", spec_path.name),
        "built_from": str(spec_path),
    }
    index = write_store(out, tracks, meta)

    print(f"magasin écrit en {out}")
    print(f"  assemblage  {index['assembly']}   source  {index['source']}")
    total = 0
    for name, n in index["tracks"].items():
        print(f"  {name:<14} {n:>9,} intervalles")
        total += n
    print(f"  {'total':<14} {total:>9,}")
    return 0


def cmd_query(args: argparse.Namespace) -> int:
    store = Store(Path(args.store))
    try:
        t0 = time.perf_counter_ns()
        rep = report(store, args.region, args.track or None, args.flank)
        elapsed = (time.perf_counter_ns() - t0) / 1e6
    except (RegionError, KeyError) as exc:
        print(f"erreur : {exc}", file=sys.stderr)
        return 2
    finally:
        pass

    if args.json:
        payload = rep.to_dict()
        payload["elapsed_ms"] = round(elapsed, 4)
        print(json.dumps(payload, indent=2, ensure_ascii=False))
    else:
        print(render(rep))
        if args.time:
            print(f"  {elapsed:.3f} ms")
    store.close()
    return 0 if not rep.is_empty() else 1


def cmd_tracks(args: argparse.Namespace) -> int:
    store = Store(Path(args.store))
    print(f"{store.path}   assemblage {store.assembly}   source {store.source}")
    for name, track in store.tracks.items():
        chroms = ", ".join(track.chroms) if len(track.chroms) <= 6 else f"{len(track.chroms)} chromosomes"
        print(f"  {name:<14} {len(track):>9,} intervalles   {chroms}")
    store.close()
    return 0


def _synthetic(n: int, seed: int = 0):
    """Intervalles synthétiques à l'échelle réelle, y compris quelques très longs.

    Les gènes longs sont le cas qui casse un index naïf ; le benchmark serait
    malhonnête sans eux.
    """
    rng = np.random.default_rng(seed)
    sizes = {"chr1": 248_956_422, "chr2": 242_193_529, "chr7": 159_345_973, "chr17": 83_257_441}
    per = n // len(sizes)
    for chrom, size in sizes.items():
        starts = rng.integers(0, size - 3_000_000, per)
        lengths = rng.integers(200, 40_000, per)
        # 0,2 % d'intervalles très longs : DMD fait 2,2 Mb, CNTNAP2 2,3 Mb.
        long_ix = rng.choice(per, max(1, per // 500), replace=False)
        lengths[long_ix] = rng.integers(500_000, 2_400_000, long_ix.size)
        for s, ln in zip(starts.tolist(), lengths.tolist()):
            yield chrom, s, s + ln, {"name": "syn"}


def cmd_bench(args: argparse.Namespace) -> int:
    out = Path(args.store)
    stamp = out / ".n"
    if not (out / "index.json").exists() or not stamp.exists() or stamp.read_text() != str(args.n):
        print(f"construction d'un magasin synthétique de {args.n:,} intervalles…")
        t0 = time.perf_counter()
        write_store(
            out,
            {"syn": _synthetic(args.n)},
            {"assembly": "synthetic", "source": "bench"},
        )
        stamp.write_text(str(args.n))
        print(f"  construit en {time.perf_counter() - t0:.1f} s")

    track = Track(out / "syn")
    rng = random.Random(1)
    chroms = track.chroms
    sizes = {"chr1": 248_956_422, "chr2": 242_193_529, "chr7": 159_345_973, "chr17": 83_257_441}

    # La latence suit le nombre de features ramenées, pas la taille du magasin.
    # Un benchmark qui mélange les tailles de fenêtre masque exactement ça.
    windows = [
        (3_000, "un gène"),
        (30_000, "un gène + flancs"),
        (300_000, "un TAD"),
        (3_000_000, "un compartiment"),
    ]

    print(f"\n  {len(track):,} intervalles · {args.queries} requêtes par ligne\n")
    print(f"  {'fenêtre':<25} {'features':>7}   {'index seul':>22}   {'index + attributs':>22}")
    print(f"  {'-' * 25} {'-' * 7}   {'-' * 22}   {'-' * 22}")

    worst_interactive = 0.0
    for width, label in windows:
        queries = []
        for _ in range(args.queries):
            c = rng.choice(chroms)
            s = rng.randrange(0, sizes[c] - width - 1)
            queries.append((c, s, s + width))

        cells, hits = [], 0
        for is_query, fn in ((False, track.count), (True, track.query)):
            for c, s, e in queries[:50]:  # chauffe
                fn(c, s, e)
            times = []
            for c, s, e in queries:
                t0 = time.perf_counter_ns()
                r = fn(c, s, e)
                times.append((time.perf_counter_ns() - t0) / 1e6)
                hits += r if isinstance(r, int) else len(r)
            times.sort()
            p99 = times[min(len(times) - 1, int(len(times) * 0.99))]
            cells.append(f"méd {statistics.median(times):6.3f}  p99 {p99:6.3f}")
            if is_query and width <= 30_000:
                worst_interactive = max(worst_interactive, p99)

        n = hits // (2 * len(queries))
        print(f"  {width // 1000:>4} kb  {label:<17} {n:>7,}   {cells[0]:>22}   {cells[1]:>22}")

    track.close()
    print(
        f"\n  Cible semaine 2 : < 1 ms sur une requête de locus.\n"
        f"  p99 à l'échelle d'un locus (≤ 30 kb), attributs compris : {worst_interactive:.3f} ms."
    )
    print(
        "  Au-delà, le coût suit le nombre de features : ~3 µs chacune, dominés par le\n"
        "  parsing JSON des attributs. Une requête qui ramène un compartiment entier est un\n"
        "  export, pas une interaction — la fiche de la semaine 13 clique UNE bille."
    )
    return 0


def cmd_hic(args: argparse.Namespace) -> int:
    """Plante une structure connue, fait tourner les callers, mesure l'accord.

    C'est la seule configuration où « le caller est correct » est vérifiable :
    sur des données réelles, un désaccord avec Rao 2014 ne dit pas lequel des
    deux a tort.
    """
    import logging as _logging
    import warnings

    _logging.disable(_logging.INFO)
    # Bruit de cooltools 0.7 sur pandas 2 : `.idxmin()` sur colonne toute-NA.
    # C'est exactement ce qui lève une ValueError sous pandas 3 — d'où l'épinglage.
    warnings.simplefilter("ignore", FutureWarning)
    try:
        import cooler
    except ImportError:
        print(
            "la pile Hi-C n'est pas installée dans cet interpréteur.\n"
            "  → voir docs/SETUP.md ; make hic-validate utilise pipeline/.venv",
            file=sys.stderr,
        )
        return 2

    from .hic import (
        balance,
        compartments,
        gc_track,
        loops,
        match_pairs,
        match_positions,
        plant,
        tad_boundaries,
        write_cool,
    )

    out = Path(args.out)
    print(f"plantation  {args.bins:,} bins × {args.resolution // 1000} kb, graine {args.seed}")
    bins, pixels, truth = plant(n_bins=args.bins, resolution=args.resolution, seed=args.seed)
    write_cool(out, bins, pixels)
    n_contacts = int(pixels["count"].sum())
    print(
        f"            {truth.span / 1e6:.1f} Mb · {n_contacts:,} contacts · "
        f"{len(pixels):,} pixels · {len(truth.boundaries)} TADs · {len(truth.loops)} boucles"
    )
    print(f"            écrit en {out}\n")

    clr = cooler.Cooler(str(out))
    balance(clr)
    clr = cooler.Cooler(str(out))

    cv = lambda m: float(np.nanstd(np.nansum(m, 1)) / np.nanmean(np.nansum(m, 1)))  # noqa: E731
    cv_raw, cv_bal = cv(clr.matrix(balance=False)[:]), cv(clr.matrix(balance=True)[:])
    w = clr.bins()["weight"][:].to_numpy()
    ok = np.isfinite(w)
    r_bias = float(np.corrcoef(np.log(w[ok]), -np.log(truth.bias[ok]))[0, 1])
    print(f"  équilibrage ICE     CV des marginales {cv_raw:.3f} → {cv_bal:.3f}"
          f"   ·  poids vs biais planté r = {r_bias:+.3f}")

    table = compartments(clr, phasing_track=gc_track(truth))
    e1 = table["E1"].to_numpy()
    ok = np.isfinite(e1)
    acc = float((np.where(e1[ok] > 0, 1, -1) == truth.compartment[ok]).mean())
    print(f"  compartiments A/B   accord {acc:.1%} sur {ok.sum():,} bins"
          f"   ·  signe orienté par la piste GC")

    ins = tad_boundaries(clr, args.window)
    called = np.flatnonzero(ins["is_boundary"].fillna(False).to_numpy())
    a = match_positions(called, truth.boundaries, tol=1)
    print(f"  frontières de TAD   {a}   ·  fenêtre {args.window // 1000} kb, ±1 bin")

    d = loops(clr)
    pairs = np.c_[d["start1"].to_numpy() // clr.binsize, d["start2"].to_numpy() // clr.binsize]
    b = match_pairs(pairs, truth.loops, tol=2)
    print(f"  boucles             {b}   ·  ±2 bins")

    worst = min(acc, a.f1, b.f1)
    print(f"\n  Le plus faible des accords : {worst:.0%}. "
          f"En dessous de 80 %, c'est un bug du caller, pas un mauvais jour.")
    return 0 if worst >= 0.8 else 1


def cmd_recon(args: argparse.Namespace) -> int:
    """Balaie l'exposant de conversion contact → distance contre une géométrie connue.

    Le Hi-C réel ne permet pas cette mesure : la structure 3D y est précisément
    l'inconnue. Ici on fabrique la conformation, on en dérive les contacts par un
    modèle direct d'exposant `gamma`, et on regarde quel `alpha` la restitue.
    """
    import logging as _logging
    import warnings

    _logging.disable(_logging.INFO)
    warnings.simplefilter("ignore")
    try:
        import scipy  # noqa: F401
    except ImportError:
        print("scipy absent — voir docs/SETUP.md", file=sys.stderr)
        return 2

    from .hic.polymer import chain, contacts
    from .hic.reconstruct import sweep

    alphas = np.round(np.arange(args.lo, args.hi + 1e-9, args.step), 4)
    inv = 1.0 / args.gamma

    conf0 = chain(n=args.n, seed=0)
    r = conf0.radial()
    print(f"conformation   {args.n} billes, {args.conformations} tirages")
    print(
        f"               radial A {r[conf0.compartment == 1].mean():.3f} vs "
        f"B {r[conf0.compartment == -1].mean():.3f}  "
        f"— A central, B périphérique (Cremer & Cremer)"
    )
    print(f"modèle direct  f ∝ d^(-{args.gamma})   →   inversion exacte : alpha = {inv:.3f}\n")

    print("  profondeur     densité   alpha* médian   étendue        nRMSD min   plateau +5%")
    print("  ------------   -------   -------------   ------------   ---------   -----------")

    for total in args.depths:
        stars, mins, widths, dens = [], [], [], []
        for s in range(args.conformations):
            conf = chain(n=args.n, seed=s)
            counts = contacts(conf, gamma=args.gamma, total=total, seed=100 + s)
            dens.append((counts > 0).sum() / (counts.size - len(counts)))
            nr = np.array([f.nrmsd for f in sweep(counts, conf.coords, alphas)])
            stars.append(float(alphas[nr.argmin()]))
            mins.append(float(nr.min()))
            ok = alphas[nr <= nr.min() * 1.05]
            widths.append(float(ok.max() - ok.min()))
        print(
            f"  {total:>12,}   {np.mean(dens):>6.0%}   {np.median(stars):>13.3f}   "
            f"[{min(stars):.2f}–{max(stars):.2f}]{'':4}   {np.median(mins):>9.3f}   "
            f"{np.median(widths):>11.3f}"
        )

    print(
        f"\n  alpha* décroît vers {inv:.3f} avec la profondeur, sans jamais l'atteindre à\n"
        f"  profondeur finie. Et le plateau s'élargit quand les données se creusent :\n"
        f"  là où il faudrait le plus calibrer alpha, c'est là qu'il est le moins\n"
        f"  déterminé. Reprendre alpha = 1/3 d'un article sans regarder sa profondeur\n"
        f"  de séquençage n'est pas une convention, c'est une approximation non chiffrée."
    )
    return 0


def cmd_nucleus(args: argparse.Namespace) -> int:
    """Construit un noyau diploïde entier et rend compte de ce qu'il vaut.

    La semaine 5 reconstruisait une chaîne. Ici : 46 chaînes, deux mètres d'ADN
    diploïde, une sphère de dix micromètres, et rien qui se traverse.
    """
    try:
        import scipy  # noqa: F401
    except ImportError:
        print("scipy absent — voir docs/SETUP.md", file=sys.stderr)
        return 2

    from .nucleus import build, capacity_at, gm12878, save

    karyotype = gm12878()
    print(f"caryotype      {karyotype}")
    if karyotype.provenance == "builtin":
        print("               ↑ longueurs de la table interne : le fichier officiel n'a "
              "jamais pu être récupéré (voir docs/DATA_SOURCES.md § 8)")

    print(f"\n  La coquille de contact d'un modèle à billes ne tient qu'une monocouche, "
          f"et\n  son volume suit le rayon des billes. Ce que « être à la lamina » peut "
          f"vouloir\n  dire dépend donc de la résolution du modèle, à densité nucléaire "
          f"fixée ({args.phi:.0%}) :\n")
    print("  résolution     billes      rayon    densité uniforme   borne d'empilement")
    print("  ----------   ---------   --------   ----------------   ------------------")
    for bp in (3_000_000, args.bp_per_bead, 250_000, 100_000, 10_000):
        n, r_nm, cap = capacity_at(
            karyotype.total_bp, bp, nuclear_radius=args.nuclear_radius, phi=args.phi
        )
        # Part de volume exacte de la coquille de contact dans la boule
        # accessible aux centres — pas son approximation 1,5·r/(R−r), pour que
        # ce tableau et celui d'ARCHITECTURE.md § 10 donnent les mêmes chiffres.
        free = (args.nuclear_radius - r_nm / 1000.0) ** 3
        u = (free - (args.nuclear_radius - 1.5 * r_nm / 1000.0) ** 3) / free
        print(f"  {bp // 1000:>6} kb    {n:>9,}   {r_nm:>6.1f} nm   {u:>15.1%}   {cap:>18.0%}")

    t0 = time.perf_counter()
    nucleus = build(
        karyotype,
        bp_per_bead=args.bp_per_bead,
        nuclear_radius=args.nuclear_radius,
        phi=args.phi,
        outward_gain=args.lamina,
        seed=args.seed,
    )
    elapsed = time.perf_counter() - t0
    b = nucleus.beads

    print(
        f"\nnoyau          {b.n:,} billes · {int(np.median(b.end - b.start)) // 1000} kb "
        f"chacune · rayon {b.radius.mean() * 1000:.0f} nm · phi {b.phi:.0%} · "
        f"R = {b.nuclear_radius:.1f} µm\n"
        f"               construit en {elapsed:.0f} s, graine {nucleus.seed}"
    )

    if not nucleus.sealed:
        print(
            "\n  ! le tube de chaîne n'est pas étanche pendant l'inflation "
            "(inflate_from ≤ stretch/2) :\n"
            "    une chaîne peut traverser une liaison, et le défaut est difficile à défaire."
        )

    print(f"\n  conditions          {nucleus.final}")
    print(f"  territorialité      {nucleus.after}")
    print(f"                      à l'initialisation {nucleus.before.index:.1f}× — "
          f"le recuit en conserve {nucleus.after.index / nucleus.before.index:.0%}")
    print(f"  périphérie          {nucleus.rim}")

    if args.out:
        print(f"\n  écrit               {save(nucleus, Path(args.out))}")

    # Deux verdicts distincts, et seul le premier décide du code de retour.
    # Les conditions dures valent à toute résolution ; la fourchette de billes est
    # une exigence de la *semaine 6*, pas une propriété d'un noyau valide. Les
    # confondre faisait échouer `make nucleus N=250000`, qui produit pourtant un
    # noyau irréprochable — simplement plus fin que ce que la feuille de route
    # demandait cette semaine-là.
    ok = nucleus.final.acceptable(args.tol)
    print(
        f"\n  Conditions dures, jugées à la même tolérance relative ({args.tol:.0%}) — une\n"
        f"  chaîne tendue au-delà de sa limite viole une condition autant qu'un\n"
        f"  chevauchement :\n"
        f"  chevauchement maximal {nucleus.final.max_overlap:.3%} sur "
        f"{nucleus.final.n_overlapping:,} paires en contact · liaison la plus\n"
        f"  tendue +{nucleus.final.bond_stretch:.2%} · {nucleus.final.outside} bille "
        f"hors du noyau. {'Tenues.' if ok else 'NON TENUES.'}"
    )

    if 6_000 <= b.n <= 10_000:
        print(
            f"\n  Critère semaine 6 — un noyau diploïde de 6 000 à 10 000 billes TAD, sans\n"
            f"  interpénétration : {b.n:,} billes. {'Atteint.' if ok else 'NON ATTEINT.'}"
        )
    else:
        print(
            f"\n  {b.n:,} billes, hors de la fourchette 6 000–10 000 de la semaine 6 : c'est\n"
            f"  une autre résolution, pas un échec. Le critère de la semaine ne s'y applique\n"
            f"  pas ; les conditions dures ci-dessus, si."
        )
    print(
        f"\n  Ce que ça ne dit pas : la territorialité est *entrée* dans le modèle par\n"
        f"  l'initialisation — une relaxation ne fait jamais se croiser deux chaînes.\n"
        f"  Le seul énoncé honnête est que le recuit la conserve. Et la stratification\n"
        f"  radiale vient d'un terme du modèle, pas d'une mesure : la piste LAD porte\n"
        f"  la source « {b.lad_source} ». La confrontation au DamID publié est le\n"
        f"  livrable de la semaine 7, et elle attend le réseau."
    )
    return 0 if ok else 1


def cmd_ensemble(args: argparse.Namespace) -> int:
    """Produit N repliements du même génome et rend compte de ce qu'ils ont en commun.

    Le principe n° 1 du projet dit qu'une structure unique est un artefact
    statistique. Ici il devient un nombre : quelle part de « la bille i est à
    telle profondeur » est un énoncé sur la bille, et quelle part sur le tirage.
    """
    try:
        import scipy  # noqa: F401
        import zarr  # noqa: F401
    except ImportError as exc:
        print(f"dépendance absente ({exc.name}) — voir docs/SETUP.md", file=sys.stderr)
        return 2

    import numpy as np

    from .ensemble import (
        REGIMES,
        contacts,
        frame,
        generate,
        medoid,
        pending,
        radial_by_quartile,
        read,
        reproducibility,
        self_consistency,
    )

    out = Path(args.out)
    if not args.report_only:
        def progress(k: int, total: int, elapsed: float, quality) -> None:
            if k % args.every == 0 or k == total:
                left = elapsed / k * (total - k)
                print(
                    f"  {k:>4}/{total}   {elapsed / 60:5.1f} min écoulées, "
                    f"~{left / 60:4.1f} restantes   (dernière : {quality.max_overlap:.3%})",
                    flush=True,
                )

        done = args.n - len(pending(out)) if (out / ".zgroup").exists() and not args.fresh else 0
        print(
            f"génération   {args.n} structures × {args.workers or 'tous les'} cœurs\n"
            f"             génome figé par lad_seed={args.lad_seed}, "
            f"conformations depuis {args.first_seed}"
            + (f"\n             {done} déjà en magasin, on complète" if done else "")
            + "\n"
        )
        generate(
            out,
            n_structures=args.n,
            workers=args.workers,
            lad_seed=args.lad_seed,
            first_seed=args.first_seed,
            fresh=args.fresh,
            progress=progress,
            bp_per_bead=args.bp_per_bead,
        )

    e = read(out)
    b = e.beads
    radial = e.radial()
    bp = float(np.median(b.end - b.start))

    print(f"\nensemble     {e}")
    print(f"             {out}")
    print(
        f"\n  conditions          chevauchement max sur l'ensemble "
        f"{e.max_overlap.max():.3%} · liaison la plus tendue +{e.bond_stretch.max():.2%} ·\n"
        f"                      {int(e.outside.sum())} bille hors du noyau · "
        f"{int(e.shakes.sum())} secousses au total, "
        f"{int((e.shakes > 0).sum())} structures concernées"
    )

    print("\n  — Il n'y a pas de repère commun —")
    print(f"  {frame(e.coords, b.nuclear_radius, pairs=args.frame_pairs)}")
    print(
        "  Deux noyaux recuits séparément ne partagent ni orientation ni placement des\n"
        "  territoires. Une variance par bille en x, y, z serait donc un nombre sans objet ;\n"
        "  tout ce qui suit ne manipule que des grandeurs invariantes par rotation."
    )

    rep = reproducibility(radial, b.nuclear_radius)
    print("\n  — Ce qui se reproduit d'un tirage à l'autre —")
    print(f"  {rep}")
    print(
        f"  Autrement dit, {rep.icc:.0%} de la variance de profondeur tient à la bille et\n"
        f"  {1 - rep.icc:.0%} au tirage. Une structure isolée porte donc les deux, sans les\n"
        f"  distinguer — c'est exactement ce que le principe n° 1 interdit d'afficher seul."
    )

    print("\n  — Distributions radiales par quartile de contenu LAD —")
    print("  quartile          fraction LAD    rayon moyen    dispersion entre billes")
    for k, (mean, sd, lo, hi) in enumerate(radial_by_quartile(radial, b.lad), start=1):
        tag = {1: "le moins LAD", 4: "le plus LAD"}.get(k, "")
        print(
            f"  Q{k} {tag:<13} {lo:>6.2f}–{hi:<6.2f} {mean:>11.3f}    {sd:>10.3f}"
        )
    print(f"  auto-cohérence avec la piste LAD d'entrée : r = {self_consistency(radial, b.lad):+.3f}")
    print(
        "  Ce r ne valide rien : la piste LAD est ce qui a fixé les rayons visés. Il dit\n"
        "  seulement que le solveur a fait ce qu'on lui demandait."
    )

    index, total = medoid(radial)
    print("\n  — Structure médoïde —")
    print(
        f"  index {index}, graine {int(e.seeds[index])} · distance cumulée "
        f"{total[index]:.1f} contre {total.max():.1f} pour la plus excentrée\n"
        f"  Médoïde **pour la distance entre profils radiaux** : deux structures aux\n"
        f"  territoires disposés tout autrement peuvent avoir le même profil."
    )

    print("\n  — Ce que l'ensemble prédit d'une expérience Hi-C —")
    print("  seuil   contacts/structure     trans   homologues   plateau P(s)")
    print("  -----   ------------------   -------   ----------   ------------")
    main: object | None = None
    for cut in args.cutoffs:
        sub = e.coords if cut == args.cutoffs[0] else e.coords[: args.sweep]
        c = contacts(sub, b.radius, b.copy_id, b.labels, bp, cutoff=cut)
        main = c if main is None else main
        mark = "" if cut == args.cutoffs[0] else f"  ({len(sub)} structures)"
        print(
            f"  {cut:>5.2f}   {c.per_structure:>18,.0f}   {c.trans_fraction:>6.1%}   "
            f"{c.homolog / max(c.trans, 1):>9.1%}   {c.plateau:>12.4f}{mark}"
        )
    print(
        f"\n  Un seuil sous {1.15:.2f} ne mesurerait rien : c'est l'allongement maximal d'une\n"
        f"  liaison, donc en dessous même deux billes voisines de chaîne ne « se touchent »\n"
        f"  pas et il ne reste que les chevauchements résiduels.\n"
    )

    print("  P(s) n'a pas une pente, elle en a trois :")
    print("  régime                        domaine      pente")
    print("  ------------------------   ------------   -------")
    for (name, sl), (_, lo, hi) in zip(main.regimes, REGIMES):
        print(f"  {name:<24}   {lo / 1e6:>4.1f}–{hi / 1e6:<5.0f} Mb   {sl:>+7.2f}")
    print(
        f"  puis un plateau à P = {main.plateau:.4f} au-delà de 15 Mb.\n\n"
        f"  Aucune matrice de contacts n'a été montrée au modèle : P(s) sort de la seule\n"
        f"  géométrie. C'est donc la seule grandeur de cet ensemble qu'une expérience Hi-C\n"
        f"  puisse contredire — et elle la contredit. Le Hi-C réel décroît en ~s^-1 de façon\n"
        f"  continue de la centaine de kb à la dizaine de Mb. Ce modèle s'en approche au\n"
        f"  régime polymère, puis **s'aplatit** au-delà de 15 Mb au lieu de continuer à\n"
        f"  décroître. Ses territoires sont trop bien mélangés à l'intérieur : il reproduit\n"
        f"  le *fait* des territoires, pas leur organisation interne. C'est précisément ce\n"
        f"  que la semaine 8 — extrusion de boucles, polymère fin — doit apporter.\n\n"
        f"  Le hasard pur donnerait une fraction trans de {1 - 1 / len(b.labels):.0%} et des\n"
        f"  homologues à {1 / (len(b.labels) - 1):.1%} des trans."
    )

    print("\n  — Position radiale contre DamID mesuré —")
    if args.damid:
        from .ensemble import damid, read_bedgraph

        track, rows = read_bedgraph(args.damid)
        d = damid(radial, b.start, b.end, b.copy_id, b.labels, track, str(args.damid))
        print(f"  {rows:,} intervalles lus sur {len(track)} chromosomes")
        print(f"  {d}")
        if d.covered < 0.5 * d.total:
            print(
                "  ! moins de la moitié des billes sont couvertes : la corrélation porte sur\n"
                "    un sous-ensemble, et il faut dire lequel avant de la citer."
            )
    else:
        print(
            "  La feuille de route demande la corrélation entre position radiale modélisée\n"
            "  et LADs DamID **publiés**. Aucune entrée DamID du manifeste n'a pu être\n"
            "  récupérée (docs/DATA_SOURCES.md § 8), et aucun fichier n'a été fourni.\n"
            "  Le calcul est écrit, testé sur une piste construite exprès, et se lance par\n"
            "  `make ensemble DAMID=chemin.bedGraph` — il lui faut un fichier, pas une\n"
            "  ligne de code de plus."
        )
    return 0


def cmd_fine(args: argparse.Namespace) -> int:
    """Descend sous le TAD : extrusion de boucles sur une région, et le raccord au noyau.

    La semaine 7 a laissé une P(s) plate au-delà de 15 Mb là où le Hi-C réel
    décroît. Cette commande construit l'échelle où cette organisation naît, et
    répond à trois questions que le noyau seul ne pouvait pas poser : les
    barrières orientées produisent-elles des domaines, produisent-elles des
    points d'angle **aux seules paires convergentes**, et la chaîne fine
    tombe-t-elle sur la même `R(s)` que le noyau là où les deux se recouvrent ?
    """
    try:
        import openmm  # noqa: F401
        import scipy  # noqa: F401
    except ImportError as exc:
        print(f"dépendance absente ({exc.name}) — voir docs/SETUP.md", file=sys.stderr)
        return 2

    import numpy as np

    from .fine import (
        CONVERGENT,
        HIC_REFERENCE,
        KINDS,
        REGIMES,
        Region,
        call_boundaries,
        contact_map,
        dots_by_kind,
        insulation,
        junction,
        ps,
        regime,
        separation_curve,
        stationarity,
    )
    from .fine.polymer import Field
    from .fine.run import Setup, build, knockout, load, region_of, save
    from .hic.features import match_positions

    out = Path(args.out)
    region = Region(args.chrom, args.start, args.end, args.bp_per_bead)
    setup = Setup(
        region=region,
        separation=args.separation,
        processivity=args.processivity,
        release=args.release,
        md_per_step=args.md_per_step,
        relax=args.relax,
        snapshots=args.snapshots,
        stride=args.stride,
        replicates=args.replicates,
        cutoff=args.cutoff,
        field=Field(trunc=args.trunc, stiffness=args.stiffness),
    )

    print(f"région       {region}")
    print(
        f"             confinement {setup.confine:.2f} sigma "
        f"({setup.confine * setup.sigma_nm:.0f} nm de rayon) à phi = {setup.phi}\n"
        f"             sigma = {setup.sigma_nm:.1f} nm — loi de la semaine 6, "
        f"pas un réglage de la semaine 8"
    )

    if args.report_only:
        coords, sites, truth, meta = load(str(out))
        region = region_of(meta)
        ko_coords = None
        if Path(str(out).replace(".npz", "") + "-ko.npz").exists():
            ko_coords, _, _, _ = load(str(out).replace(".npz", "") + "-ko.npz")
        elapsed = meta.get("elapsed_s", 0.0)
        occupancy = meta.get("lef_occupancy", float("nan"))
    else:
        def progress(k: int, total: int) -> None:
            print(f"  réplicat {k}/{total}", flush=True)

        given = None
        if args.ctcf:
            from .fine import read_ctcf

            given = read_ctcf(args.ctcf, region)
            print(f"\nsites CTCF   {given.k} motifs orientés lus dans {args.ctcf}")

        print(f"\ngénération   {args.replicates} réplicats × {args.snapshots} instantanés")
        fine = build(
            setup,
            sites=given,
            site_seed=args.site_seed,
            first_seed=args.first_seed,
            workers=args.workers,
            progress=progress,
            **(
                {}
                if given is not None
                else dict(
                    mean_domain=args.mean_domain,
                    min_domain=args.min_domain,
                    control=args.control,
                )
            ),
        )
        out.parent.mkdir(parents=True, exist_ok=True)
        save(fine, str(out))
        coords, sites, truth = fine.coords, fine.sites, fine.truth
        elapsed, occupancy = fine.elapsed, fine.lef_occupancy

        ko_coords = None
        if not args.no_knockout:
            print("\ntémoin       les mêmes sites, tous inoccupés")
            ko = knockout(
                setup, sites, first_seed=args.first_seed, workers=args.workers,
                progress=progress,
            )
            save(ko, str(out).replace(".npz", "") + "-ko.npz")
            ko_coords = ko.coords
            elapsed += ko.elapsed

    print(
        f"\n             {len(coords)} conformations en {elapsed / 60:.1f} min · "
        f"cohésines chargées {occupancy:.0%} du temps\n             {out}"
    )

    binned = max(1, args.bin_kb * 1_000 // region.bp_per_bead)
    cmap = contact_map(
        coords, cutoff=args.cutoff, bin_beads=binned, bp_per_bead=region.bp_per_bead
    )
    print(
        f"\ncarte        {cmap.m}×{cmap.m} casiers de {cmap.bp_per_bin // 1000} kb · "
        f"seuil {cmap.cutoff} diamètre (convention de la semaine 7)"
    )

    # — Frontières —
    window = max(2, args.insulation_kb * 1_000 // cmap.bp_per_bin)
    ins = insulation(cmap, window)
    called = call_boundaries(ins, prominence=args.prominence)
    if truth is not None:
        planted = cmap.bin_of(truth.boundary[truth.strong])
        weak = cmap.bin_of(truth.boundary[~truth.strong])
        print("\n  — Les barrières font-elles des domaines ? —")
        print(
            f"  fenêtre d'insulation {window} casiers "
            f"({window * cmap.bp_per_bin // 1000} kb), tolérance ±1 casier"
        )
        print(
            f"  frontières à deux sens bloqués ({len(planted)}/{len(truth.boundary)}) : "
            f"{match_positions(called, planted, tol=1)}"
        )
        if len(weak):
            found = match_positions(called, weak, tol=1)
            print(
                f"  frontières à un seul sens bloqué  ({len(weak)}/{len(truth.boundary)}) : "
                f"{found.n_matched}/{found.n_true} retrouvées"
            )
        print(
            "  Il n'existe pas de frontière non bloquante dans ce dispositif : elle porte\n"
            "  deux ancres, et une ancre arrête toujours l'un des deux sens de marche. Ce qui\n"
            "  les distingue est le nombre de sens arrêtés, pas l'existence d'une barrière."
        )
        if args.sweep and ko_coords is not None:
            # Le seuil d'appel est un choix, et il se justifie sur le **témoin** : à
            # seuil égal, ce que retrouve une carte sans aucun CTCF occupé est ce que
            # le hasard du dispositif donne. Un réglage dont le témoin retrouve 42 %
            # des frontières est disqualifié quel que soit son rappel.
            ko_map = contact_map(
                ko_coords, cutoff=args.cutoff, bin_beads=binned,
                bp_per_bead=region.bp_per_bead,
            )
            print("\n  — Le seuil d'appel se choisit sur le témoin, pas sur le rappel —")
            print(
                f"  {'fenêtre':>8} {'proém.':>7} {'appelées':>9} {'rappel':>7} {'préc.':>7}"
                f"   {'témoin app.':>12} {'témoin rap.':>12}"
            )
            all_planted = cmap.bin_of(truth.boundary)
            for kb in (30, 50, 100, 200):
                w = max(2, kb * 1_000 // cmap.bp_per_bin)
                a_ins = insulation(cmap, w)
                b_ins = insulation(ko_map, w)
                for prom in (0.05, 0.08, 0.12, 0.18):
                    a = match_positions(
                        call_boundaries(a_ins, prominence=prom), all_planted, tol=1
                    )
                    b = match_positions(
                        call_boundaries(b_ins, prominence=prom), all_planted, tol=1
                    )
                    print(
                        f"  {kb:>6} kb {prom:>7.2f} {a.n_called:>9} {a.recall:>7.0%} "
                        f"{a.precision:>7.0%}   {b.n_called:>12} {b.recall:>12.0%}"
                    )

    # — Points d'angle —
    if truth is not None:
        print("\n  — Les points d'angle sont-ils aux paires convergentes ? —")
        dots = dots_by_kind(cmap, truth, inner=args.dot_inner, outer=args.dot_outer)
        print(f"  {'classe':<12} {'n':>4} {'médiane':>9} {'moyenne':>9}")
        for name, n, med, mean in dots.summary():
            print(f"  {name:<12} {n:>4} {med:>9.2f} {mean:>9.2f}")
        print(
            f"  convergents / divergents : {dots.separation('divergent'):.2f}×   ·   "
            f"convergents / hasard : {dots.separation('hasard'):.2f}×"
        )
        print(
            "  Le témoin propre est **divergent** : aucune de ses deux ancres n'arrête la\n"
            "  jambe qui l'atteint. Un domaine en tandem n'est témoin qu'à moitié — si le\n"
            "  domaine suivant commence par un `−`, ce `−` se trouve à un monomère de son\n"
            "  ancre droite et y ancre une vraie boucle, que la tolérance du score ramasse."
        )
        if ko_coords is not None:
            ko_map = contact_map(
                ko_coords, cutoff=args.cutoff, bin_beads=binned,
                bp_per_bead=region.bp_per_bead,
            )
            ko_dots = dots_by_kind(ko_map, truth, inner=args.dot_inner, outer=args.dot_outer)
            ko_ins = insulation(ko_map, window)
            ko_called = call_boundaries(ko_ins, prominence=args.prominence)
            print("\n  — Le témoin qui tranche : les mêmes sites, tous inoccupés —")
            print(
                f"  points d'angle aux convergents : "
                f"{dots.median(KINDS[CONVERGENT]):.2f} → "
                f"**{ko_dots.median(KINDS[CONVERGENT]):.2f}**"
            )
            print(f"  frontières : {match_positions(ko_called, planted, tol=1)}")
            print(
                "  Rien d'autre ne change entre les deux exécutions — mêmes graines, mêmes\n"
                "  positions de sites, même champ de force. La différence est donc bien\n"
                "  l'arrêt de l'extrusion par CTCF, et rien d'autre. Ce que le témoin retrouve\n"
                "  encore de frontières mesure le hasard du dispositif, et il faut le soustraire\n"
                "  mentalement du rappel ci-dessus."
            )

    # — P(s) —
    s, p = ps(cmap)
    print("\n  — P(s) sous le TAD —")
    print(f"  {'régime':<14} {'fenêtre':<20} {'pente':>7}")
    for name, lo, hi in REGIMES:
        print(f"  {name:<14} {lo / 1e3:,.0f}–{hi / 1e3:,.0f} kb{'':<6} {regime(s, p, lo, hi):>+7.2f}")
    print(
        "  Seul le régime « domaine » est comparable au Hi-C publié : Lieberman-Aiden 2009\n"
        f"  rapporte s^{HIC_REFERENCE:.2f} sur 500 kb–7 Mb — une fenêtre qu'une région de 4 Mb ne\n"
        "  couvre pas, donc les deux nombres se ressemblent sans être mesurés pareil.\n"
        "  Au-delà d'environ 1,5 Mb, une région\n"
        "  de 4 Mb confinée dans sa propre sphère mesure la sphère — c'est un artefact de la\n"
        "  taille de région, pas un résultat, et le nommer évite de republier la P(s) plate de\n"
        "  la semaine 7 (−0,10 au-delà de 15 Mb) comme si c'en était une prédiction."
    )

    # — Raccord —
    if args.nucleus and Path(args.nucleus).exists():
        from .ensemble import read as read_ensemble

        e = read_ensemble(args.nucleus)
        take = min(args.nucleus_structures, e.n_structures)
        coarse_bp = int(np.median(e.beads.end - e.beads.start))
        try:
            j = junction(
                coords,
                region,
                e.coords[:take],
                coarse_copy_id=e.beads.copy_id,
                coarse_bp_per_bead=coarse_bp,
                sigma_nm=setup.sigma_nm,
            )
        except ValueError as exc:
            j = None
            print(f"\n  — Raccord non mesurable —\n  {exc}")
        if j is not None:
            print(f"\n  — Le raccord avec le noyau entier ({take} structures) —")
            print(
                f"  {'séparation':>12} {'fin (nm)':>10} {'noyau (nm)':>12} "
                f"{'rapport':>9} {'paires/conf.':>13}"
            )
            for bp, f, c, r, np_, thin in zip(
                j.bp, j.fine_nm, j.coarse_nm, j.ratio, j.fine_pairs, j.thin
            ):
                mark = " ·" if thin else ""
                print(f"  {bp / 1e6:>9.2f} Mb {f:>10.0f} {c:>12.0f} {r:>9.2f} {np_:>13,}{mark}")
            if j.thin.any():
                print(
                    "  · peu de paires : à la séparation s il n'en reste que n − s,\n"
                    "    donc le dernier point de recouvrement est le plus fragile."
                )
            print(
                f"\n  Aucun paramètre n'a été réglé sur l'autre modèle : le rayon d'un monomère\n"
                f"  sort de la loi de la semaine 6, `phi` est le même, et la semaine 8 n'a rien\n"
                f"  ajusté. Écart le plus grand {j.worst:.2f}× — "
                + ("le raccord tient." if j.holds(args.junction_tol) else
                   f"au-delà du facteur {args.junction_tol} admis.")
            )
            fine_e, coarse_e = j.exponents()
            print(
                f"\n  Ce que le tableau dit vraiment est une différence de **pente** :\n"
                f"  R(s) ∝ s^{fine_e:.2f} pour le modèle fin, s^{coarse_e:.2f} pour le noyau.\n"
                f"  Repères : 0,50 pour une marche aléatoire idéale, 0,59 pour une marche\n"
                f"  auto-évitante gonflée, 0,25 à 0,33 pour ce que le traçage de chromatine\n"
                f"  mesure au-dessus du mégabase. Les deux modèles encadrent la mesure au lieu\n"
                f"  de la reproduire, et ils l'encadrent par les deux bouts."
            )
            print(
                f"\n  Une réserve, et elle porte sur le haut du tableau : le confinement du\n"
                f"  modèle fin est la sphère que la loi de volume de la semaine 6 alloue à\n"
                f"  {region.n * region.bp_per_bead / 1e6:.0f} Mb, soit {setup.confine * setup.sigma_nm:.0f} nm de rayon. "
                f"Au-delà d'environ 1 Mb, `R(s)` y sature\n"
                f"  donc contre une paroi dont la taille vient du modèle grossier lui-même, et\n"
                f"  l'accord y devient partiellement circulaire. Les lignes informatives sont\n"
                f"  celles nettement sous cette échelle."
            )
    elif args.nucleus:
        print(
            f"\n  — Raccord non mesuré —\n"
            f"  {args.nucleus} n'existe pas. `make ensemble` d'abord, puis `make fine`."
        )

    # — Stationnarité : la trajectoire a-t-elle oublié son point de départ ? —
    st = stationarity(
        coords,
        per_replicate=max(1, len(coords) // max(1, args.replicates)),
        bp_per_bead=region.bp_per_bead,
        unit_nm=setup.sigma_nm,
        seps=np.array([10, 100, 400, 1000, int(0.9 * region.n)]),
    )
    print("\n  — La trajectoire a-t-elle oublié sa conformation initiale ? —")
    print(f"  {'séparation':>12} {'1re moitié':>12} {'2e moitié':>12} {'dérive':>8}")
    for bp, a, b, d in zip(st.bp, st.first_nm, st.second_nm, st.drift):
        unit = f"{bp / 1e6:.2f} Mb" if bp >= 1e6 else f"{bp / 1e3:.0f} kb"
        print(f"  {unit:>12} {a:>9.0f} nm {b:>9.0f} nm {d:>+8.1%}")
    print(
        f"  Dérive la plus forte {st.worst:.1%} — "
        + ("stationnaire sur la fenêtre d'échantillonnage."
           if st.holds() else
           "la conformation initiale décide encore ; `--relax` est trop court.")
    )
    print(
        "  C'est `R(s)` aux grandes séparations qu'il faut regarder, pas `Rg` : une\n"
        "  conformation initiale comprimée a déjà le bon `Rg` par construction."
    )

    # — R(s) de la chaîne fine seule, pour mémoire —
    bp, mean, sd = separation_curve(
        coords, bp_per_bead=region.bp_per_bead, unit_nm=setup.sigma_nm
    )
    print("\n  — R(s) de la chaîne fine —")
    print(f"  {'séparation':>12} {'moyenne':>9} {'écart-type':>11}")
    for x, m, d in list(zip(bp, mean, sd))[::4]:
        unit = f"{x / 1e6:.2f} Mb" if x >= 1e6 else f"{x / 1e3:.0f} kb"
        print(f"  {unit:>12} {m:>7.0f} nm {d:>9.0f} nm")

    print(
        "\n  Sites CTCF : "
        + (
            f"**plantés** ({sites.source}) — aucune piste publiée n'est lisible ici,\n"
            "  le réseau est fermé. Comme à la semaine 3, la vérité plantée valide la\n"
            "  méthode et rien de biologique. `make fine CTCF=motifs.bed` attend un fichier."
            if sites.source.startswith("planted")
            else sites.source
        )
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="geno", description="Socle 1D du génome — magasin d'intervalles.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build", help="construit le magasin depuis un tracks.json")
    b.add_argument("--tracks", default=str(DEFAULT_TRACKS))
    b.add_argument("--out", default=str(DEFAULT_STORE))
    b.set_defaults(fn=cmd_build)

    q = sub.add_parser("query", help="interroge une région (1-based, bornes incluses)")
    q.add_argument("region", help="ex. chr7:5,527,000-5,530,600")
    q.add_argument("--store", default=str(DEFAULT_STORE))
    q.add_argument("--track", action="append", help="limiter à cette piste (répétable)")
    q.add_argument(
        "--flank",
        type=int,
        default=0,
        metavar="PB",
        help="élargit la fenêtre de PB de chaque côté — une boucle CTCF encadre son gène, "
        "ses ancres sont donc hors des bornes",
    )
    q.add_argument("--json", action="store_true")
    q.add_argument("--time", action="store_true", help="affiche la latence")
    q.set_defaults(fn=cmd_query)

    t = sub.add_parser("tracks", help="liste les pistes du magasin")
    t.add_argument("--store", default=str(DEFAULT_STORE))
    t.set_defaults(fn=cmd_tracks)

    h = sub.add_parser(
        "hic", help="plante une structure Hi-C connue et valide les callers dessus"
    )
    h.add_argument("--bins", type=int, default=1_000)
    h.add_argument("--resolution", type=int, default=10_000)
    h.add_argument("--window", type=int, default=100_000, help="fenêtre d'insulation")
    h.add_argument("--seed", type=int, default=3)
    h.add_argument("--out", default=str(ROOT / "data" / "synthetic" / "planted.cool"))
    h.set_defaults(fn=cmd_hic)

    rc = sub.add_parser(
        "recon", help="balaie l'exposant contact → distance contre une géométrie connue"
    )
    rc.add_argument("--n", type=int, default=300, help="billes de la chaîne")
    rc.add_argument("--gamma", type=float, default=3.0, help="exposant du modèle direct")
    rc.add_argument("--conformations", type=int, default=3)
    rc.add_argument("--lo", type=float, default=0.15)
    rc.add_argument("--hi", type=float, default=0.80)
    rc.add_argument("--step", type=float, default=0.025)
    rc.add_argument(
        "--depths",
        type=int,
        nargs="+",
        default=[500_000, 2_000_000, 8_000_000, 40_000_000, 200_000_000],
    )
    rc.set_defaults(fn=cmd_recon)

    nu = sub.add_parser("nucleus", help="construit un noyau diploïde complet de billes TAD")
    nu.add_argument("--bp-per-bead", type=int, default=750_000,
                    help="taille génomique d'une bille — 750 kb est l'échelle TAD de Dixon")
    nu.add_argument("--nuclear-radius", type=float, default=5.0, metavar="µm")
    nu.add_argument("--phi", type=float, default=0.30,
                    help="fraction du volume nucléaire occupée par les billes")
    nu.add_argument("--lamina", type=float, default=0.05,
                    help="force du rappel radial vers la périphérie (0 = aucun)")
    nu.add_argument("--tol", type=float, default=0.01, help="chevauchement maximal toléré")
    nu.add_argument("--seed", type=int, default=0)
    nu.add_argument("--out", default=str(ROOT / "data" / "nucleus" / "gm12878.npz"))
    nu.set_defaults(fn=cmd_nucleus)

    en = sub.add_parser(
        "ensemble", help="produit N repliements du même génome et les caractérise"
    )
    en.add_argument("--n", type=int, default=200, help="nombre de structures")
    en.add_argument("--workers", type=int, default=0, help="0 = tous les cœurs")
    en.add_argument("--lad-seed", type=int, default=0,
                    help="graine du *génome* — fixe pour tout l'ensemble")
    en.add_argument("--first-seed", type=int, default=1_000,
                    help="première graine de *conformation*")
    en.add_argument("--bp-per-bead", type=int, default=750_000)
    en.add_argument("--fresh", action="store_true",
                    help="efface un magasin existant et recommence ; par défaut on le complète")
    en.add_argument("--report-only", action="store_true",
                    help="ne produit rien, se contente de relire et rendre compte")
    en.add_argument("--cutoffs", type=float, nargs="+", default=[1.5, 1.25, 2.0],
                    help="seuils de contact ; le premier sert à l'ensemble entier")
    en.add_argument("--sweep", type=int, default=50,
                    help="structures utilisées pour les seuils secondaires")
    en.add_argument("--frame-pairs", type=int, default=20)
    en.add_argument("--every", type=int, default=10, help="cadence des lignes d'avancement")
    en.add_argument("--damid", default=None, metavar="BEDGRAPH",
                    help="piste DamID mesurée, pour la corrélation demandée par la semaine 7")
    en.add_argument("--out", default=str(ROOT / "data" / "ensemble" / "gm12878.zarr"))
    en.set_defaults(fn=cmd_ensemble)

    fi = sub.add_parser(
        "fine", help="extrusion de boucles sur une région, et son raccord au noyau"
    )
    # La région par défaut encadre ACTB (chr7:5,527,151-5,530,601), le gène que
    # `make query` affiche depuis la semaine 2 : le modèle fin et la fiche de
    # locus regardent alors le même endroit du génome.
    fi.add_argument("--chrom", default="chr7")
    fi.add_argument("--start", type=int, default=4_000_000)
    fi.add_argument("--end", type=int, default=8_000_000)
    fi.add_argument("--bp-per-bead", type=int, default=2_000,
                    help="résolution du modèle fin, 1–5 kb")
    fi.add_argument("--separation", type=int, default=200_000, help="pb par cohésine")
    fi.add_argument("--processivity", type=int, default=200_000,
                    help="taille de boucle sans obstacle, en pb")
    fi.add_argument("--release", type=float, default=0.003,
                    help="probabilité par pas qu'une jambe quitte un site CTCF")
    fi.add_argument("--replicates", type=int, default=4)
    fi.add_argument("--snapshots", type=int, default=120)
    fi.add_argument("--stride", type=int, default=4, help="pas d'extrusion entre instantanés")
    fi.add_argument("--md-per-step", type=int, default=200)
    fi.add_argument("--relax", type=int, default=150_000,
                    help="pas de mise en place avant le premier instantané")
    fi.add_argument("--trunc", type=float, default=3.0,
                    help="coût kT d'un recouvrement complet — au-delà, les chaînes ne se croisent plus")
    fi.add_argument("--stiffness", type=float, default=1.5, help="terme angulaire, en kT")
    fi.add_argument("--cutoff", type=float, default=1.5,
                    help="seuil de contact en diamètres — convention de la semaine 7")
    fi.add_argument("--bin-kb", type=int, default=10, help="casier de la carte de contacts")
    # La semaine 3 avait montré qu'une fenêtre **trop large** tue le rappel : sur
    # des TADs de 450 kb, 100 % à une fenêtre de 100 kb contre 42 % à 600 kb. J'en
    # avais déduit un rapport de ~0,22 et posé 50 kb pour des domaines de 200 kb.
    # Le balayage mesuré dit l'inverse à ce bout-là : à 100 kb (rapport 0,5), le
    # rappel monte à 89 % pour 94 % de précision, contre 95 %/56 % à 50 kb, et
    # surtout le témoin CTCF inoccupé tombe de 42 % à 5 %. Une fenêtre étroite
    # moyenne moins de pixels et appelle du bruit. La contrainte de la semaine 3
    # tient — rester sous l'échelle des domaines — mais son optimum n'était pas
    # transposable tel quel. Voir VALIDATION.md § S8.
    fi.add_argument("--insulation-kb", type=int, default=100,
                    help="fenêtre d'insulation ; sous l'échelle des domaines, mais pas trop")
    fi.add_argument("--prominence", type=float, default=0.05,
                    help="seuil choisi sur l'écart au témoin, pas sur le rappel")
    fi.add_argument("--dot-inner", type=int, default=1)
    fi.add_argument("--dot-outer", type=int, default=6)
    fi.add_argument("--mean-domain", type=int, default=200_000,
                    help="taille moyenne des domaines plantés — médiane Rao 2014")
    fi.add_argument("--min-domain", type=int, default=80_000)
    fi.add_argument("--control", type=float, default=0.35,
                    help="fraction de domaines témoins, divergents ou en tandem")
    fi.add_argument("--ctcf", default=None, metavar="BED",
                    help="motifs CTCF orientés mesurés ; à défaut on plante une vérité connue")
    fi.add_argument("--site-seed", type=int, default=1,
                    help="graine des *sites* — commune à tous les réplicats")
    fi.add_argument("--first-seed", type=int, default=4_000,
                    help="première graine de *conformation*")
    fi.add_argument("--workers", type=int, default=0)
    fi.add_argument("--no-knockout", action="store_true",
                    help="saute le témoin à CTCF inoccupé — il double le temps de calcul")
    fi.add_argument("--nucleus", default=str(ROOT / "data" / "ensemble" / "gm12878.zarr"),
                    help="ensemble de la semaine 7, pour le raccord R(s)")
    fi.add_argument("--nucleus-structures", type=int, default=25)
    fi.add_argument("--junction-tol", type=float, default=1.25)
    fi.add_argument("--sweep", action="store_true",
                    help="balaye fenêtre × proéminence contre le témoin CTCF inoccupé")
    fi.add_argument("--report-only", action="store_true")
    fi.add_argument("--out", default=str(ROOT / "data" / "fine" / "actb.npz"))
    fi.set_defaults(fn=cmd_fine)

    n = sub.add_parser("bench", help="mesure la latence de requête à l'échelle réelle")
    n.add_argument("--n", type=int, default=1_000_000)
    n.add_argument("--queries", type=int, default=2_000)
    n.add_argument("--store", default=str(BENCH_STORE))
    n.set_defaults(fn=cmd_bench)

    args = ap.parse_args(argv)
    try:
        return args.fn(args)
    except FileNotFoundError as exc:
        print(f"erreur : {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
