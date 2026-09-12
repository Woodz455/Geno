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
