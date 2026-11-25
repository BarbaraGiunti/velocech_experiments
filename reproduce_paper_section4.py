"""
reproduce_paper_section4.py
===========================

Re-runs the experiments of Sections 4.1, 4.2 and 4.3 of

    E. Carlsson and J. Carlsson, "Computing the alpha complex using dual
    active set quadratic programming", arXiv:2310.00536 (2023)

with our implementation (alpha_dualqp.py) and prints a side-by-side
comparison with the numbers reported in the paper.

    python reproduce_paper_section4.py                      # both experiments, both solvers
    python reproduce_paper_section4.py --solvers native --repeats 5 --seeds 10
    python reproduce_paper_section4.py --only 4.1 --check-gudhi
    python reproduce_paper_section4.py --only 4.3 --betti      # + Betti numbers
    python reproduce_paper_section4.py --cpp                    # time the C++ program
    python reproduce_paper_section4.py --max-dim 3              # also tetrahedra

Section 4.1 (Stanford bunny, 3D)
    5000 points from the CGAL website (file `bunny_5000` of CGAL's
    Alpha_shapes_3 example data; downloaded automatically if missing).
    The paper says "radius approximately 1/15 of the diameter" with a maximum
    Cech degree of 200. The radius that reproduces the paper's counts EXACTLY
    (all four of them, checked with GUDHI) and the degree 200 is r = 0.01
    (= diameter / 19.8), so that is the value used here.
    Paper: (|X_k|) = (5000, 23830, 29995, 11163) for k = 0..3, ~12 s.

Section 4.2 (random points, 10D)
    1000 points with coordinates uniform in [-1, 1]^10, Alpha(S, 1.0).
    Paper: max Cech degree 451, (|X_k|) = (1000, 64785, 560459, 1783194),
    ~5 minutes. The paper's random sample is unknown, so we run several
    seeds and report mean / std / range, and how many standard deviations
    the paper's value lies from our mean (|z| < 2: consistent). The counts
    vary a lot between samples (edges +-5%, triangles +-15%).

Section 4.3 (two persistence data sets from Henry Adams' tutorial)
    a) Cyclooctane conformations: 6040 points in R^24, Alpha(S, 0.25).
       The full 6040-point set is Ripser's example file pointsCycloOctane.csv
       (the tutorial ships a 1000-point subset of it). Paper: (6040, 24646,
       28352, 11858), ~7 s, Betti numbers (1, 1, 2).
    b) Optical image patches X(300,30): 1000 points, Alpha(S, 0.5), up to
       triangles only, so here the TIME is directly comparable. File
       optical_k300_points.txt of the tutorial, in DCT coordinates (R^8; an
       isometric copy of the paper's R^9 patches). Paper: (1000, 44080,
       454843), ~38 s, Betti numbers (1, 1).
       Expect a tiny count difference (about 1 edge, 20 triangles): this data
       set lies on the sphere S^7 and has many nearly degenerate Voronoi
       faces, some thinner (~1e-6) than the precision of the data file (~5
       significant digits), so whether they count as present depends on each
       implementation's numerical tolerances. Checked with an exact linear
       program: the native solver's extra triangles do have nonempty (but
       extremely thin) Voronoi faces; DAQP misses 3 of them.
    With --betti, beta_0 and beta_1 of our complexes are computed with GUDHI
    (exact from the 2-skeleton; beta_2 would need tetrahedra).

What is (and is not) comparable
    * Our code stops at triangles (k <= 2); the paper's times include the
      tetrahedra (k = 3), which in Section 4.2 are the largest part of the
      output (1.78 million). Compare COUNTS for k <= 2 exactly; treat the
      time comparison as indicative only.
    * The paper's implementation is in MAPLE (with a compiled MAPLE version
      of DAQP), on unknown hardware.
    * Timing protocol as in run_experiments.py: one timer around one call of
      alpha_complex_dualqp (input points -> ordered filtration), after an
      untimed warm-up that compiles the Numba code.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import sys
import time
import urllib.request

import numpy as np

import alpha_dualqp as aq

# Where the bunny comes from (CGAL's own copy on GitHub) and its checksum,
# so that we know we are using exactly the same file every time.
BUNNY_URL = ("https://raw.githubusercontent.com/CGAL/cgal/master/"
             "Alpha_shapes_3/examples/Alpha_shapes_3/data/bunny_5000")
BUNNY_SHA256 = "697f1d97e9a46e54a0a13ad0954aa1079fa0ab11ea20fad9e88a9eaa10816741"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")          # the data sets of Section 4
BUNNY_FILE = os.path.join(DATA_DIR, "bunny_5000.txt")

# Section 4.3 data: the public files the paper refers to (Henry Adams'
# tutorial, ref. [1]; the full cyclooctane set is in Ripser's examples).
CYCLO_URL = ("https://raw.githubusercontent.com/Ripser/ripser/master/"
             "examples/pointsCycloOctane.csv")
CYCLO_SHA256 = "484e9dbdb8b45204edaba7414a66e21345b13709370126b1635fddfbcc6e7e8d"
CYCLO_FILE = os.path.join(DATA_DIR, "cyclooctane_6040.csv")
OPTICAL_URL = ("https://raw.githubusercontent.com/ds4m/topological-data-analysis/"
               "master/persistent-homology/point_clouds/optical_k300_points.txt")
OPTICAL_SHA256 = "302287756127a557b452b0c85ffe12bf7a8201ab27218c71f08e41891f1afaeb"
OPTICAL_FILE = os.path.join(DATA_DIR, "optical_k300_points.txt")

# Numbers reported in the paper (counts for k = 0, 1, 2, 3; time in seconds).
PAPER = {
    "4.1": dict(counts=(5000, 23830, 29995, 11163), max_degree=200, time_s=12.0,
                time_note="approx., MAPLE, includes tetrahedra"),
    "4.2": dict(counts=(1000, 64785, 560459, 1783194), max_degree=451, time_s=300.0,
                time_note="approx. 5 min, MAPLE, includes tetrahedra"),
    "4.3a": dict(counts=(6040, 24646, 28352, 11858), max_degree=None, time_s=7.0,
                 betti=(1, 1, 2), time_note="approx., MAPLE, includes tetrahedra"),
    "4.3b": dict(counts=(1000, 44080, 454843), max_degree=None, time_s=38.0,
                 betti=(1, 1), time_note="approx., MAPLE, up to triangles (comparable)"),
}


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
def load_bunny(path=BUNNY_FILE):
    """Load the 5000-point bunny (first line = number of points, then x y z),
    downloading it from CGAL's GitHub repository if the file is missing."""
    if not os.path.exists(path):
        print(f"Downloading {BUNNY_URL} ...")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        urllib.request.urlretrieve(BUNNY_URL, path)
    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
    if digest != BUNNY_SHA256:
        print(f"  ! warning: {path} differs from the file used for validation "
              f"(sha256 {digest[:12]}...)")
    X = np.loadtxt(path, skiprows=1)
    assert X.shape == (5000, 3), X.shape
    return X


def load_downloaded(path, url, sha256, delimiter=","):
    """Load a public data file, downloading it if missing, and warn if its
    checksum differs from the file used to validate this script."""
    if not os.path.exists(path):
        print(f"Downloading {url} ...")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        urllib.request.urlretrieve(url, path)
    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
    if digest != sha256:
        print(f"  ! warning: {path} differs from the validated file "
              f"(sha256 {digest[:12]}...)")
    return np.loadtxt(path, delimiter=delimiter)


def betti_01(res):
    """beta_0 and beta_1 (coefficients Z/2) of our complex, with GUDHI.
    Exact from the 2-skeleton (beta_1 only needs simplices up to dim 2)."""
    import gudhi
    st = gudhi.SimplexTree()
    o = res["ordered"]
    for s, k, v in zip(o["simplices"], o["dims"], o["values"]):
        st.insert([int(i) for i in s[:k + 1]], filtration=float(v))
    st.compute_persistence(homology_coeff_field=2, persistence_dim_max=True)
    return tuple(st.betti_numbers()[:2])


def random_cube(seed, n=1000, dim=10):
    """Section 4.2: n points with i.i.d. coordinates uniform in [-1, 1]."""
    return np.random.default_rng(seed).uniform(-1.0, 1.0, size=(n, dim))


# ---------------------------------------------------------------------------
# One timed run (same protocol as run_experiments.py)
# ---------------------------------------------------------------------------
# --cpp: run the C++ program alphadual instead of alpha_dualqp.py. The time is
# then measured INSIDE the program (reading the input + computing the ordered
# filtration), see alphadual_io / alphadual.cpp.
USE_CPP = False
_LAST_CPP_OUTPUT = None
MAX_DIM = 2      # --max-dim: largest simplex dimension computed (paper: up to 3)


class Rows(list):
    """List of result rows that also appends every row to the CSV file at
    once (flushed), so the results survive a crash of the script."""
    def __init__(self, path):
        super().__init__()
        self.path = path
        self.fields = None
        open(path, "w").close()

    def append(self, row):
        super().append(row)
        row = {k: (";".join(map(str, v)) if isinstance(v, list) else v)
               for k, v in row.items()}
        with open(self.path, "a", newline="") as fh:
            if self.fields is None:
                self.fields = list(row)
                csv.DictWriter(fh, fieldnames=self.fields).writeheader()
            csv.DictWriter(fh, fieldnames=self.fields,
                           extrasaction="ignore").writerow(row)
            fh.flush()
            os.fsync(fh.fileno())


def compare_counts(x, p):
    """Our counts and the paper's, over the dimensions both have."""
    n = min(MAX_DIM, len(p["counts"]) - 1) + 1
    return tuple(x["counts"][:n]), tuple(p["counts"][:n])


def counts_str(c):
    return "/".join(str(v) for v in c)


def timed_run_cpp(X, r, solver):
    import tempfile
    import alphadual_io as io
    global _LAST_CPP_OUTPUT
    tmp = tempfile.mkdtemp()
    inp, out = os.path.join(tmp, "in.bin"), os.path.join(tmp, "out.alpha")
    io.write_input(inp, X)                                    # untimed
    st = io.run_alphadual(inp, out, r=r, solver=solver, max_dim=MAX_DIM)
    _LAST_CPP_OUTPUT = out
    c = st["n_simplices"]
    return dict(time_s=st["alpha_time_s"], read_time_s=st["read_time_s"],
                compute_time_s=st["compute_time_s"],
                counts=c, V=c[0], E=c[1] if len(c) > 1 else 0, T=c[2] if len(c) > 2 else 0,
                max_degree=st["max_cech_degree"], cech_edges=st["cech_edges"],
                qps=sum(st[f"dim{k}"]["candidates"] for k in range(len(c))),
                failed=sum(st[f"dim{k}"]["failed"] for k in range(len(c))))


def timed_run(X, r, solver):
    if USE_CPP:
        return timed_run_cpp(X, r, solver)
    t0 = time.perf_counter()
    res = aq.alpha_complex_dualqp(X, r=r, solver=solver, max_dim=MAX_DIM)  # timed call
    elapsed = time.perf_counter() - t0
    st = aq.compute_stats(res)                               # bookkeeping, untimed
    c = [len(res["simplices"][k]) for k in range(MAX_DIM + 1)]
    return dict(time_s=elapsed,
                counts=c, V=c[0], E=c[1] if len(c) > 1 else 0, T=c[2] if len(c) > 2 else 0,
                max_degree=st["max_cech_degree"], cech_edges=st["cech_edges"],
                qps=sum(st[f"dim{k}"]["candidates"] for k in range(MAX_DIM + 1)),
                failed=sum(st[f"dim{k}"]["failed"] for k in range(MAX_DIM + 1)))


def gudhi_counts(X, r):
    """Independent reference: GUDHI alpha complex counts for k = 0..3 at r."""
    import gudhi
    st = gudhi.AlphaComplex(points=X).create_simplex_tree(max_alpha_square=r * r)
    c = [0, 0, 0, 0]
    for s, f in st.get_skeleton(3):
        if f <= r * r:
            c[len(s) - 1] += 1
    return c


# ---------------------------------------------------------------------------
# Experiments
# ---------------------------------------------------------------------------
def experiment_4_1(solvers, repeats, check_gudhi, rows):
    print("\n=== Section 4.1: Stanford bunny, 5000 points in R^3, r = 0.01 ===")
    X = load_bunny()
    p = PAPER["4.1"]
    if check_gudhi:
        try:
            g = gudhi_counts(X, 0.01)
            print(f"GUDHI reference counts (k=0..3): {g}   paper: {list(p['counts'])}")
        except ImportError:
            print("  (GUDHI not installed: skipping the reference check)")
    for solver in solvers:
        # Same data every repeat: the spread only reflects timing noise.
        runs = [timed_run(X, 0.01, solver) for _ in range(repeats)]
        t = np.array([x["time_s"] for x in runs])
        x = runs[0]
        same = all(y["counts"] == x["counts"] for y in runs)
        ours, ref = compare_counts(x, p)
        fail_note = f"  [{x['failed']} QP failures]" if x["failed"] else ""
        print(f"[{solver:6s}] counts = {counts_str(x['counts'])} "
              f"({'= paper' if ours == ref else 'DIFFERS from paper ' + counts_str(ref)})"
              f"  max degree {x['max_degree']} (paper {p['max_degree']})"
              f"  time {t.mean():.2f} +- {t.std(ddof=1) if len(t) > 1 else 0:.2f} s"
              f" (paper ~{p['time_s']:.0f} s for dims 0-3)"
              f"{'' if same else '  [counts varied between repeats!]'}"
              f"{fail_note}")
        for i, y in enumerate(runs, 1):
            rows.append(dict(experiment="4.1 bunny", solver=solver, seed="", repeat=i, **y))


def experiment_4_2(solvers, repeats, seeds, rows):
    print(f"\n=== Section 4.2: 1000 uniform points in [-1,1]^10, r = 1.0, "
          f"{seeds} random samples ===")
    p = PAPER["4.2"]
    for solver in solvers:
        res = []
        for seed in range(seeds):
            X = random_cube(seed)                              # untimed
            for i in range(repeats if seed == 0 else 1):       # timing repeats on seed 0
                y = timed_run(X, 1.0, solver)
                rows.append(dict(experiment="4.2 cube10", solver=solver, seed=seed,
                                 repeat=i + 1, **y))
                if i == 0:
                    res.append(y)
            print(f"  [{solver}] seed {seed}: counts {counts_str(res[-1]['counts'])} "
                  f"max degree {res[-1]['max_degree']}  {res[-1]['time_s']:.2f} s",
                  flush=True)
        stats = [(k, name) for k, name in ((1, "edges"), (2, "triangles"), (3, "tetrahedra"))
                 if k <= MAX_DIM] + [(None, "max degree")]
        for k, name in stats:
            v = np.array([y["counts"][k] if k is not None else y["max_degree"] for y in res])
            ref = p["counts"][k] if k is not None else p["max_degree"]
            sd = v.std(ddof=1) if len(v) > 1 else 0.0
            # How many of our standard deviations the paper's value is from our
            # mean: |z| < 2 means the paper's single sample is typical of ours.
            z = (ref - v.mean()) / sd if sd > 0 else float("nan")
            verdict = ("consistent" if abs(z) < 2 else "check") if sd > 0 else "need >1 seed"
            print(f"[{solver:6s}] {name:10s}: mean {v.mean():9.0f}  std {sd:7.0f}  "
                  f"range [{v.min()}, {v.max()}]   paper {ref} = mean {z:+.1f} std"
                  f"  ({verdict})")
        t = np.array([y["time_s"] for y in res])
        print(f"[{solver:6s}] time      : mean {t.mean():.2f} s over {len(t)} samples "
              f"(paper ~{p['time_s']:.0f} s for dims 0-3, incl. 1.78M tetrahedra)")


def experiment_4_3(solvers, repeats, betti, rows):
    cases = [
        ("4.3a", "cyclooctane, 6040 points in R^24", CYCLO_FILE, CYCLO_URL,
         CYCLO_SHA256, 0.25),
        ("4.3b", "optical patches X(300,30), 1000 points (R^8 = DCT coords)",
         OPTICAL_FILE, OPTICAL_URL, OPTICAL_SHA256, 0.5),
    ]
    for key, title, path, url, sha, r in cases:
        print(f"\n=== Section {key[:3]}: {title}, r = {r} ===")
        X = load_downloaded(path, url, sha)
        p = PAPER[key]
        for solver in solvers:
            runs = []
            for i in range(repeats):
                y = timed_run(X, r, solver)
                runs.append(y)
                rows.append(dict(experiment=f"{key} {title.split(',')[0]}",
                                 solver=solver, seed="", repeat=i + 1, **y))
            t = np.array([y["time_s"] for y in runs])
            x = runs[0]
            ours, ref = compare_counts(x, p)
            diff = "= paper" if ours == ref else (
                "differs from paper " + counts_str(ref) + " by "
                + "/".join(f"{a - b:+d}" for a, b in zip(ours, ref)))
            print(f"[{solver:6s}] counts = {counts_str(x['counts'])} ({diff})"
                  f"  max degree {x['max_degree']}"
                  f"  time {t.mean():.2f} +- {t.std(ddof=1) if len(t) > 1 else 0:.2f} s"
                  f" (paper ~{p['time_s']:.0f} s, {p['time_note']})"
                  f"{'  [' + str(x['failed']) + ' QP failures]' if x['failed'] else ''}")
        if betti:
            try:
                if USE_CPP:
                    import alphadual_io as io
                    timed_run_cpp(X, r, solvers[0])
                    res = {"ordered": io.read_output(_LAST_CPP_OUTPUT)}
                else:
                    res = aq.alpha_complex_dualqp(X, r=r, solver=solvers[0])
                print(f"beta_0, beta_1 of our complex ({solvers[0]}): {betti_01(res)}"
                      f"   paper Betti numbers: {p['betti']}")
            except ImportError:
                print("  (GUDHI not installed: skipping Betti numbers)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", choices=["4.1", "4.2", "4.3"], help="run just one section")
    ap.add_argument("--solvers", default="native,daqp",
                    help="comma-separated: native, daqp (default: both)")
    ap.add_argument("--repeats", type=int, default=3,
                    help="timing repetitions (4.1; and seed 0 of 4.2)")
    ap.add_argument("--seeds", type=int, default=10, help="random samples for 4.2")
    ap.add_argument("--check-gudhi", action="store_true",
                    help="also count the bunny complex with GUDHI (k = 0..3)")
    ap.add_argument("--betti", action="store_true",
                    help="4.3: also compute beta_0, beta_1 of our complexes (needs GUDHI)")
    ap.add_argument("--max-dim", type=int, default=2,
                    help="largest simplex dimension (3 compares tetrahedra with the paper)")
    ap.add_argument("--cpp", action="store_true",
                    help="time the C++ program alphadual (must be built, same folder)")
    ap.add_argument("--csv", default="paper_reproduction_results.csv",
                    help="where to write every individual run (written as it finishes)")
    args = ap.parse_args()

    global USE_CPP, MAX_DIM
    USE_CPP = args.cpp
    MAX_DIM = args.max_dim
    solvers = [s.strip() for s in args.solvers.split(",") if s.strip()]
    for s in solvers:
        if s not in aq.SOLVERS:
            sys.exit(f"unknown solver {s}; choose from {aq.SOLVERS}")
        if s == "daqp" and not aq.HAVE_DAQP and not USE_CPP:
            sys.exit("DAQP is not installed (pip install daqp), or use --solvers native")

    if USE_CPP:
        import alphadual_io
        if not os.path.exists(alphadual_io.DEFAULT_EXE):
            sys.exit("alphadual not found next to this script: run make first")
        print("Timing the C++ program alphadual (time measured inside the program)")
    else:
        print("Warming up (Numba compilation, not timed)...")
        for s in solvers:
            aq.warmup(s)

    rows = Rows(args.csv)          # every run is written to the CSV immediately
    if args.only in (None, "4.1"):
        experiment_4_1(solvers, args.repeats, args.check_gudhi, rows)
    if args.only in (None, "4.2"):
        experiment_4_2(solvers, args.repeats, args.seeds, rows)
    if args.only in (None, "4.3"):
        experiment_4_3(solvers, args.repeats, args.betti, rows)

    print(f"\nEvery run written to {os.path.abspath(args.csv)}")


if __name__ == "__main__":
    main()
