"""
validate_against_gudhi.py
Checks alpha_dualqp (both solvers: DAQP and native) -- or, with --cpp, the C++
program alphadual -- against GUDHI's AlphaComplex (full Delaunay, so only
feasible in low ambient dimension): same simplices up to the chosen dimension and the
same filtration values (squared radii). Run:
    python validate_against_gudhi.py          # Python implementation
    python validate_against_gudhi.py --cpp    # C++ program (after make)
"""
import os
import sys
import numpy as np
import tadasets
from alpha_dualqp import alpha_complex_dualqp, compute_stats, warmup, SOLVERS, HAVE_DAQP

try:
    import gudhi
except ImportError:
    sys.exit("GUDHI is needed for validation:  pip install gudhi")


# Reference answer from GUDHI (which builds the full Delaunay triangulation
# with CGAL, then assigns alpha values). Returns {dim: {simplex: value}}
# for dims 0-2, keeping simplices with value <= a1. GUDHI's filtration
# values are squared radii, the same convention as ours.
def gudhi_simplices(X, a1, weights=None, max_dim=2):
    ac = (gudhi.AlphaComplex(points=X) if weights is None
          else gudhi.AlphaComplex(points=X, weights=list(weights)))
    st = (ac.create_simplex_tree() if not np.isfinite(a1)
          else ac.create_simplex_tree(max_alpha_square=a1))
    out = {k: {} for k in range(max_dim + 1)}
    for s, f in st.get_skeleton(max_dim):
        if f <= a1:
            out[len(s) - 1][tuple(sorted(s))] = f
    return out


# Our answer, in the same {dim: {simplex: value}} format.
# --cpp: check the C++ program alphadual instead of alpha_dualqp.py. Each case
# is written as an alphadual input file, run through the executable, and its
# ordered filtration is read back; the comparisons with GUDHI are unchanged.
USE_CPP = False


def ours_cpp(X, a1, weights, solver, max_dim):
    import tempfile
    import alphadual_io as io
    tmp = tempfile.mkdtemp()
    inp, out = os.path.join(tmp, "in.bin"), os.path.join(tmp, "out.alpha")
    io.write_input(inp, X, weights)
    stats = io.run_alphadual(inp, out, a1=a1, solver=solver, max_dim=max_dim)
    o = io.read_output(out)
    res = {"ordered": o, "stats": stats}
    per_dim = {k: {} for k in range(max_dim + 1)}
    for s, d, f in zip(o["simplices"], o["dims"], o["values"]):
        per_dim[int(d)][tuple(sorted(int(v) for v in s[:d + 1]))] = float(f)
    return per_dim, res


def ours(X, a1, weights=None, solver="daqp", max_dim=2):
    if USE_CPP:
        return ours_cpp(X, a1, weights, solver, max_dim)
    res = alpha_complex_dualqp(X, a1=a1, weights=weights, solver=solver, max_dim=max_dim)
    out = {}
    for k in range(max_dim + 1):
        out[k] = {tuple(sorted(map(int, s))): float(f)
                  for s, f in zip(res["simplices"][k], res["filtration"][k])}
    return out, res


# The ordered output must be a valid filtration: values never decrease
# along the list, and each simplex appears after all of its faces.
def check_order(res):
    """Ordered output: values non-decreasing and every face before its cofaces."""
    o = res["ordered"]
    if np.any(np.diff(o["values"]) < 0):
        return False
    # Walk the list, remembering what has appeared; for every edge/triangle,
    # check that all its faces (drop one vertex at a time) were seen before.
    seen = set()
    for s, d in zip(o["simplices"], o["dims"]):
        s = tuple(int(i) for i in s[:d + 1])
        if d >= 1 and any(f not in seen for f in
                          (tuple(s[:i] + s[i + 1:]) for i in range(d + 1))):
            return False
        seen.add(s)
    return True


# Compare ours vs GUDHI on one point cloud and print a one-line report.
# A case passes when, in every dimension: no simplex is missing, none is
# extra, filtration values agree to 1e-8 (relative to the largest value),
# the ordered output is a valid filtration, and no QP failed.
def compare(name, X, a1, weights=None, max_dim=2, solver="daqp"):
    g = gudhi_simplices(X, a1, weights, max_dim)
    o, res = ours(X, a1, weights, solver, max_dim)
    ok = True
    msg = []
    for k in range(max_dim + 1):
        gs, os_ = set(g[k]), set(o[k])
        miss, extra = gs - os_, os_ - gs
        common = gs & os_
        err = max((abs(g[k][s] - o[k][s]) for s in common), default=0.0)
        # Relative tolerance: the complete filtration has values up to ~1e5 on
        # the convex hull, where absolute agreement to 1e-8 is not meaningful.
        tol = 1e-8 * max(1.0, max(g[k].values(), default=1.0))
        if miss or extra or err > tol:
            ok = False
        msg.append(f"dim{k}: {len(os_)} (gudhi {len(gs)}), missing {len(miss)},"
                   f" extra {len(extra)}, max|dfilt| {err:.1e}")
    stats = res["stats"] if USE_CPP else compute_stats(res)
    failed = sum(stats[f"dim{k}"]["failed"] for k in range(max_dim + 1))
    order_ok = check_order(res)
    ok = ok and order_ok
    msg.append(f"order {'valid' if order_ok else 'INVALID'}")
    impl = "C++ " if USE_CPP else ""
    print(f"[{'OK ' if ok and not failed else 'BAD'}] {impl}{solver:6s} d={max_dim} {name}: "
          + "; ".join(msg) + f"; QP failures {failed}")
    return ok and not failed


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cpp", action="store_true",
                    help="validate the C++ program alphadual (built with make, same "
                         "folder) instead of the Python implementation")
    USE_CPP = ap.parse_args().cpp
    if USE_CPP:
        import alphadual_io
        exe = alphadual_io.DEFAULT_EXE         # alphadual (alphadual.exe on Windows)
        if not os.path.exists(exe):
            sys.exit("alphadual not found next to this script: run make first")
        print("Validating the C++ program alphadual against GUDHI")
    else:
        print("Validating the Python implementation (alpha_dualqp.py) against GUDHI")
        # Compile first, so the test output is not delayed by Numba.
        warmup()
    rng = np.random.default_rng(1)
    # Test cases: (name, points, a1 = r^2 or inf, weights). They cover 2D-4D,
    # several tadasets shapes, a weighted case, and the complete filtration.
    # (name, points, a1 = r^2 or inf, weights, max_dim)
    tests = [
        ("uniform 2D, n=400", rng.random((400, 2)), 0.05 ** 2, None, 2),
        ("uniform 3D, n=500", rng.random((500, 3)), 0.12 ** 2, None, 2),
        ("uniform 4D, n=300", rng.random((300, 4)), 0.25 ** 2, None, 2),
        ("torus 3D, n=800", tadasets.torus(800, noise=0.02, seed=3), 0.3 ** 2, None, 2),
        ("sphere 3D, n=600", tadasets.sphere(600, seed=4), 0.2 ** 2, None, 2),
        ("swiss roll 3D, n=600", tadasets.swiss_roll(600, seed=5), 1.5 ** 2, None, 2),
        ("infty sign 2D, n=300", tadasets.infty_sign(300, noise=0.01, seed=6), 0.1 ** 2, None, 2),
        ("weighted 3D, n=300", rng.random((300, 3)), 0.01, rng.random(300) * 0.01, 2),
        # complete filtration (all radii): r = inf
        ("ALL RADII uniform 2D, n=400", rng.random((400, 2)), np.inf, None, 2),
        ("ALL RADII torus 3D, n=500", tadasets.torus(500, noise=0.02, seed=7), np.inf, None, 2),
        ("ALL RADII uniform 4D, n=200", rng.random((200, 4)), np.inf, None, 2),
        # higher dimensions (and the trivial d = 0, 1; d = 5 > ambient dim 4)
        ("uniform 3D, n=500", rng.random((500, 3)), 0.12 ** 2, None, 3),
        ("torus 3D, n=800", tadasets.torus(800, noise=0.02, seed=8), 0.3 ** 2, None, 3),
        ("uniform 4D, n=300", rng.random((300, 4)), 0.25 ** 2, None, 4),
        ("weighted 3D, n=300", rng.random((300, 3)), 0.01, rng.random(300) * 0.01, 3),
        ("ALL RADII uniform 3D, n=200", rng.random((200, 3)), np.inf, None, 3),
        ("ALL RADII uniform 4D, n=120", rng.random((120, 4)), np.inf, None, 4),
        ("uniform 4D, n=200", rng.random((200, 4)), 0.3 ** 2, None, 5),
        ("uniform 3D, n=300", rng.random((300, 3)), 0.12 ** 2, None, 0),
        ("uniform 3D, n=300", rng.random((300, 3)), 0.12 ** 2, None, 1),
    ]
    # (DAQP is built into the C++ program; Python needs the daqp package)
    solvers = SOLVERS if (HAVE_DAQP or USE_CPP) else ("native",)
    if not HAVE_DAQP and not USE_CPP:
        print("DAQP not installed (pip install daqp): testing the native solver only")
    # the complete filtration (r = inf) is only supported with DAQP
    results = [compare(*t, solver=s) for s in solvers for t in tests
               if s == "daqp" or np.isfinite(t[2])]
    print(f"\n{sum(results)}/{len(results)} test cases agree with GUDHI"
          f" ({'C++ alphadual' if USE_CPP else 'Python implementation'}).")
    sys.exit(0 if all(results) else 1)
