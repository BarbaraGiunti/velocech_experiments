"""
compare_cpp_python.py
Checks that the C++ program `alphadual` produces the same ordered alpha
filtration as the validated Python implementation (alpha_dualqp.py), for both
solvers, on a range of inputs (low/high dimension, weighted, complete
filtration, and the data sets of Section 4 of the paper).
Run from the folder containing alphadual, alpha_dualqp.py and the data files:
    python compare_cpp_python.py
"""
import os
import sys
import tempfile

import numpy as np
import tadasets

import alpha_dualqp as aq
import alphadual_io as io

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def as_dict(o):
    """{(dim, v0, ..): value} for an ordered filtration."""
    return {(int(k),) + tuple(int(v) for v in s[:k + 1]): float(f)
            for s, k, f in zip(o["simplices"], o["dims"], o["values"])}


def order_valid(o):
    """Values non-decreasing and every face listed before its cofaces."""
    if np.any(np.diff(o["values"]) < 0):
        return False
    seen = set()
    for s, d in zip(o["simplices"], o["dims"]):
        s = tuple(int(i) for i in s[:d + 1])
        if d >= 1 and any(tuple(s[:i] + s[i + 1:]) not in seen for i in range(d + 1)):
            return False
        seen.add(s)
    return True


# Pass criterion: identical simplex sets, filtration values within VALUE_TOL
# (relative), and a valid filtration order. The tolerance is 1e-6 because on
# nearly degenerate data (the optical patches of Section 4.3) the native
# solver's values are only accurate to ~1e-7 on sliver faces -- in Python as
# in C++ (Python's two solvers differ by ~1e-6 there). Elsewhere agreement is
# ~1e-12 or better. Exact ties may be ordered differently after rounding, so
# the order itself is checked for validity rather than for identity.
VALUE_TOL = 1e-6


def compare(name, X, r=None, a1=None, w=None, solvers=("native", "daqp"), max_dim=2):
    tmp = tempfile.mkdtemp()
    inp = os.path.join(tmp, "in.bin")
    out = os.path.join(tmp, "out.bin")
    io.write_input(inp, X, w)
    ok_all = True
    for s in solvers:
        if r is not None and np.isinf(r) and s == "native":
            continue                      # complete filtration needs DAQP
        py = aq.alpha_complex_dualqp(X, r=r, a1=a1, weights=w, solver=s,
                                     max_dim=max_dim)["ordered"]
        st = io.run_alphadual(inp, out, r=r, a1=a1, solver=s, max_dim=max_dim)
        cp = io.read_output(out)
        P, C = as_dict(py), as_dict(cp)
        miss, extra = set(P) - set(C), set(C) - set(P)
        common = set(P) & set(C)
        err = max((abs(P[k] - C[k]) / max(1.0, abs(P[k])) for k in common), default=0.0)
        valid = order_valid(cp)
        ok = not miss and not extra and err < VALUE_TOL and valid
        ok_all &= ok
        print(f"[{'OK ' if ok else 'BAD'}] {name:32s} {s:6s} d={max_dim}: {len(C):8d} simplices, "
              f"missing {len(miss)}, extra {len(extra)}, max rel |dvalue| {err:.1e}, "
              f"order {'valid' if valid else 'INVALID'}, C++ {st['alpha_time_s']:.3f}s "
              f"(read {st['read_time_s'] * 1e3:.1f} ms)", flush=True)
    return ok_all


if __name__ == "__main__":
    aq.warmup("native")
    if aq.HAVE_DAQP:
        aq.warmup("daqp")
    rng = np.random.default_rng(3)
    np.random.seed(0)                     # tadasets.embed uses the global RNG
    cases = [
        ("uniform 2D n=400", rng.random((400, 2)), 0.05, None, None),
        ("uniform 4D n=300", rng.random((300, 4)), 0.25, None, None),
        ("torus 3D n=2000", tadasets.torus(2000, noise=0.02, seed=1), 0.4, None, None),
        ("torus -> R^100 n=1000", tadasets.embed(tadasets.torus(1000, seed=2), 100),
         0.4, None, None),
        ("weighted 3D n=300", rng.random((300, 3)), None, 0.01, rng.random(300) * 0.01),
        ("ALL RADII torus n=300", tadasets.torus(300, noise=0.02, seed=5), np.inf, None, None),
    ]
    # higher dimensions: (name, points, r, a1, weights, max_dim)
    high = [
        ("uniform 3D n=600", rng.random((600, 3)), 0.12, None, None, 3),
        ("uniform 4D n=300", rng.random((300, 4)), 0.25, None, None, 4),
        ("torus -> R^20 n=800", tadasets.embed(tadasets.torus(800, seed=3), 20), 0.4, None, None, 3),
        ("weighted 3D n=300", rng.random((300, 3)), None, 0.01, rng.random(300) * 0.01, 3),
        ("ALL RADII uniform 4D n=120", rng.random((120, 4)), np.inf, None, None, 4),
        ("uniform 3D n=300 d=0", rng.random((300, 3)), 0.12, None, None, 0),
        ("uniform 3D n=300 d=1", rng.random((300, 3)), 0.12, None, None, 1),
    ]
    data = [("bunny (paper 4.1) r=0.01", os.path.join(DATA, "bunny_5000.txt"),
             lambda p: np.loadtxt(p, skiprows=1), 0.01),
            ("cube10 (paper 4.2) r=1", None,
             lambda p: np.random.default_rng(0).uniform(-1, 1, (1000, 10)), 1.0),
            ("cyclooctane (paper 4.3) r=0.25", os.path.join(DATA, "cyclooctane_6040.csv"),
             lambda p: np.loadtxt(p, delimiter=","), 0.25),
            ("optical (paper 4.3) r=0.5", os.path.join(DATA, "optical_k300_points.txt"),
             lambda p: np.loadtxt(p, delimiter=","), 0.5)]
    for name, path, load, r in data:
        if path is None or os.path.exists(path):
            cases.append((name, load(path), r, None, None))
        else:
            print(f"(skipping {name}: {path} not found)")
    solvers = ("native", "daqp") if aq.HAVE_DAQP else ("native",)
    res = [compare(n, X, r, a1, w, solvers) for n, X, r, a1, w in cases]
    res += [compare(n, X, r, a1, w, solvers, d) for n, X, r, a1, w, d in high]
    print(f"\n{sum(res)}/{len(res)} cases: C++ agrees with Python")
    sys.exit(0 if all(res) else 1)
