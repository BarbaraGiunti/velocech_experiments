"""
alpha_at_half_diameter.py
=========================

For each given point cloud: compute its diameter D (the largest distance
between two points), then the alpha filtration up to radius r = D / 2, with
the implementation of this folder (C++ program alphadual by default, or the
Python implementation). The TIME MEASURED is the diameter computation plus
the alpha construction -- reading the input is NOT included.

    python alpha_at_half_diameter.py myrun/points                 # all point clouds of a folder
    python alpha_at_half_diameter.py cloud1.csv cloud2.csv        # given files
    python alpha_at_half_diameter.py myrun/points --max-dim 3 --repeats 3
    python alpha_at_half_diameter.py data/bunny_5000.txt --impl python

Options:
  --impl cpp|python     implementation (default cpp: the program alphadual,
                        which must be built, e.g. by running run_benchmark.sh once)
  --solver daqp|native  dual QP solver (default daqp)
  --max-dim D           largest simplex dimension (default 2: up to triangles)
  --repeats K           runs per input (default 1)
  --output NAME         results file name, without extension (default
                        alpha_half_diameter_results): NAME.csv and NAME.xlsx

Inputs: .csv (comma-separated, one point per line: the point clouds saved by
the benchmark in points/), .txt / .xyz (comma- or space-separated; a first line
holding only the number of points, as in data/bunny_5000.txt, is skipped),
.bin (alphadual input files, e.g. inputs/ of a benchmark), .npy. Files named
index.csv are skipped. Folders are searched (not recursively) in name order.

TIMING
  * The input is read first (not timed).
  * --impl cpp: the clock measures the diameter computation (in this
    process); then the cloud is passed to alphadual with r = D/2, and its own
    measurement of the computation alone (compute_time_s: from the end of
    reading its input to the finished ordered filtration) is added. Both parts
    exclude reading, so total_time_s = diameter_time_s + alpha_time_s is
    exactly "diameter + alpha construction". (Two intervals, because the
    existing C++ program, which is not modified, cannot compute the diameter.)
  * --impl python: ONE interval from the start of the diameter computation to
    the finished filtration, in this process; one extra clock reading at the
    end of the diameter computation splits it into the two columns.
  * The diameter routine is warmed up once before the first input (its
    one-time start-up cost is not measured).
  * The diameter is exact and computed single-threaded (like the alpha
    construction): for points in dimension <= 3 via the convex hull (qhull)
    and the largest distance between hull vertices; otherwise all pairwise
    distances, block by block (scipy, compiled C).

COST: at r = D/2 every two points are at most 2r apart, so the Cech graph is
complete: every pair of points is a candidate edge and every point constrains
every QP, as for the complete filtration -- roughly cubic in the number of
points, and memory ~ N^2 per vertex. Practical up to a few thousand points.
The result is (almost) the whole Delaunay complex: only simplices whose entry
radius exceeds D/2 (flat simplices on the hull) are left out.

Results are written after every run (CSV appended, Excel rewritten), so they
survive a crash. Needs the other files of this folder (it imports
alphadual_io, alpha_dualqp and run_experiments) but does not modify them.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import sys
import tempfile
import time

import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import alphadual_io as io   # noqa: E402  (C++ program: input/output files)

EXTENSIONS = (".csv", ".txt", ".xyz", ".bin", ".npy")


# ---------------------------------------------------------------------------
# Reading (never timed)
# ---------------------------------------------------------------------------
def read_points(path):
    """N x m float64 array from a .csv/.txt/.xyz/.bin/.npy file."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".npy":
        X = np.load(path)
    elif ext == ".bin":                          # alphadual input format
        buf = open(path, "rb").read()
        if buf[:8] != b"ALPHAIN1":
            raise ValueError("not an alphadual input file")
        N, m, _ = (int(v) for v in np.frombuffer(buf, "<i8", 3, 8))
        X = np.frombuffer(buf, "<f8", N * m, 32).reshape(N, m)
    else:                                        # text: comma or whitespace
        rows = []
        with open(path) as fh:
            for line in fh:
                parts = line.replace(",", " ").split()
                if parts:
                    rows.append([float(v) for v in parts])
        if len(rows) > 1 and len(rows[0]) == 1 and len(rows[1]) > 1:
            rows = rows[1:]                      # header line with the point count
        X = np.array(rows, dtype=np.float64)
    X = np.ascontiguousarray(X, dtype=np.float64)
    if X.ndim != 2 or len(X) < 2:
        raise ValueError(f"expected at least 2 points in an N x m array, got shape {X.shape}")
    return X


def expand_inputs(paths):
    files = []
    for p in paths:
        if os.path.isdir(p):
            files += [os.path.join(p, f) for f in sorted(os.listdir(p))
                      if f.lower().endswith(EXTENSIONS) and f.lower() != "index.csv"]
        elif os.path.isfile(p):
            files.append(p)
        else:
            sys.exit(f"not found: {p}")
    if not files:
        sys.exit("no input files found")
    return files


# ---------------------------------------------------------------------------
# Diameter (timed): exact, single-threaded
# ---------------------------------------------------------------------------
def diameter(X, block=2048):
    """Largest distance between two points. Dimension <= 3: only convex-hull
    vertices can realise it, so compute the hull (qhull) and check those
    pairs; otherwise (or if the hull fails, e.g. all points on a line), all
    pairs in blocks with scipy's compiled cdist."""
    N, m = X.shape
    P = X
    if m <= 3 and N > 64:
        try:
            P = X[ConvexHull(X).vertices]
        except Exception:  # noqa: BLE001  (degenerate: flat / collinear points)
            P = X
    best = 0.0
    for s in range(0, len(P), block):
        D = cdist(P[s:s + block], P[s:], "sqeuclidean")
        best = max(best, float(D.max()))
    return float(np.sqrt(best))


# ---------------------------------------------------------------------------
# One run
# ---------------------------------------------------------------------------
def run_cpp(X, solver, max_dim, tmpdir):
    t0 = time.perf_counter()
    D = diameter(X)                                          # timed (this process)
    t_diam = time.perf_counter() - t0
    inp = os.path.join(tmpdir, "in.bin")
    out = os.path.join(tmpdir, "out.alpha")
    io.write_input(inp, X)                                   # not timed
    st = io.run_alphadual(inp, out, r=D / 2, solver=solver, max_dim=max_dim)
    counts = st["n_simplices"]
    failed = sum(st[f"dim{k}"]["failed"] for k in range(len(counts)))
    # alphadual's own measurement of the computation alone (reading excluded)
    return D, t_diam, st["compute_time_s"], counts, st["max_cech_degree"], failed


def run_python(X, solver, max_dim, aq):
    t0 = time.perf_counter()                                 # ONE interval:
    D = diameter(X)
    t1 = time.perf_counter()                                 # (split point only)
    res = aq.alpha_complex_dualqp(X, r=D / 2, solver=solver, max_dim=max_dim)
    t2 = time.perf_counter()
    st = aq.compute_stats(res)                               # not timed
    counts = [len(res["simplices"][k]) for k in range(max_dim + 1)]
    failed = sum(st[f"dim{k}"]["failed"] for k in range(max_dim + 1))
    return D, t1 - t0, t2 - t1, counts, st["max_cech_degree"], failed


# ---------------------------------------------------------------------------
# Output: CSV (appended run by run) and Excel (rewritten after every run)
# ---------------------------------------------------------------------------
COLS = ["file", "n_points", "ambient_dim", "diameter", "r", "impl", "solver", "max_dim",
        "repeat", "status", "diameter_time_s", "alpha_time_s", "total_time_s",
        "max_cech_degree", "failed_QPs"]


def write_excel(path, rows, max_dim, info):
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter as col
    F, FB = Font(name="Arial", size=10), Font(name="Arial", size=10, bold=True, color="FFFFFF")
    FILL = PatternFill("solid", fgColor="1F4E78")
    cols = COLS + [f"n_simplices_dim{k}" for k in range(max_dim + 1)]
    wb = Workbook()
    wb.calculation.fullCalcOnLoad = True
    ws = wb.active
    ws.title = "Runs"
    for j, h in enumerate(cols, 1):
        c = ws.cell(row=1, column=j, value=h)
        c.font, c.fill = FB, FILL
        c.alignment = Alignment(horizontal="center", wrap_text=True)
        ws.column_dimensions[col(j)].width = max(11, min(28, len(h) + 2))
    ws.column_dimensions["A"].width = 30
    L = {h: col(j) for j, h in enumerate(cols, 1)}
    for i, r in enumerate(rows, 2):
        for h in cols:
            v = r.get(h)
            if h == "total_time_s" and r.get("status") == "ok":
                v = f"={L['diameter_time_s']}{i}+{L['alpha_time_s']}{i}"   # sum of the parts
            c = ws.cell(row=i, column=cols.index(h) + 1, value=v)
            c.font = F
            if h.endswith("_s"):
                c.number_format = "0.000000"
            elif h in ("diameter", "r"):
                c.number_format = "0.000000"
    ws.freeze_panes = "B2"

    # Summary: per input file, mean and sample std over the successful repeats
    sm = wb.create_sheet("Summary")
    heads = ["file", "runs_ok"] + [f"{s}_{t}" for t in ("diameter_time_s", "alpha_time_s",
                                                        "total_time_s") for s in ("mean", "std")]
    for j, h in enumerate(heads, 1):
        c = sm.cell(row=1, column=j, value=h)
        c.font, c.fill = FB, FILL
        sm.column_dimensions[col(j)].width = 30 if j == 1 else 20
    last = len(rows) + 1
    rng = lambda h: f"Runs!${L[h]}$2:${L[h]}${last}"  # noqa: E731
    files = list(dict.fromkeys(r["file"] for r in rows))
    for i, f in enumerate(files, 2):
        vals = [f, f"=COUNTIFS({rng('file')},$A{i},{rng('status')},\"ok\")"]
        j = 3
        for h in ("diameter_time_s", "alpha_time_s", "total_time_s"):
            mean = f"{col(j)}{i}"
            vals.append(f"=IFERROR(AVERAGEIFS({rng(h)},{rng('file')},$A{i},"
                        f"{rng('status')},\"ok\"),\"\")")
            vals.append(f"=IF($B{i}>1,SQRT(SUMPRODUCT(({rng('file')}=$A{i})*"
                        f"({rng('status')}=\"ok\")*({rng(h)}-{mean})^2)/($B{i}-1)),\"n/a\")")
            j += 2
        for jj, v in enumerate(vals, 1):
            c = sm.cell(row=i, column=jj, value=v)
            c.font = F
            if jj > 2:
                c.number_format = "0.000000"
    note = sm.cell(row=len(files) + 3, column=1, value=info)
    note.font = Font(name="Arial", size=9, italic=True)

    # atomic save, tolerant of a file open in Excel (shared helper)
    import run_experiments as rx
    tmp = os.path.join(os.path.dirname(os.path.abspath(path)), "." + os.path.basename(path) + ".tmp")
    wb.save(tmp)
    return rx.replace_locked_aware(tmp, path)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="point-cloud files and/or folders")
    ap.add_argument("--impl", choices=["cpp", "python"], default="cpp")
    ap.add_argument("--solver", choices=["daqp", "native"], default="daqp")
    ap.add_argument("--max-dim", type=int, default=2)
    ap.add_argument("--repeats", type=int, default=1)
    ap.add_argument("--output", default="alpha_half_diameter_results")
    args = ap.parse_args()
    if args.max_dim < 0 or args.repeats < 1:
        sys.exit("--max-dim must be >= 0 and --repeats >= 1")

    files = expand_inputs(args.inputs)
    aq = None
    if args.impl == "cpp":
        if not os.path.exists(io.DEFAULT_EXE):
            sys.exit(f"{io.DEFAULT_EXE} not found: build the C++ program first "
                     "(e.g. run ./run_benchmark.sh once, or make / cmake as in README.md)")
        version = io.subprocess.run([io.DEFAULT_EXE, "--version"], capture_output=True,
                                    text=True).stdout.strip()
    else:
        import alpha_dualqp as aq
        print("Warming up (Numba compilation, not timed)...")
        aq.warmup(args.solver)
        version = "Python implementation (alpha_dualqp.py)"
    info = (f"r = diameter/2; time = diameter computation + alpha construction "
            f"(reading excluded); impl {args.impl}, solver {args.solver}, "
            f"max_dim {args.max_dim}; {version}; "
            f"started {dt.datetime.now():%Y-%m-%d %H:%M:%S}")

    # Warm up the diameter routine (first-call costs of qhull / scipy are not
    # part of any measurement): both branches, hull and all-pairs.
    _rng = np.random.default_rng(0)
    diameter(_rng.random((200, 3)))
    diameter(_rng.random((200, 5)))

    csv_path, xlsx_path = args.output + ".csv", args.output + ".xlsx"
    with open(csv_path, "w", newline="") as fh:
        csv.writer(fh, lineterminator="\n").writerow(
            COLS + [f"n_simplices_dim{k}" for k in range(args.max_dim + 1)])
    rows = []
    tmpdir = tempfile.mkdtemp()
    print(f"{len(files)} input file(s), {args.repeats} run(s) each, impl {args.impl}, "
          f"solver {args.solver}, up to dimension {args.max_dim}\n")
    for f in files:
        try:
            X = read_points(f)                               # not timed
        except Exception as e:  # noqa: BLE001
            row = dict(file=os.path.basename(f), status=f"error reading: {e}")
            rows.append(row)
            print(f"{os.path.basename(f)}: {row['status']}")
            continue
        for k in range(args.repeats):
            row = dict(file=os.path.basename(f), n_points=X.shape[0], ambient_dim=X.shape[1],
                       impl=args.impl, solver=args.solver, max_dim=args.max_dim, repeat=k + 1)
            try:
                if args.impl == "cpp":
                    D, td, ta, counts, deg, failed = run_cpp(X, args.solver, args.max_dim, tmpdir)
                else:
                    D, td, ta, counts, deg, failed = run_python(X, args.solver, args.max_dim, aq)
                row.update(status="ok", diameter=D, r=D / 2, diameter_time_s=td,
                           alpha_time_s=ta, total_time_s=td + ta, max_cech_degree=deg,
                           failed_QPs=failed)
                row.update({f"n_simplices_dim{j}": c for j, c in enumerate(counts)})
                msg = (f"D = {D:.6g}, r = {D / 2:.6g}: diameter {td:.4f} s + alpha {ta:.4f} s "
                       f"= {td + ta:.4f} s, simplices per dim "
                       + "/".join(str(c) for c in counts)
                       + (f"  [{failed} QP failures]" if failed else ""))
            except Exception as e:  # noqa: BLE001  (e.g. out of memory in alphadual)
                row["status"] = f"error: {e}"
                msg = row["status"]
            rows.append(row)
            print(f"{os.path.basename(f)} (run {k + 1}): {msg}", flush=True)
            with open(csv_path, "a", newline="") as fh:          # saved at once
                csv.writer(fh, lineterminator="\n").writerow(
                    [row.get(h, "") for h in COLS]
                    + [row.get(f"n_simplices_dim{j}", "") for j in range(args.max_dim + 1)])
                fh.flush()
                os.fsync(fh.fileno())
            try:
                write_excel(xlsx_path, rows, args.max_dim, info)
            except Exception as e:  # noqa: BLE001
                print(f"  (warning: Excel file not updated this time: {e}; results are in {csv_path})")
    print(f"\nResults: {os.path.abspath(xlsx_path)}  and  {os.path.abspath(csv_path)}")


if __name__ == "__main__":
    main()
