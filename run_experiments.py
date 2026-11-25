"""
run_experiments.py
==================

Benchmark of the dual active-set QP alpha complex (vertices, edges, triangles;
Carlsson & Carlsson, arXiv:2310.00536) on point clouds from scikit-tda/tadasets.

    python run_experiments.py                 # interactive
    python run_experiments.py --config X.json # re-run a saved configuration

The script asks for: the tadasets shape and its parameters, the number(s) of
points, the ambient dimension(s), the maximum radius r_max ('inf' = all radii), the number of runs
per configuration, a base seed and the dual QP solver ("daqp": DAQP, the
solver of ref. [2] used by the paper's authors, default; "native": Numba
re-implementation of the same method). Lists may be given as "100,200,500" or as
ranges "100:1000:100" (end inclusive), or a mix of both.

Output: an Excel workbook with sheets
  Inputs  - every input, software/hardware info and the timing protocol
  Runs    - one row per run: alpha complex time, complex sizes, QP statistics
  Summary - per (n, ambient dim, r): number of runs, mean and sample standard
            deviation of every runtime (live Excel formulas over 'Runs')
  Chart   - mean total runtime vs. number of points (or vs. ambient dim)
and a JSON file with the configuration, re-usable with --config.

Output of each run: the alpha filtration up to r_max, i.e. every vertex, edge
and triangle entering at radius <= r_max, ordered by entry radius (optionally
saved to .npz).

Point clouds: every generated cloud is saved as CSV in <output>_points/,
named <shape>_<points>_<dimension>_<run>.csv (e.g. torus_1000_3_1.csv; see
index.csv there), to run other implementations on identical inputs.

Crash safety: every run is recorded in <output>_progress.jsonl as soon as it
finishes, and the Excel file is rewritten after every run. To continue an
interrupted experiment (crash, out of memory, Ctrl-C):
    python run_experiments.py --config <output>_config.json --resume
Finished runs are not repeated; a run that killed the process is recorded as
an error (--retry-crashed runs it again); a run stopped with Ctrl-C is re-done.

Timing protocol: a single timer around a single call, from the generated
point cloud to the ordered filtration of vertices, edges and triangles
(Algorithm 1: Cech graph + dims 0,1,2). There are no timers inside the
computation. Point generation, statistics and Excel writing are not timed, and
Numba compilation is done in a warm-up beforehand. Run r uses seed = base_seed + r for the shape sample and for the
random rotation used by tadasets.embed, so for a given run index the same
intrinsic cloud is used for every ambient dimension.
"""

from __future__ import annotations

import argparse
import ast
import datetime as dt
import inspect
import json
import os
import platform
import sys
import time

import numpy as np

try:
    import tadasets
except ImportError:
    sys.exit("tadasets is required:  pip install tadasets")

# The algorithm itself lives in alpha_dualqp.py (same folder).
import alpha_dualqp as aq

# tadasets arguments the script sets itself (so it does not ask for them):
# n = number of points, ambient = embedding dimension (we embed ourselves,
# see generate()), seed = random seed (derived from the run index).
EXCLUDED_ARGS = {"n", "ambient", "seed"}


# ---------------------------------------------------------------------------
# Input helpers
# ---------------------------------------------------------------------------
# List the shape generators offered by tadasets (torus, dsphere, ...).
# Read from the library at run time, so new tadasets versions just work.
def available_shapes():
    names = list(getattr(tadasets.shapes, "__all__", []))
    if not names:
        names = [k for k, v in vars(tadasets).items()
                 if inspect.isfunction(v) and "n" in inspect.signature(v).parameters
                 and k not in ("embed", "from_mesh", "plot3d")]
    return names


# The parameters of a tadasets function and their defaults, read from its
# signature (e.g. torus -> c, a, noise, uniform), minus EXCLUDED_ARGS.
def shape_params(fn):
    out = {}
    for name, p in inspect.signature(fn).parameters.items():
        if name in EXCLUDED_ARGS or p.kind in (p.VAR_POSITIONAL, p.VAR_KEYWORD):
            continue
        out[name] = None if p.default is inspect.Parameter.empty else p.default
    return out


# Ask one question on the terminal; Enter (or end of piped input) returns
# the default shown in [brackets].
def ask(prompt, default=None):
    d = "" if default is None else f" [{default}]"
    # print the question ourselves and flush: in some terminals (e.g. Git Bash
    # on Windows) Python's output is buffered and the question would not show
    print(f"{prompt}{d}: ", end="", flush=True)
    try:
        s = input().strip()
    except EOFError:
        s = ""
    return s if s else ("" if default is None else str(default))


# Turn a typed answer into a Python value: "0.1" -> 0.1, "True" -> True,
# "None" -> None; anything else stays a string.
def parse_value(s):
    if s.lower() == "none":
        return None
    try:
        return ast.literal_eval(s)
    except (ValueError, SyntaxError):
        return s


# Parse lists like "100,200", ranges like "100:1000:100" (end included),
# or mixes such as "3,10:50:10,100". 'inf' is accepted for floats (r_max).
def parse_list(s, cast):
    """'100,200' / '100:1000:100' (inclusive) / mixes -> sorted unique list."""
    vals = []
    for tok in s.replace(";", ",").replace(" ", ",").split(","):
        if not tok:
            continue
        if ":" in tok:
            parts = [float(x) for x in tok.split(":")]
            if len(parts) == 2:
                parts.append(1.0)
            a, b, st = parts
            if st <= 0:
                raise ValueError("range step must be positive")
            k = 0
            while a + k * st <= b + 1e-9 * abs(st):
                vals.append(cast(round(a + k * st, 12)))
                k += 1
        else:
            vals.append(cast(float(tok)))
    if not vals:
        raise ValueError("empty list")
    return sorted(set(vals))


# Small helper used by the input checks below.
def require(cond, msg):
    if not cond:
        raise ValueError(msg)


# Ask for a list (see parse_list) and re-ask until it parses and passes
# the optional check.
def ask_list(prompt, cast, default, check=None):
    while True:
        s = ask(prompt, default)
        try:
            vals = parse_list(s, cast)
            if check:
                check(vals)
            return vals
        except Exception as e:  # noqa: BLE001
            print(f"  ! could not use '{s}': {e}")


# Build one point cloud. Why not just pass ambient=... to tadasets?
# tadasets.embed rotates the data with numpy's GLOBAL random generator, not
# with `seed`, so the rotation would not be reproducible. We therefore
# sample in the native dimension and embed ourselves after seeding the
# global generator. The result is the same as tadasets' own ambient option
# (noise is added before embedding in both cases).
def generate(fn, params, n, ambient, native, seed):
    """tadasets sample in its native dimension, then (if ambient > native)
    embedded with tadasets.embed under a seeded global RNG."""
    sig = inspect.signature(fn).parameters
    kwargs = dict(params)
    kwargs["n"] = int(n)
    np.random.seed(seed)
    if "seed" in sig:
        kwargs["seed"] = int(seed)
    X = np.asarray(fn(**kwargs), dtype=np.float64)
    if ambient > native:
        np.random.seed(seed)
        X = tadasets.embed(X, ambient)
    return X


# Typical length scales of a sample (median nearest-neighbour distance and
# approximate diameter), printed to help choose r_max.
def scale_hint(fn, params, n, native, seed):
    X = generate(fn, params, n, native, native, seed)
    from scipy.spatial import cKDTree
    d, _ = cKDTree(X).query(X, k=2)
    sub = X[np.random.default_rng(0).choice(len(X), min(len(X), 1500), replace=False)]
    diam = np.sqrt(((sub[:, None] - sub[None]) ** 2).sum(-1)).max()
    return float(np.median(d[:, 1])), float(diam)


# =========================================================================
# Interactive questions. Returns a configuration dict (also saved as JSON at
# the end, so the same experiment can be re-run with --config).
# =========================================================================
def interactive_config(implementation="python"):
    # implementation = "python" (this script runs alpha_dualqp.py) or "cpp"
    # (make_inputs.py: the C++ program alphadual runs the computation).
    cpp = implementation == "cpp"
    print("\n=== Alpha filtration benchmark: dual active-set QP"
          + (" [C++ alphadual]" if cpp else "") + " ===\n")
    shapes = available_shapes()
    for i, s in enumerate(shapes, 1):
        print(f"  {i}. {s}")
    while True:
        s = ask("Shape (number or name)", "1")
        name = shapes[int(s) - 1] if s.isdigit() and 1 <= int(s) <= len(shapes) else s
        if name in shapes:
            break
        print("  ! unknown shape")
    fn = getattr(tadasets, name)
    doc = (inspect.getdoc(fn) or "").splitlines()
    print(f"\n{name}: {doc[0] if doc else ''}")
    print("Shape parameters (Enter keeps the tadasets default, 'None' for None):")
    params = {}
    for p, default in shape_params(fn).items():
        params[p] = parse_value(ask(f"  {p}", repr(default)))

    native = generate(fn, params, 10, 0, 0, 0).shape[1]
    # Native dimension = number of coordinates tadasets returns (e.g. 3 for
    # the torus); ambient dimensions must be at least this.
    print(f"\nNative dimension of '{name}' samples: {native}")

    n_list = ask_list("Number of points (e.g. 1000 or 500,1000 or 500:5000:500)",
                      int, "1000",
                      lambda v: require(all(x >= 3 for x in v), "need at least 3 points"))
    has_amb = "ambient" in inspect.signature(fn).parameters
    note = "" if has_amb else " (no 'ambient' arg in tadasets; tadasets.embed is used)"
    amb_list = ask_list(f"Ambient dimension(s){note}, >= {native}", int, str(native),
                        lambda v: require(all(a >= native for a in v),
                                          f"ambient dims must be >= {native}"))

    # Suggest r_max = 2.5 x median nearest-neighbour distance as a default.
    nn, diam = scale_hint(fn, params, n_list[0], native, 0)
    sugg = float(f"{2.5 * nn:.3g}")
    print(f"\nScale of a sample with n={n_list[0]}: median nearest-neighbour "
          f"distance {nn:.4g}, diameter ~{diam:.4g}.")
    print("Each run outputs the alpha filtration: every vertex, edge and triangle that\n"
          "enters at some radius between 0 and r_max, ordered by entry radius.\n"
          "r_max = inf gives the complete filtration (all radii, the whole Delaunay\n"
          "2-skeleton): exact, but every pair of points becomes a candidate, so the cost\n"
          "grows roughly like n^3 (e.g. ~40 s for n=1000 on a 3D torus); it needs DAQP.")
    r_list = ask_list("Maximum radius r_max (one or a list; 'inf' for all radii)", float,
                      str(sugg),
                      lambda v: require(all(x > 0 for x in v), "r_max must be positive"))
    print("\nMaximum simplex dimension d: 0 = vertices, 1 = edges, 2 = triangles,\n"
          "3 = tetrahedra, ... The cost grows quickly with d; in general position there\n"
          "are no simplices of dimension above the (affine) dimension of the data.")
    while True:
        try:
            max_dim = int(ask("Maximum simplex dimension d", "2"))
            require(max_dim >= 0, "d must be >= 0")
            break
        except ValueError as e:
            print(f"  ! {e}")
    runs = int(ask("Runs per configuration", "5"))
    base_seed = int(ask("Base seed (run i uses seed = base + i)", "0"))
    print("\nDual QP solver:\n"
          "  daqp   - DAQP (Arnstrom, Bemporad & Axehill 2022), the solver used by the\n"
          "           paper's authors (ref. [2]), via its "
          + ("C API" if cpp else "official Python bindings") + "\n"
          "  native - " + ("C++" if cpp else "Numba")
          + " re-implementation of the same dual active-set scheme")
    while True:
        # Re-ask until the choice is valid: the complete filtration (inf) needs
        # DAQP, and DAQP must be installed to be chosen.
        solver = ask("Solver", aq.DEFAULT_SOLVER).lower()
        if solver in aq.SOLVERS:
            if solver == "native" and any(np.isinf(r_list)):
                print("  ! r_max = inf (complete filtration) is only supported with daqp")
                continue
            if solver == "daqp" and not aq.HAVE_DAQP and not cpp:
                print("  ! DAQP is not installed (pip install daqp); choose 'native' or install it")
                continue
            break
        print(f"  ! choose one of {aq.SOLVERS}")
    if cpp:
        # The C++ program always writes the ordered filtration (binary file);
        # the witness map is not implemented there.
        witness = False
        save_filt = ask("Keep the alphadual output files (ordered filtrations) after "
                        "the Excel file is made? y/n", "y").lower().startswith("y")
    else:
        witness = ask("Also compute the witness map Phi (line 20 of Algorithm 1)? y/n",
                      "n").lower().startswith("y")
        save_filt = ask("Save the ordered filtration of every run to .npz files "
                        "(written after the timer stops)? y/n", "n").lower().startswith("y")
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    out = ask("Output Excel file",
              f"alpha_benchmark_{'cpp_' if cpp else ''}"
              f"{shape_label(dict(shape=name, shape_params=params))}_{stamp}.xlsx")
    if not out.lower().endswith(".xlsx"):
        out += ".xlsx"
    return dict(shape=name, shape_params=params, n_points=n_list,
                ambient_dims=amb_list, radii=r_list, runs=runs,
                base_seed=base_seed, solver=solver, compute_witness=witness,
                save_filtrations=save_filt, output=out, implementation=implementation,
                max_dim=max_dim)


# ---------------------------------------------------------------------------
# Experiment loop
# ---------------------------------------------------------------------------
# =========================================================================
# The experiment loop. For every (ambient dim, r_max, n) combination and
# every run: generate the cloud, time ONE call of the algorithm, then record
# statistics (outside the timer). Returns one dict per run ("rows").
# =========================================================================
# ---------------------------------------------------------------------------
# Saved point clouds. Every cloud generated with tadasets is written as a CSV
# file (one point per line, coordinates separated by commas, 17 significant
# digits: the doubles are reproduced EXACTLY when read back), so that other
# implementations can be run on identical inputs. index.csv describes them.
# Files are written after generation and before any timer starts.
# ---------------------------------------------------------------------------
def save_points_csv(path, X):
    np.savetxt(path, np.asarray(X, dtype=np.float64), fmt="%.17g", delimiter=",")


def shape_label(cfg):
    """Shape part of the file names: the tadasets function name, except that
    dsphere includes its dimension d (2sphere, 3sphere, ...)."""
    shape = cfg["shape"]
    if shape == "dsphere":
        d = cfg.get("shape_params", {}).get("d", 2)
        try:
            d = int(d) if float(d) == int(float(d)) else d
        except (TypeError, ValueError):
            pass
        return f"{d}sphere"
    return shape


def cloud_name(cfg, n, amb, run):
    """Name of a generated point cloud: <shape>_<points>_<ambient dimension>_<run>,
    e.g. torus_1000_3_1 or 2sphere_1000_5_1 (dsphere with d = 2). A cloud
    depends only on these (its seed is base_seed + run - 1), so runs with
    different r_max share the same cloud."""
    return f"{shape_label(cfg)}_{n}_{amb}_{run}"


def points_file_name(cfg, n, amb, run):
    return cloud_name(cfg, n, amb, run) + ".csv"


def run_name(cfg, n, amb, r, run):
    """Name for the RESULTS of one run: the cloud name, plus _r<r_max> when
    several r_max values are run (e.g. torus_1000_3_1_r0.4, ..._rinf)."""
    name = cloud_name(cfg, n, amb, run)
    if len(cfg["radii"]) > 1:
        name += "_r" + ("inf" if np.isinf(r) else f"{r:g}")
    return name


def write_points_index(folder, cfg, configs):
    """index.csv: one line per point cloud (file name, configuration, seed;
    r_max lists every r_max value run on that cloud)."""
    import csv
    params = ";".join(f"{k}={v}" for k, v in cfg["shape_params"].items())
    clouds = {}                                  # (n, amb, run) -> r_max values
    for amb, r, n in configs:
        for i in range(cfg["runs"]):
            clouds.setdefault((n, amb, i + 1), []).append("inf" if np.isinf(r) else f"{r:g}")
    with open(os.path.join(folder, "index.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        # shape = the label used in the file names (e.g. 2sphere); tadasets_function
        # = the generator (e.g. dsphere), with its parameters in shape_params
        w.writerow(["file", "shape", "tadasets_function", "shape_params", "n_points",
                    "ambient_dim", "run",
                    "seed", "r_max", "max_dim"])
        for (n, amb, run), radii in clouds.items():
            w.writerow([points_file_name(cfg, n, amb, run), shape_label(cfg), cfg["shape"],
                        params, n, amb, run,
                        cfg["base_seed"] + run - 1, ";".join(radii), cfg.get("max_dim", 2)])


# ---------------------------------------------------------------------------
# Crash safety. Every run is recorded in <output>_progress.jsonl: a "start"
# event before it, the full result row ("done") after it, both flushed to
# disk immediately. After every run the Excel file is rewritten (atomically).
# With --resume, finished runs are read back instead of being re-run; a run
# that was started but never finished killed the process (e.g. out of memory)
# and is recorded as an error instead of being retried (--retry-crashed
# retries it); a run interrupted with Ctrl-C is simply run again.
# ---------------------------------------------------------------------------
def _json_default(o):
    return o.item() if hasattr(o, "item") else str(o)   # numpy scalars


def append_progress(path, event):
    """Append one event to the progress file and force it to disk."""
    if path is None:
        return
    with open(path, "a") as fh:
        fh.write(json.dumps(event, default=_json_default) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def load_progress(path):
    """Read a progress file: (rows of finished runs by (config_id, run),
    runs started but never finished, time the experiment began)."""
    done, started, begin = {}, set(), None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                ev = json.loads(line)
            except json.JSONDecodeError:      # last line cut by a crash
                continue
            if ev["event"] == "begin":
                begin = begin or ev["time"]
            key = (ev.get("config_id"), ev.get("run"))
            if ev["event"] == "start":
                started.add(key)
            elif ev["event"] == "done":
                if ev["row"].get("crashed"):    # recorded crash: stays retryable
                    started.add(key)
                    done.pop(key, None)
                else:
                    done[key] = ev["row"]
                    started.discard(key)
            elif ev["event"] == "interrupted":
                started.discard(key)
    return done, started, begin


def run(cfg, progress=None, previous=None, retry_crashed=False, on_row=None):
    """progress: path of the progress file (None: no recording);
    previous: result of load_progress when resuming; on_row(rows, configs,
    t_begin, finished): called after every run (saves the Excel file)."""
    fn = getattr(tadasets, cfg["shape"])
    params = cfg["shape_params"]
    native = generate(fn, params, 10, 0, 0, 0).shape[1]
    cfg["native_dim"] = native
    # Validate the configuration (also needed for configs loaded from JSON,
    # which skip the interactive checks).
    cfg.setdefault("solver", aq.DEFAULT_SOLVER)       # configs saved before the option
    if cfg["solver"] not in aq.SOLVERS:
        raise ValueError(f"solver must be one of {aq.SOLVERS}")
    if cfg["solver"] == "daqp" and not aq.HAVE_DAQP:
        raise ImportError("solver 'daqp' needs the DAQP package: pip install daqp")
    # JSON may store radii as ints or as Infinity; make them floats.
    cfg["radii"] = [float(r) for r in cfg["radii"]]
    if cfg["solver"] == "native" and any(np.isinf(cfg["radii"])):
        raise ValueError("r_max = inf (complete filtration) needs solver 'daqp'")
    cfg.setdefault("save_filtrations", False)
    cfg.setdefault("max_dim", 2)                       # configs saved before the option
    # Folder for the optional .npz files with the ordered filtrations.
    filt_dir = os.path.splitext(cfg["output"])[0] + "_filtrations"
    points_dir = os.path.splitext(cfg["output"])[0] + "_points"      # saved clouds
    if cfg["save_filtrations"]:
        os.makedirs(filt_dir, exist_ok=True)
    bad = [a for a in cfg["ambient_dims"] if a < native]
    if bad:
        raise ValueError(f"ambient dims {bad} are below the native dimension {native}")

    # Compile everything once on a tiny example, so that the one-off Numba
    # compilation time is never counted in a measurement.
    print("\nWarming up (Numba JIT compilation, not timed)...", flush=True)
    t = time.perf_counter()
    aq.warmup(cfg["solver"])
    print(f"  done in {time.perf_counter() - t:.1f} s"
          + ("" if aq.HAVE_NUMBA else "  [Numba NOT installed: pure Python, very slow]"))

    # All combinations, ordered ambient -> r_max -> n, so that rows for the
    # same (ambient, r_max) are consecutive (used by the chart).
    configs = [(a, r, n) for a in cfg["ambient_dims"] for r in cfg["radii"]
               for n in cfg["n_points"]]
    os.makedirs(points_dir, exist_ok=True)
    write_points_index(points_dir, cfg, configs)
    rows = []
    done, crashed, begin = previous if previous else ({}, set(), None)
    t_begin = dt.datetime.fromisoformat(begin) if begin else dt.datetime.now()
    if not previous:
        append_progress(progress, {"event": "begin", "time": t_begin.isoformat()})
    total = len(configs) * cfg["runs"]
    for cid, (amb, r, n) in enumerate(configs, 1):
        for i in range(cfg["runs"]):
            key = (cid, i + 1)
            if key in done:                         # finished in an earlier session
                rows.append(done[key])
                print(f"[config {cid}/{len(configs)}] run {i + 1}: already done "
                      f"(read from the progress file)", flush=True)
                continue
            if key in crashed and not retry_crashed:
                row = dict(config_id=cid, n_points=n, ambient_dim=amb, radius=r,
                           run=i + 1, seed=cfg["base_seed"] + i,
                           status="error: the process was killed during this run "
                                  "(e.g. out of memory); not retried "
                                  "(use --retry-crashed to try again)",
                           crashed=True)
                rows.append(row)
                append_progress(progress, {"event": "done", "config_id": cid,
                                           "run": i + 1, "row": row})
                print(f"[config {cid}/{len(configs)}] run {i + 1}: crashed in an "
                      f"earlier session, recorded as an error", flush=True)
                if on_row:
                    on_row(rows, configs, t_begin, len(rows) == total)
                continue
            append_progress(progress, {"event": "start", "config_id": cid, "run": i + 1})
            # Run i always uses seed base+i: the same intrinsic cloud for every
            # ambient dimension, so the dimensions are compared on equal footing.
            seed = cfg["base_seed"] + i
            row = dict(config_id=cid, n_points=n, ambient_dim=amb, radius=r,
                       run=i + 1, seed=seed)
            try:
                # Point generation happens BEFORE the timer starts (not timed).
                X = generate(fn, params, n, amb, native, seed)
                pfile = os.path.join(points_dir, points_file_name(cfg, n, amb, i + 1))
                if not os.path.exists(pfile):      # same cloud for every r_max
                    save_points_csv(pfile, X)
                witness = cfg["compute_witness"]
                # ---- single timer: generated input -> all simplices (V, E, T)
                t0 = time.perf_counter()
                res = aq.alpha_complex_dualqp(X, r=r, compute_witness=witness,
                                              solver=cfg["solver"], max_dim=cfg["max_dim"])
                elapsed = time.perf_counter() - t0
                # ---- everything below is bookkeeping, outside the timer
                st = aq.compute_stats(res)
                row.update(status="ok", total_time_s=elapsed,
                           cech_edges=st["cech_edges"],
                           max_cech_degree=st["max_cech_degree"])
                fill_dimension_columns(row, [len(res["simplices"][k])
                                             for k in range(cfg["max_dim"] + 1)], st)
                # Largest entry radius = when the last simplex enters (values are
                # squared radii, sorted increasingly).
                vals = res["ordered"]["values"]
                row["max_entry_radius"] = (float(np.sqrt(max(vals[-1], 0.0)))
                                           if len(vals) else None)
                # Optionally save the ordered filtration (after the timer: not timed).
                # Named <shape>_<points>_<dim>_<run>[_r<r_max>].npz (see run_name).
                if cfg["save_filtrations"]:
                    fname = run_name(cfg, n, amb, r, i + 1) + ".npz"
                    o = res["ordered"]
                    np.savez_compressed(
                        os.path.join(filt_dir, fname),
                        simplices=o["simplices"], dims=o["dims"],
                        alpha_squared=o["values"],
                        radius=np.sqrt(np.maximum(o["values"], 0.0)),
                        points=X, seed=seed, r_max=r)
                    row["filtration_file"] = os.path.join(os.path.basename(filt_dir), fname)
                msg = (f"{elapsed:.4f} s  simplices per dim = {counts_string(row)}  "
                       f"max Cech degree {row['max_cech_degree']}")
                if row["failed"]:
                    msg += f"  [{row['failed']} QP failures]"
                del res, X
            # A failing run (e.g. out of memory) is recorded as an error row and the
            # experiment continues with the next run.
            except Exception as e:  # noqa: BLE001  (e.g. MemoryError)
                row["status"] = f"error: {type(e).__name__}: {e}"
                msg = row["status"]
            except KeyboardInterrupt:
                # Ctrl-C: save what is finished; this run is re-done on resume.
                append_progress(progress, {"event": "interrupted", "config_id": cid,
                                           "run": i + 1})
                if on_row:
                    on_row(rows, configs, t_begin, False)
                raise SystemExit(
                    "\nInterrupted. Finished runs are saved; continue with:\n"
                    f"  python run_experiments.py --config "
                    f"{os.path.splitext(cfg['output'])[0]}_config.json --resume")
            rows.append(row)
            append_progress(progress, {"event": "done", "config_id": cid,
                                       "run": i + 1, "row": row})
            print(f"[config {cid}/{len(configs)}] n={n} ambient={amb} r={r} "
                  f"run {i + 1}/{cfg['runs']} (seed {seed}): {msg}", flush=True)
            if on_row:                              # save the Excel file now
                on_row(rows, configs, t_begin, len(rows) == total)
    t_end = dt.datetime.now()
    return rows, configs, t_begin, t_end


# ---------------------------------------------------------------------------
# Excel output
# ---------------------------------------------------------------------------
# Line colours for the chart (one per series).
PALETTE = ["1F77B4", "D62728", "2CA02C", "FF7F0E", "9467BD", "8C564B",
           "E377C2", "7F7F7F", "BCBD22", "17BECF"]

# Column helpers. The simplex counts and QP statistics get one column per
# dimension 0..d (vertices, edges and triangles keep their usual names).
def count_key(k):
    """Row key of the number of k-simplices."""
    return {0: "n_vertices", 1: "n_edges", 2: "n_triangles"}.get(k, f"n_dim{k}")


def count_header(k):
    """Column header of the number of k-simplices."""
    return {0: "n_vertices", 1: "n_edges", 2: "n_triangles"}.get(k, f"n_simplices_dim{k}")


def qp_header(k):
    return {0: "QPs_vertices", 1: "QPs_edges", 2: "QPs_triangles"}.get(k, f"QPs_dim{k}")


def counts_string(row):
    """'5000/23830/29995' : the simplex counts of a row, dimension 0 upwards."""
    out, k = [], 0
    while count_key(k) in row and row[count_key(k)] is not None:
        out.append(str(row[count_key(k)]))
        k += 1
    return "/".join(out)


def fill_dimension_columns(row, counts, st):
    """Write simplex counts and QP statistics of dimensions 0..d into a row.
    counts[k] = number of k-simplices; st = statistics dict with 'dimK'
    entries (from alpha_dualqp.compute_stats or the C++ stats file).
    Dimensions 0-2 always get a column (None if not computed)."""
    d = len(counts) - 1
    for k in range(max(d, 2) + 1):
        row[count_key(k)] = counts[k] if k <= d else None
        row[f"qps_dim{k}"] = st[f"dim{k}"]["candidates"] if k <= d else None
    row["rejected"] = sum(st[f"dim{k}"]["rejected"] for k in range(d + 1))
    row["failed"] = sum(st[f"dim{k}"]["failed"] for k in range(d + 1))
    row["qp_iterations"] = sum(st[f"dim{k}"]["iterations"] for k in range(d + 1))


def run_cols(max_dim=2):
    """Columns of the 'Runs' sheet: (key in the row dict, column header, Excel
    number format), with one count and one QP column per dimension 0..max(d,2).
    Keys euler, qps_total and iters_per_qp are not stored in the rows: they
    are written as Excel formulas computed from other cells."""
    D = max(max_dim, 2)
    return ([
        ("config_id", "config_id", "0"), ("n_points", "n_points", "0"),
        ("ambient_dim", "ambient_dim", "0"), ("radius", "max_radius_r", "General"),
        ("run", "run", "0"), ("seed", "seed", "0"), ("status", "status", "@"),
        ("total_time_s", "alpha_complex_time_s", "0.0000"),
        # C++ only: part of alpha_complex_time_s spent reading the input file
        ("read_time_s", "input_read_time_s", "0.000000"),
        # C++ only: the computation alone (end of reading -> ordered filtration)
        ("compute_time_s", "compute_time_s", "0.0000")]
        + [(count_key(k), count_header(k), "#,##0") for k in range(D + 1)]
        + [("max_entry_radius", "largest_entry_radius", "0.0000"),
           ("euler", "euler_characteristic", "#,##0"),
           ("cech_edges", "cech_graph_edges", "#,##0"),
           ("max_cech_degree", "max_cech_degree", "#,##0")]
        + [(f"qps_dim{k}", qp_header(k), "#,##0") for k in range(D + 1)]
        + [("qps_total", "QPs_total", "#,##0"),
           ("rejected", "QPs_rejected", "#,##0"),
           ("failed", "QPs_failed", "#,##0"),
           ("qp_iterations", "active_set_iterations", "#,##0"),
           ("iters_per_qp", "iterations_per_QP", "0.00"),
           ("filtration_file", "filtration_file", "@")])


RUN_COLS = run_cols(2)    # default layout (d = 2)


# Excel column letter of the timing column in the 'Runs' sheet
TIME_COL = chr(ord("A") + [k for k, _, _ in RUN_COLS].index("total_time_s"))


# =========================================================================
# Write the Excel workbook: Inputs, Runs, Summary (live formulas), Chart.
# =========================================================================
def write_excel(cfg, rows, configs, t_begin, t_end, path, status=None):
    """status: shown in the Inputs sheet (e.g. "in progress: 3 of 10 runs");
    None means complete."""
    from openpyxl import Workbook
    from openpyxl.chart import Reference, ScatterChart, Series
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.utils import get_column_letter as col

    # Fonts and styles. Blue font (FI) marks the user's inputs.
    F = Font(name="Arial", size=10)
    FB = Font(name="Arial", size=10, bold=True)
    FH = Font(name="Arial", size=10, bold=True, color="FFFFFF")
    FT = Font(name="Arial", size=12, bold=True)
    FI = Font(name="Arial", size=10, color="0000FF")      # user inputs
    FN = Font(name="Arial", size=9, italic=True, color="555555")
    HFILL = PatternFill("solid", fgColor="1F4E78")
    thin = Border(bottom=Side(style="thin", color="999999"))

    # Write a styled header row (white bold text on dark blue).
    def header(ws, r, labels):
        for j, h in enumerate(labels, 1):
            c = ws.cell(row=r, column=j, value=h)
            c.font, c.fill, c.border = FH, HFILL, thin
            c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)

    wb = Workbook()
    # Ask Excel/LibreOffice to compute all formulas when the file is opened
    # (openpyxl writes formulas but cannot compute their values).
    wb.calculation.fullCalcOnLoad = True
    cpp = cfg.get("implementation", "python") == "cpp"   # results from alphadual?

    # ---------------- Inputs ----------------
    ws = wb.active
    ws.title = "Inputs"
    ws["A1"] = "Alpha complex benchmark (vertices, edges, triangles) - dual active-set QP"
    ws["A1"].font = FT
    ws["A2"] = ("Carlsson & Carlsson, 'Computing the alpha complex using dual active set "
                "quadratic programming', arXiv:2310.00536, Algorithm 1 with d = 2.")
    ws["A2"].font = FN
    header(ws, 4, ["Parameter", "Value", "Notes"])
    from importlib import metadata
    vers = {}
    # Installed versions of the libraries, for reproducibility.
    for mod in ("numpy", "scipy", "numba", "daqp", "tadasets", "openpyxl"):
        try:
            vers[mod] = metadata.version(mod)
        except Exception:  # noqa: BLE001
            vers[mod] = "not installed"
    fmt = lambda v: ", ".join(str(x) for x in v)  # noqa: E731
    # One row per item: (label, value, note, is_user_input).
    n_total = len(configs) * cfg["runs"]
    items = [
        ("Experiment status", status or f"complete ({len(rows)} of {n_total} runs)",
         "the file is saved after every run, so it is usable even if the experiment "
         "stops early", False),
        ("Shape (tadasets function)", cfg["shape"], "user input", True),
        *[(f"Shape parameter: {k}", repr(v) if v is None or isinstance(v, str) else v,
           "user input (tadasets argument)", True) for k, v in cfg["shape_params"].items()],
        ("Native dimension of samples", cfg["native_dim"], "columns returned by tadasets", False),
        ("Number of points n", fmt(cfg["n_points"]), "user input", True),
        ("Ambient dimension(s)", fmt(cfg["ambient_dims"]),
         "user input; ambient > native uses tadasets.embed (zero-pad + random rotation, "
         "an isometry: same complex, different cost)", True),
        ("Maximum radius r_max", fmt(cfg["radii"]),
         "user input; output = alpha filtration: all simplices entering at radius <= r_max "
         "(a1 = r_max^2, p = 0), ordered by entry radius (ties: lower dimension first); "
         "inf = complete filtration (whole Delaunay 2-skeleton)", True),
        (("alphadual output files kept", "yes" if cfg.get("save_filtrations") else "no",
          "user input; binary ordered filtration per run (read with "
          "alphadual_io.read_output), in <work dir>/outputs", True) if cpp else
         ("Ordered filtration saved", "yes" if cfg.get("save_filtrations") else "no",
          "user input; .npz per run (simplices, dims, alpha_squared, radius, points), "
          "written after the timer stops; folder <output>_filtrations", True)),
        ("Runs per configuration", cfg["runs"], "user input", True),
        ("Base seed", cfg["base_seed"],
         "user input; run i uses seed = base + i (sample and embedding rotation)", True),
        ("Witness map computed", "yes" if cfg["compute_witness"] else "no",
         "not implemented in the C++ program" if cpp else
         "user input (line 20 of Algorithm 1)", not cpp),
        ("Maximum simplex dimension d", cfg.get("max_dim", 2),
         "user input; simplices of dimensions 0..d (0 vertices, 1 edges, 2 triangles, "
         "3 tetrahedra, ...)", True),
        ("Configurations (n x ambient x r)", len(configs), "", False),
        (("Timed section (alpha_complex_time_s)",
          "inside the C++ program: from opening the input file to the ordered filtration "
          "up to dimension d "
          "(reading the input INCLUDED, writing the output excluded)",
          "one interval; a single extra clock reading after the input is read splits it "
          "into input_read_time_s (reading) and compute_time_s (the computation alone); "
          "no timers inside the computation; point generation and Excel writing are not "
          "timed", False) if cpp else
         ("Timed section (alpha_complex_time_s)",
          "one timer around one call: input point cloud -> all simplices up to dimension d",
          "no timers inside the computation; point generation/embedding, statistics and "
          "Excel writing are NOT timed; Numba JIT compiled in an untimed warm-up", False)),
        *([("Implementation", "C++ (alphadual)", "Python/bash only generate the inputs "
            "and collect the results", False),
           ("alphadual build", cfg.get("cpp_version", "unknown"),
            "from alphadual --version", False)] if cpp else
          [("Implementation", "Python (alpha_dualqp.py, Numba)", "", False)]),
        ("Standard deviation", "sample std (n-1)", "'n/a' when fewer than 2 successful runs", False),
        (("Dual QP solver", cfg["solver"],
          "user input; daqp = DAQP v0.10.3 (Arnstrom, Bemporad & Axehill, IEEE TAC 2022), "
          "ref. [2] of the paper, called through its C API (built without -DPROFILING, i.e. "
          "no timers inside DAQP); native = C++ port of our dual active-set scheme", True)
         if cpp else
         ("Dual QP solver", cfg["solver"],
          "user input; daqp = DAQP (Arnstrom, Bemporad & Axehill, IEEE TAC 2022), ref. [2] "
          "of the paper, official Python bindings; native = Numba re-implementation of the "
          "same dual active-set scheme. With daqp the time includes the Python->C call "
          "overhead of the DAQP bindings (one workspace per vertex, ~5-15 us per QP)", True)),
        *([("DAQP primal tolerance (rescaled problem)", aq.DAQP_PRIMAL_TOL,
            "each local QP is rescaled to unit size; early termination = DAQP fval_bound = c1",
            False),
           ("Gram factor tolerance", aq.GRAM_FACTOR_TOL,
            "B = G G^T by pivoted Cholesky (dot products only); stops at residual <= tol*max",
            False)] if cfg["solver"] == "daqp" else
          [("Dual QP pivot tolerance (relative)", aq.DEFAULT_PIV_TOL, "", False),
           ("Dual QP feasibility tolerance (relative)", aq.DEFAULT_FEAS_TOL, "", False),
           ("Max active-set iterations per QP", f"{aq.DEFAULT_ITER_FACTOR}*neighbours+50",
            "", False)]),
        ("QPs_rejected (Runs sheet)", "ruled out by the bound c1 or infeasible",
         "DAQP reports both with the same exit code, so they are counted together", False),
        *([] if cpp else
          [("Numba JIT used", "yes" if aq.HAVE_NUMBA else "no (pure Python)", "", False)]),
        ("Python", sys.version.split()[0], "", False),
        *[(f"{k} version", v, "", False) for k, v in vers.items()],
        ("Platform", platform.platform(), "", False),
        ("Processor", platform.processor() or platform.machine(), "", False),
        ("Logical CPUs", os.cpu_count(), "computation is single-threaded except BLAS", False),
        ("Started", t_begin.strftime("%Y-%m-%d %H:%M:%S"), "", False),
        ("Finished", t_end.strftime("%Y-%m-%d %H:%M:%S"), "", False),
        ("Wall-clock time of the whole script (s)", round((t_end - t_begin).total_seconds(), 2),
         "NOT a benchmark figure: includes point generation, Excel writing, etc.", False),
        ("Sum of alpha complex times over all runs (s)", f"=SUM(Runs!{TIME_COL}:{TIME_COL})",
         "sum of column alpha_complex_time_s in Runs (no point generation)", False),
    ]
    # Write the items starting at row 5 (rows 1-4: title, reference, header).
    for i, (k, v, note, is_input) in enumerate(items, 5):
        ws.cell(row=i, column=1, value=k).font = F
        c = ws.cell(row=i, column=2, value=v)
        c.font = FI if is_input else F
        c.alignment = Alignment(horizontal="left")
        ws.cell(row=i, column=3, value=note).font = FN
    ws.cell(row=len(items) + 6, column=1,
            value="Blue values are the user's inputs; everything else is recorded automatically."
            ).font = FN
    ws.column_dimensions["A"].width = 40
    ws.column_dimensions["B"].width = 42
    ws.column_dimensions["C"].width = 90

    # ---------------- Runs ----------------
    wr = wb.create_sheet("Runs")
    max_dim = cfg.get("max_dim", 2)
    RUN_COLS_ = run_cols(max_dim)          # columns for dimensions 0..max(d, 2)
    header(wr, 1, [h for _, h, _ in RUN_COLS_])
    # cidx: key -> column number; L: key -> column letter (for formulas).
    cidx = {k: j for j, (k, _, _) in enumerate(RUN_COLS_, 1)}
    L = {k: col(j) for k, j in cidx.items()}
    for i, row in enumerate(rows, 2):
        ok = row["status"] == "ok"
        for k, _, nf in RUN_COLS_:
            # Derived columns are Excel formulas (blank for failed runs).
            if k == "euler":
                # Euler characteristic: alternating sum of the counts, k = 0..d
                v = ("=" + "".join(f"{'+' if k % 2 == 0 else '-'}{L[count_key(k)]}{i}"
                                   for k in range(max_dim + 1)).lstrip("+")) if ok else None
            elif k == "qps_total":
                v = ("=" + "+".join(f"{L[f'qps_dim{k}']}{i}" for k in range(max_dim + 1))
                     if ok else None)
            elif k == "iters_per_qp":
                v = (f"=IF({L['qps_total']}{i}=0,\"\",{L['qp_iterations']}{i}/{L['qps_total']}{i})"
                     if ok else None)
            else:
                v = row.get(k)
                # Excel has no infinity: write the text "inf" instead.
                if isinstance(v, float) and np.isinf(v):
                    v = "inf"
            c = wr.cell(row=i, column=cidx[k], value=v)
            c.font, c.number_format = F, nf
    for j, (_, h, _) in enumerate(RUN_COLS_, 1):
        wr.column_dimensions[col(j)].width = max(10, min(22, len(h) + 2))
    wr.column_dimensions[L["status"]].width = 14
    wr.row_dimensions[1].height = 30
    wr.freeze_panes = "B2"
    wr.auto_filter.ref = f"A1:{col(len(RUN_COLS_))}{max(2, len(rows) + 1)}"

    # ---------------- Summary ----------------
    # One row per configuration. Every statistic is a formula over the Runs
    # sheet, restricted to that config_id and to successful runs, so the
    # sheet stays correct if you edit or delete rows in Runs.
    sm = wb.create_sheet("Summary")
    last = len(rows) + 1
    # rng(k) = absolute Excel range of column k in Runs, e.g. Runs!$H$2:$H$31.
    rng = lambda k: f"Runs!${L[k]}$2:${L[k]}${last}"  # noqa: E731
    timed = [("total_time_s", "alpha_complex_time")]
    if any(r.get("read_time_s") is not None for r in rows):      # C++ results
        timed.append(("read_time_s", "input_read_time"))
    if any(r.get("compute_time_s") is not None for r in rows):   # C++ results
        timed.append(("compute_time_s", "compute_time"))
    counts = [(count_key(k), "mean_" + count_header(k)) for k in range(max_dim + 1)]
    counts += [("max_cech_degree", "mean_max_cech_degree"),
              ("qps_total", "mean_QPs_total"), ("iters_per_qp", "mean_iterations_per_QP")]
    heads = ["config_id", "n_points", "ambient_dim", "max_radius_r", "successful_runs"]
    for _, h in timed:
        heads += [f"mean_{h}_s", f"std_{h}_s"]
    heads += [h for _, h in counts]
    header(sm, 1, heads)
    # Criteria shared by the AVERAGEIFS formulas: same config_id, status ok.
    crit = f"{rng('config_id')},$A{{r}},{rng('status')},\"ok\""
    for i, (amb, r, n) in enumerate(configs, 2):
        vals = [i - 1, n, amb, "inf" if np.isinf(r) else r,
                f"=COUNTIFS({rng('config_id')},$A{i},{rng('status')},\"ok\")"]
        j = 6
        for k, _ in timed:
            # Mean (AVERAGEIFS) and sample standard deviation. Excel has no
            # 'STDEV with conditions', so the std is written out explicitly:
            # sqrt( sum over matching runs of (time - mean)^2 / (count - 1) ),
            # with SUMPRODUCT doing the conditional sum. 'n/a' if fewer than 2 runs.
            mcell = f"{col(j)}{i}"
            vals.append(f"=IFERROR(AVERAGEIFS({rng(k)},{crit.format(r=i)}),\"\")")
            vals.append(
                f"=IF($E{i}>1,SQRT(SUMPRODUCT(({rng('config_id')}=$A{i})*"
                f"({rng('status')}=\"ok\")*({rng(k)}-{mcell})^2)/($E{i}-1)),\"n/a\")")
            j += 2
        for k, _ in counts:
            vals.append(f"=IFERROR(AVERAGEIFS({rng(k)},{crit.format(r=i)}),\"\")")
        for jj, v in enumerate(vals, 1):
            c = sm.cell(row=i, column=jj, value=v)
            c.font = F
            h = heads[jj - 1]
            c.number_format = ("0.000000" if "read_time" in h else
                               "0.0000" if h.endswith("_s") else
                               "0.00" if "iterations" in h else
                               "General" if h == "max_radius_r" else "#,##0")
    for j, h in enumerate(heads, 1):
        sm.column_dimensions[col(j)].width = max(11, min(20, len(h) + 1))
    sm.row_dimensions[1].height = 42
    sm.freeze_panes = "F2"
    sm.cell(row=len(configs) + 3, column=1,
            value="Means and sample standard deviations are live formulas over the 'Runs' "
                  "sheet (successful runs only, matched on config_id).").font = FN

    # ---------------- Chart ----------------
    # Chart: mean time vs n (one line per ambient dim and r_max) if several n
    # were given; otherwise vs ambient dimension if several were given;
    # otherwise no chart. 'groups' maps each line to its Summary rows.
    xs = None
    if len(cfg["n_points"]) > 1:
        xs, xcol, xlabel = "n", 2, "number of points n"
        groups = {}
        for i, (amb, r, n) in enumerate(configs, 2):
            groups.setdefault((amb, r), []).append(i)
        label = lambda g: f"ambient {g[0]}, r={g[1]}"  # noqa: E731
    elif len(cfg["ambient_dims"]) > 1:
        xs, xcol, xlabel = "ambient", 3, "ambient dimension"
        groups = {}
        for i, (amb, r, n) in enumerate(configs, 2):
            groups.setdefault((n, r), []).append(i)
        label = lambda g: f"n={g[0]}, r={g[1]}"  # noqa: E731
    if xs:
        ch = wb.create_sheet("Chart")
        ch["A1"] = "Mean alpha complex time (s), up to triangles, excluding point generation"
        ch["A1"].font = FT
        sc = ScatterChart()
        sc.title = "Mean alpha complex time"
        sc.scatterStyle = "lineMarker"
        sc.x_axis.title = xlabel
        sc.y_axis.title = "seconds"
        sc.height, sc.width = 11, 20
        for g, idx in groups.items():
            # The rows of one line must be consecutive in Summary (true by the
            # loop order in run()); each line plots column 6 (mean time).
            idx = sorted(idx, key=lambda k: configs[k - 2][0 if xs == "ambient" else 2])
            if idx != list(range(idx[0], idx[-1] + 1)):
                continue  # non-contiguous block (should not happen)
            xref = Reference(sm, min_col=xcol, min_row=idx[0], max_row=idx[-1])
            yref = Reference(sm, min_col=6, min_row=idx[0], max_row=idx[-1])
            s = Series(yref, xref, title=label(g))
            s.marker.symbol = "circle"
            s.smooth = False
            color = PALETTE[len(sc.series) % len(PALETTE)]
            s.graphicalProperties.line.solidFill = color
            s.marker.graphicalProperties.solidFill = color
            s.marker.graphicalProperties.line.solidFill = color
            sc.series.append(s)
        sc.x_axis.delete = False
        sc.y_axis.delete = False
        ch.add_chart(sc, "A3")

    # Atomic save: write a temporary file, then rename it over the old one,
    # so that a crash while saving never leaves a corrupted workbook.
    tmp = os.path.join(os.path.dirname(os.path.abspath(path)),
                       "." + os.path.basename(path) + ".tmp")
    wb.save(tmp)
    return replace_locked_aware(tmp, path)


def replace_locked_aware(tmp, path, attempts=10, wait_s=0.5):
    """Move tmp onto path. On Windows this fails while the file is open in
    another program (Excel, or briefly OneDrive / an antivirus scan): retry for
    a few seconds; if it stays locked, save as <name>_latest.xlsx instead and
    say so. Returns the path actually written."""
    for _ in range(attempts):
        try:
            os.replace(tmp, path)
            return path
        except PermissionError:
            time.sleep(wait_s)
    alt = os.path.splitext(path)[0] + "_latest.xlsx"
    try:
        os.replace(tmp, alt)
    except PermissionError:
        os.remove(tmp)
        raise PermissionError(
            f"{os.path.basename(path)} and {os.path.basename(alt)} are both locked "
            "(open in Excel?)") from None
    print(f"  note: {os.path.basename(path)} is locked (open in Excel?), so the latest "
          f"results were saved to {os.path.basename(alt)}. Close the file to let it "
          f"be updated again.", flush=True)
    return alt


# ---------------------------------------------------------------------------
# Short summary table on the terminal (same numbers as the Summary sheet).
def print_summary(rows, configs):
    print("\nSummary (successful runs):")
    print("(times = alpha complex computation only, point generation excluded)")
    print(f"{'n':>8} {'ambient':>8} {'r':>10} {'runs':>5} {'mean s':>10} {'std s':>10}"
          f"  mean simplices per dimension")
    for cid, (amb, r, n) in enumerate(configs, 1):
        rr = [x for x in rows if x["config_id"] == cid and x["status"] == "ok"]
        if not rr:
            print(f"{n:>8} {amb:>8} {r:>10.4g}  (no successful runs)")
            continue
        t = np.array([x["total_time_s"] for x in rr])
        sd = f"{t.std(ddof=1):.4f}" if len(t) > 1 else "n/a"
        dims = []                       # computed dimensions 0..d
        while count_key(len(dims)) in rr[0] and rr[0][count_key(len(dims))] is not None:
            dims.append(len(dims))
        means = "/".join(f"{np.mean([x[count_key(k)] for x in rr]):.0f}" for k in dims)
        print(f"{n:>8} {amb:>8} {r:>10.4g} {len(t):>5} {t.mean():>10.4f} {sd:>10}  {means}")


# Entry point: interactive questions (or a JSON config), run, write the
# workbook and the JSON config next to it, print a summary.
def main():
    try:                     # show every message at once (Git Bash buffers Python's output)
        sys.stdout.reconfigure(line_buffering=True)
    except AttributeError:
        pass
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", help="JSON configuration saved by a previous run")
    ap.add_argument("--output", help="override the output .xlsx path")
    ap.add_argument("--resume", action="store_true",
                    help="continue an interrupted experiment (with --config "
                         "<output>_config.json): finished runs are not repeated")
    ap.add_argument("--retry-crashed", action="store_true",
                    help="with --resume: also re-run runs that killed the process")
    args = ap.parse_args()
    if args.resume and not args.config:
        ap.error("--resume needs --config (the *_config.json saved at the start)")

    if args.config:
        with open(args.config) as fh:
            cfg = json.load(fh)
    else:
        cfg = interactive_config()
    if args.output:
        cfg["output"] = args.output

    out = cfg["output"]
    base = os.path.splitext(out)[0]
    cfg_path = base + "_config.json"
    progress = base + "_progress.jsonl"
    previous = None
    if args.resume:
        if not os.path.exists(progress):
            sys.exit(f"nothing to resume: {progress} not found")
        previous = load_progress(progress)
        print(f"Resuming: {len(previous[0])} run(s) already done"
              + (f", {len(previous[1])} crashed" if previous[1] else ""))
    else:
        open(progress, "w").close()                 # a new experiment
        # Save the configuration at the START, so --resume works after a crash.
        with open(cfg_path, "w") as fh:
            json.dump({k: v for k, v in cfg.items() if k != "native_dim"}, fh, indent=2)

    def save(rows, configs, t_begin, finished):
        n_total = len(configs) * cfg["runs"]
        try:
            write_excel(cfg, rows, configs, t_begin, dt.datetime.now(), out,
                        status=None if finished else
                        f"in progress: {len(rows)} of {n_total} runs done")
        except Exception as e:  # noqa: BLE001  (e.g. the file is open in Excel)
            print(f"  (warning: the Excel report was not updated this time: {e}; "
                  f"no results are lost, they are in {progress})", flush=True)

    rows, configs, t0, t1 = run(cfg, progress=progress, previous=previous,
                                retry_crashed=args.retry_crashed, on_row=save)
    written = write_excel(cfg, rows, configs, t0, t1, out)    # final version
    print_summary(rows, configs)
    print(f"\nExcel report : {os.path.abspath(written)}  (also saved after every run)")
    print(f"Config (JSON): {os.path.abspath(cfg_path)}  (re-run with --config)")
    print(f"Progress     : {os.path.abspath(progress)}  (resume with --config ... --resume)")


if __name__ == "__main__":
    main()
