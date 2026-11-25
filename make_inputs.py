"""
make_inputs.py
==============

Stage 1 of the C++ benchmark pipeline (see run_benchmark.sh):
ask the same questions as run_experiments.py (shape, parameters, numbers of
points, ambient dimensions, r_max, maximum simplex dimension d, runs, seed,
solver), generate every point
cloud with tadasets, and write them as binary input files for `alphadual`.

    python make_inputs.py WORKDIR                     # interactive
    python make_inputs.py WORKDIR --config cfg.json   # re-use a configuration

Creates in WORKDIR:
    inputs/<shape>_<points>_<dim>_<run>.bin   one point cloud per run (alphadual
                              input format), e.g. inputs/torus_1000_3_1.bin
    points/<shape>_<points>_<dim>_<run>.csv   the same clouds as CSV (exact values, one
                              point per line) + points/index.csv describing them
    outputs/                  where alphadual will write its results
    jobs.tsv                  one line per run: job id, input, output, r, solver, max_dim
    manifest.json             configuration + description of every job
    config.json               the configuration alone (for --config)

Point clouds are generated exactly as in run_experiments.py (same seeds, same
embedding), so the Python and C++ benchmarks see identical inputs.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import subprocess
import sys

import numpy as np

import alphadual_io as io
import run_experiments as rx   # prompts and point generation are shared

HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------------------
# File names and index.csv -- defined HERE (not taken from run_experiments.py),
# so that the names do not depend on the version of the other files.
# ---------------------------------------------------------------------------
def shape_label(cfg):
    """Shape part of the names: the tadasets function name, except that
    dsphere includes its dimension d: 2sphere, 3sphere, ..."""
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
    """<shape>_<points>_<ambient dimension>_<run>, e.g. 2sphere_1000_5_1. A cloud
    depends only on these (seed = base_seed + run - 1): runs with different
    r_max share it."""
    return f"{shape_label(cfg)}_{n}_{amb}_{run}"


def run_name(cfg, n, amb, r, run):
    """Results of one run: the cloud name, plus _r<r_max> if several r_max."""
    name = cloud_name(cfg, n, amb, run)
    if len(cfg["radii"]) > 1:
        name += "_r" + ("inf" if np.isinf(r) else f"{r:g}")
    return name


def save_points_csv(path, X):
    """One point per line, coordinates separated by commas, 17 significant
    digits (the numbers read back are exactly the same doubles)."""
    np.savetxt(path, np.asarray(X, dtype=np.float64), fmt="%.17g", delimiter=",")


def write_points_index(folder, cfg, configs):
    """index.csv: one line per point cloud (name, shape, tadasets function and
    its parameters, points, ambient dimension, run, seed, r_max values, d)."""
    import csv
    params = ";".join(f"{k}={v}" for k, v in cfg["shape_params"].items())
    clouds = {}
    for amb, r, n in configs:
        for i in range(cfg["runs"]):
            clouds.setdefault((n, amb, i + 1), []).append("inf" if np.isinf(r) else f"{r:g}")
    with open(os.path.join(folder, "index.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["file", "shape", "tadasets_function", "shape_params", "n_points",
                    "ambient_dim", "run", "seed", "r_max", "max_dim"])
        for (n, amb, run), radii in clouds.items():
            w.writerow([cloud_name(cfg, n, amb, run) + ".csv", shape_label(cfg), cfg["shape"],
                        params, n, amb, run, cfg["base_seed"] + run - 1, ";".join(radii),
                        cfg.get("max_dim", 2)])


def main():
    try:                     # show every message at once (Git Bash buffers Python's output)
        sys.stdout.reconfigure(line_buffering=True)
    except AttributeError:
        pass
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("workdir", help="folder for inputs, outputs and results")
    ap.add_argument("--config", help="JSON configuration (skips the questions)")
    args = ap.parse_args()

    started = dt.datetime.now()
    if args.config:
        with open(args.config) as fh:
            cfg = json.load(fh)
    else:
        cfg = rx.interactive_config(implementation="cpp")
    cfg["implementation"] = "cpp"
    cfg["compute_witness"] = False
    cfg.setdefault("save_filtrations", True)
    cfg.setdefault("solver", "daqp")
    cfg.setdefault("max_dim", 2)                 # configs saved before the option
    if int(cfg["max_dim"]) < 0:
        sys.exit("max_dim must be >= 0")
    cfg["radii"] = [float(r) for r in cfg["radii"]]
    if cfg["solver"] not in ("daqp", "native"):
        sys.exit("solver must be daqp or native")
    if cfg["solver"] == "native" and any(np.isinf(cfg["radii"])):
        sys.exit("r_max = inf (complete filtration) needs solver daqp")

    # Record which build of alphadual will run (written into the Excel file).
    exe = io.DEFAULT_EXE                      # alphadual (alphadual.exe on Windows)
    try:
        cfg["cpp_version"] = subprocess.run([exe, "--version"], capture_output=True,
                                            text=True, check=True).stdout.strip()
    except Exception:  # noqa: BLE001
        sys.exit(f"cannot run {exe} --version: build it first (make)")

    import tadasets
    fn = getattr(tadasets, cfg["shape"])
    params = cfg["shape_params"]
    native = rx.generate(fn, params, 10, 0, 0, 0).shape[1]
    cfg["native_dim"] = native
    bad = [a for a in cfg["ambient_dims"] if a < native]
    if bad:
        sys.exit(f"ambient dims {bad} are below the native dimension {native}")

    work = os.path.abspath(args.workdir)
    os.makedirs(os.path.join(work, "inputs"), exist_ok=True)
    os.makedirs(os.path.join(work, "outputs"), exist_ok=True)
    os.makedirs(os.path.join(work, "points"), exist_ok=True)

    # Same order of configurations as run_experiments.run(): ambient -> r -> n.
    configs = [(a, r, n) for a in cfg["ambient_dims"] for r in cfg["radii"]
               for n in cfg["n_points"]]
    jobs = []
    print(f"\nWriting the input files to {os.path.join(work, 'inputs')} and "
          f"{os.path.join(work, 'points')} ...")
    for cid, (amb, r, n) in enumerate(configs, 1):
        for i in range(cfg["runs"]):
            seed = cfg["base_seed"] + i
            X = rx.generate(fn, params, n, amb, native, seed)
            # Files named after the cloud: <shape>_<points>_<dim>_<run>; results
            # also get _r<r_max> when several r_max values are run.
            stem = cloud_name(cfg, n, amb, i + 1)
            job = run_name(cfg, n, amb, r, i + 1)
            inp = os.path.join("inputs", stem + ".bin")
            out = os.path.join("outputs", job + ".alpha")
            if not os.path.exists(os.path.join(work, inp)):    # same cloud for every r_max
                io.write_input(os.path.join(work, inp), X)
                # the same cloud as CSV (exact values), for other implementations
                save_points_csv(os.path.join(work, "points", stem + ".csv"), X)
            jobs.append(dict(job=job, config_id=cid, n_points=n, ambient_dim=amb,
                             radius=r, run=i + 1, seed=seed, input=inp, output=out))

    write_points_index(os.path.join(work, "points"), cfg, configs)
    # jobs.tsv is what the bash loop reads: paths are relative to WORKDIR.
    # newline="\n": Unix line endings even on Windows (read by bash)
    with open(os.path.join(work, "jobs.tsv"), "w", newline="\n") as fh:
        for j in jobs:
            r = "inf" if np.isinf(j["radius"]) else repr(float(j["radius"]))
            fh.write(f"{j['job']}\t{j['input']}\t{j['output']}\t{r}\t{cfg['solver']}"
                     f"\t{int(cfg['max_dim'])}\n")
    cfg_clean = {k: v for k, v in cfg.items() if k not in ("native_dim", "cpp_version")}
    with open(os.path.join(work, "config.json"), "w") as fh:
        json.dump(cfg_clean, fh, indent=2)
    with open(os.path.join(work, "manifest.json"), "w") as fh:
        json.dump(dict(config=cfg, configs=configs, jobs=jobs,
                       started=started.isoformat(timespec="seconds")), fh, indent=2)
    print(f"Done: point clouds in {os.path.join(work, 'points')} "
          f"({len(jobs)} jobs listed in {os.path.join(work, 'jobs.tsv')})")


if __name__ == "__main__":
    main()
