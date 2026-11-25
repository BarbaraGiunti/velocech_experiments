"""
collect_results.py
==================

Stage 3 of the C++ benchmark pipeline (see run_benchmark.sh): read what
`alphadual` wrote for every job (statistics JSON + ordered filtration) and
build the same Excel report as run_experiments.py (sheets Inputs, Runs,
Summary, Chart), plus the configuration JSON next to it.

    python collect_results.py WORKDIR

The times come from INSIDE the C++ program (alpha_time_s: reading the input
+ computing the ordered filtration; read_time_s: the reading part alone).
Nothing measured here in Python enters the benchmark figures.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys

import numpy as np

import alphadual_io as io
import run_experiments as rx   # Excel writer and terminal summary are shared


def collect(work, partial=False):
    """partial=True: called while jobs are still running; jobs not run yet are
    left out of the report (their status files do not exist yet)."""
    with open(os.path.join(work, "manifest.json")) as fh:
        man = json.load(fh)
    cfg = man["config"]
    cfg.setdefault("max_dim", 2)
    configs = [tuple(c) for c in man["configs"]]
    rows = []
    for j in man["jobs"]:
        row = {k: j[k] for k in ("config_id", "n_points", "ambient_dim", "radius",
                                 "run", "seed")}
        out = os.path.join(work, j["output"])
        err = out + ".error"
        if partial and not (os.path.exists(out + ".json") or os.path.exists(err)):
            continue                                     # not run yet
        try:
            if os.path.exists(err):                      # written by the bash loop
                raise RuntimeError(open(err).read().strip() or "alphadual failed")
            st = io.read_stats(out + ".json")
            row.update(
                status="ok",
                total_time_s=st["alpha_time_s"],         # read + compute (C++ clock)
                read_time_s=st["read_time_s"],
                compute_time_s=st["compute_time_s"],     # computation alone
                cech_edges=st["cech_edges"], max_cech_degree=st["max_cech_degree"])
            # simplex counts and QP statistics for dimensions 0..d
            rx.fill_dimension_columns(row, st["n_simplices"], st)
            # Largest entry radius = last value of the ordered filtration (in the
            # stats file since alphadual 1.1; older versions: read the output).
            if "max_value" in st:
                mv = st["max_value"]
            else:
                vals = io.read_output(out)["values"]
                mv = float(vals[-1]) if len(vals) else None
            row["max_entry_radius"] = float(np.sqrt(max(mv, 0.0))) if mv is not None else None
            if cfg.get("save_filtrations", True):
                row["filtration_file"] = j["output"]
            elif os.path.exists(out):
                os.remove(out)                           # user chose not to keep it
        except Exception as e:  # noqa: BLE001
            row["status"] = "error: " + " | ".join(str(e).split("\n"))
        rows.append(row)
    return cfg, configs, rows, man


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("workdir")
    ap.add_argument("--partial", action="store_true",
                    help="jobs are still running: report the finished ones, quietly")
    args = ap.parse_args()
    work = os.path.abspath(args.workdir)

    cfg, configs, rows, man = collect(work, partial=args.partial)
    t_begin = dt.datetime.fromisoformat(man["started"])
    t_end = dt.datetime.now()
    xlsx = os.path.join(work, os.path.basename(cfg["output"]))
    n_total = len(man["jobs"])
    status = (f"in progress: {len(rows)} of {n_total} jobs done"
              if len(rows) < n_total else None)
    try:
        written = rx.write_excel(cfg, rows, configs, t_begin, t_end, xlsx, status=status)
    except Exception as e:  # noqa: BLE001
        # Say why (the bash loop then prints its warning and continues; every
        # job's results stay on disk and are included in the next update).
        print(f"  could not update the Excel report: {type(e).__name__}: {e}", flush=True)
        sys.exit(1)
    if args.partial:
        sys.exit(0)
    xlsx = written or xlsx
    cfg_path = os.path.splitext(xlsx)[0] + "_config.json"
    with open(cfg_path, "w") as fh:
        json.dump({k: v for k, v in cfg.items() if k not in ("native_dim", "cpp_version")},
                  fh, indent=2)
    rx.print_summary(rows, configs)
    ok = [r for r in rows if r["status"] == "ok"]
    if ok:
        share = max(r["read_time_s"] / r["total_time_s"] for r in ok if r["total_time_s"] > 0)
        print(f"\nInput reading: at most {100 * share:.2f}% of alpha_complex_time "
              f"over all runs (column input_read_time_s)")
    failed = len(rows) - len(ok)
    if failed:
        print(f"\n{failed} run(s) failed: see the status column of the Runs sheet")
    print(f"\nExcel report : {xlsx}")
    print(f"Config (JSON): {cfg_path}  (re-run with: run_benchmark.sh DIR --config ...)")
    sys.exit(0)


if __name__ == "__main__":
    main()
