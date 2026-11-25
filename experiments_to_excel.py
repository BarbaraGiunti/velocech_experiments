"""
experiments_to_excel.py
=======================

Builds the Excel summary of ONE input subfolder from the records written by
run_all_experiments.sh (called by it after every run; can be re-run by hand):

    python experiments_to_excel.py RESULTS/<subfolder>

Reads in that folder:
    runs.jsonl     one line per run (measure_run.py): algorithm, input,
                   status, wall time, peak memory measured from outside
    reported.tsv   values reported by each code itself (run id, key=value),
                   e.g. time_s, read_time_s, n_dim0, ...
    inputs.json    number of points and ambient dimension of each input
and ../settings.json (the experiment's settings). Writes summary.xlsx:
    Settings          settings of the experiment and meaning of the columns
    Runs              one row per (algorithm, input)
    Summary           one row per algorithm (live formulas over Runs)
    Time by input     inputs x algorithms: time_s (the code's own time)
    Memory by input   inputs x algorithms: peak memory (from outside)
"""

from __future__ import annotations

import json
import os
import sys

PRIORITY = ["time_s", "read_time_s", "compute_time_s", "alpha_time_s", "diameter_time_s",
            "diameter", "r"]


def load(folder):
    runs = []
    p = os.path.join(folder, "runs.jsonl")
    if os.path.exists(p):
        for line in open(p):
            line = line.strip()
            if line:
                try:
                    runs.append(json.loads(line))
                except json.JSONDecodeError:          # a line cut by a crash
                    pass
    reported = {}
    p = os.path.join(folder, "reported.tsv")
    if os.path.exists(p):
        for line in open(p):
            parts = line.rstrip("\n").split("\t", 1)
            if len(parts) == 2 and "=" in parts[1]:
                k, _, v = parts[1].partition("=")
                reported.setdefault(parts[0], {})[k.strip()] = v.strip()
    inputs = {}
    p = os.path.join(folder, "inputs.json")
    if os.path.exists(p):
        inputs = json.load(open(p))
    settings = {}
    p = os.path.join(os.path.dirname(os.path.abspath(folder)), "settings.json")
    if os.path.exists(p):
        settings = json.load(open(p))
    return runs, reported, inputs, settings


# Windows exit codes of programs that crashed or could not start
WINDOWS_CODES = {
    3221225781: "the program could not start: a DLL it needs was not found "
                "(e.g. the MinGW runtime libraries are not on the PATH)",
    3221225785: "the program could not start: a DLL of the wrong version was found "
                "(e.g. Anaconda's libstdc++/libgomp before MinGW's on the PATH)",
    3221225595: "the program could not start: it is not a valid program for this Windows "
                "(32/64-bit mismatch of the program or a DLL)",
    3221225477: "the program crashed (access violation)",
    3221226505: "the program crashed (stack buffer overrun / abort)",
    3221225725: "the program crashed (stack overflow)",
    3221225794: "the program could not start (initialization failed)",
}


def exit_code_message(code):
    try:
        c = int(code)
    except (TypeError, ValueError):
        return ""
    return WINDOWS_CODES.get(c % (1 << 32), "")


def stderr_tail(folder, algo, name, n=2):
    """Last non-empty lines of a failed run's saved error output."""
    p = os.path.join(folder, "outputs", str(algo), f"{name}.stderr")
    try:
        lines = [l.strip() for l in open(p, errors="replace")
                 if l.strip() and l.strip().strip("^~ ")]      # skip traceback markers
        return " | ".join(lines[-n:])[:300]
    except OSError:
        return ""


def num(v):
    try:
        return float(v) if v not in (None, "") else None
    except (TypeError, ValueError):
        return v


def build(folder):
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter as col

    runs, reported, inputs, settings = load(folder)
    sub = os.path.basename(os.path.abspath(folder))
    # Rows: the last record of each (algorithm, input) -- a resumed run replaces
    # an earlier interrupted one.
    by_id = {}
    for r in runs:
        by_id[r.get("id", f"{r.get('algorithm')}|{r.get('input')}")] = r
    rows = []
    extra_keys = set()
    for rid, r in by_id.items():
        rep = reported.get(rid, {})
        row = dict(subfolder=sub, input=r.get("input"), algorithm=r.get("algorithm"),
                   status=r.get("status"), exit_code=r.get("exit_code"),
                   wall_time_s=r.get("wall_time_s"), peak_memory_mb=r.get("peak_memory_mb"),
                   message=r.get("message", ""))
        info = inputs.get(r.get("input"), {})
        row["n_points"], row["ambient_dim"] = info.get("n_points"), info.get("ambient_dim")
        row["input_group"] = info.get("group")
        # a code may report its own failure (e.g. an error caught inside it)
        if row["status"] == "ok" and rep.get("status", "ok") != "ok":
            row["status"] = "error"
            row["message"] = rep.get("status")
        if row["status"] != "ok" and not row["message"]:
            row["message"] = (rep.get("message", "")
                              or stderr_tail(folder, r.get("algorithm"), r.get("input"))
                              or exit_code_message(r.get("exit_code"))
                              or (f"exit code {r.get('exit_code')}" if r.get("exit_code") is not None
                                  else ""))
        for k, v in rep.items():
            if k in ("status", "message"):
                continue
            row[k] = num(v)
            extra_keys.add(k)
        rows.append(row)
    algos = list(dict.fromkeys(r["algorithm"] for r in rows))
    names = list(dict.fromkeys(r["input"] for r in rows))
    rows.sort(key=lambda r: (algos.index(r["algorithm"]), names.index(r["input"])))
    extra = [k for k in PRIORITY if k in extra_keys and k != "time_s"]
    dims = sorted((k for k in extra_keys if k.startswith("n_dim") and k[5:].isdigit()),
                  key=lambda k: int(k[5:]))
    extra += sorted(k for k in extra_keys if k not in PRIORITY and k not in dims) + dims
    cols = (["subfolder", "input", "input_group", "n_points", "ambient_dim", "algorithm", "status",
             "exit_code", "wall_time_s", "peak_memory_mb", "time_s"] + extra + ["message"])

    F = Font(name="Arial", size=10)
    FH = Font(name="Arial", size=10, bold=True, color="FFFFFF")
    FN = Font(name="Arial", size=9, italic=True, color="555555")
    FT = Font(name="Arial", size=12, bold=True)
    FILL = PatternFill("solid", fgColor="1F4E78")

    def header(ws, labels, row=1):
        for j, h in enumerate(labels, 1):
            c = ws.cell(row=row, column=j, value=h)
            c.font, c.fill = FH, FILL
            c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        ws.row_dimensions[row].height = 30

    wb = Workbook()
    wb.calculation.fullCalcOnLoad = True

    # ---- Settings
    ws = wb.active
    ws.title = "Settings"
    ws["A1"] = f"Experiments on input subfolder '{sub}'"
    ws["A1"].font = FT
    header(ws, ["Setting", "Value"], row=3)
    items = list(settings.items()) + [
        ("wall_time_s", "measured from outside (measure_run.py): start of the program to its "
                        "end, including process start, reading and writing"),
        ("peak_memory_mb", "measured from outside: " + next(
            (r.get("memory_method") for r in runs if r.get("memory_method")), "")
         + ". On Linux/macOS the command starts as a copy of the measuring Python "
           "process, so values below ~15 MB mean 'small' rather than an exact figure; "
           "on Windows only the command itself is counted"),
        ("time_s", "the code's OWN measurement of its computation, reading excluded: "
                   "alphadual = compute_time_s; half_diameter = diameter + alpha; "
                   "gudhi = compute_time_s; velocech = distance_time_s + h1_time_s "
                   "(its 'Time for distances' + 'Time after computing distances')"),
        ("velocech", "computes the full Cech H1 barcode (r_max and the maximum dimension do "
                     "not apply); n_h1_bars = number of H1 bars; bars on a diameter scale"),
        ("status", "ok / timeout (stopped at the time limit) / error / skipped (a smaller "
                   "input of the same group failed, see --skip-larger)"),
        ("input_group", "shape + shape parameters (noise, ...) + ambient dimension: inputs of a "
                        "group differ only in size and run"),
    ]
    for i, (k, v) in enumerate(items, 4):
        ws.cell(row=i, column=1, value=k).font = F
        c = ws.cell(row=i, column=2, value=v if not isinstance(v, (list, dict)) else json.dumps(v))
        c.font = F
        c.alignment = Alignment(wrap_text=True, vertical="top")
    ws.column_dimensions["A"].width = 26
    ws.column_dimensions["B"].width = 110

    # ---- Runs
    wr = wb.create_sheet("Runs")
    header(wr, cols)
    for i, r in enumerate(rows, 2):
        for j, h in enumerate(cols, 1):
            v = r.get(h)
            c = wr.cell(row=i, column=j, value=v)
            c.font = F
            if h.endswith("_s"):
                c.number_format = "0.000000"
            elif h == "peak_memory_mb":
                c.number_format = "0.0"
            elif h.startswith("n_"):
                c.number_format = "#,##0"
    for j, h in enumerate(cols, 1):
        wr.column_dimensions[col(j)].width = 40 if h == "message" else max(10, min(24, len(h) + 3))
    wr.freeze_panes = "C2"
    if rows:
        wr.auto_filter.ref = f"A1:{col(len(cols))}{len(rows) + 1}"
    last = max(2, len(rows) + 1)
    C = {h: col(j) for j, h in enumerate(cols, 1)}
    rng = lambda h: f"Runs!${C[h]}$2:${C[h]}${last}"  # noqa: E731

    # ---- Summary (per algorithm)
    sm = wb.create_sheet("Summary")
    heads = ["algorithm", "runs", "ok", "timeout", "error", "skipped", "mean_time_s (ok)",
             "std_time_s (ok)", "mean_wall_time_s (ok)", "mean_peak_memory_mb (ok)",
             "max_peak_memory_mb"]
    header(sm, heads)
    for i, a in enumerate(algos, 2):
        A, ok = f"$A{i}", f"{rng('status')},\"ok\""
        vals = [a,
                f"=COUNTIF({rng('algorithm')},{A})",
                f"=COUNTIFS({rng('algorithm')},{A},{ok})",
                f"=COUNTIFS({rng('algorithm')},{A},{rng('status')},\"timeout\")",
                f"=COUNTIFS({rng('algorithm')},{A},{rng('status')},\"error\")",
                f"=COUNTIFS({rng('algorithm')},{A},{rng('status')},\"skipped\")",
                f"=IFERROR(AVERAGEIFS({rng('time_s')},{rng('algorithm')},{A},{ok}),\"\")",
                f"=IF(C{i}>1,SQRT(SUMPRODUCT(({rng('algorithm')}={A})*({rng('status')}=\"ok\")*"
                f"({rng('time_s')}-G{i})^2)/(C{i}-1)),\"n/a\")",
                f"=IFERROR(AVERAGEIFS({rng('wall_time_s')},{rng('algorithm')},{A},{ok}),\"\")",
                f"=IFERROR(AVERAGEIFS({rng('peak_memory_mb')},{rng('algorithm')},{A},{ok}),\"\")",
                # max over this algorithm's runs (works in every Excel version)
                f"=SUMPRODUCT(MAX(({rng('algorithm')}={A})*{rng('peak_memory_mb')}))"]
        for j, v in enumerate(vals, 1):
            c = sm.cell(row=i, column=j, value=v)
            c.font = F
            c.number_format = "0.000000" if 7 <= j <= 9 else ("0.0" if j >= 10 else "0")
    for j, h in enumerate(heads, 1):
        sm.column_dimensions[col(j)].width = max(12, len(h) + 2)
    sm.cell(row=len(algos) + 3, column=1,
            value="time_s = the code's own time (computation, reading excluded); wall time and "
                  "memory measured from outside. Formulas over the Runs sheet.").font = FN

    # ---- Inputs x algorithms
    status_of = {(r["input"], r["algorithm"]): r["status"] for r in rows}
    for title, key, fmt, only_ok in (("Time by input", "time_s", "0.000000", True),
                                     ("Memory by input", "peak_memory_mb", "0.0", False)):
        wt = wb.create_sheet(title)
        header(wt, ["input", "n_points", "ambient_dim"] + algos)
        for i, name in enumerate(names, 2):
            info = inputs.get(name, {})
            for j, v in enumerate([name, info.get("n_points"), info.get("ambient_dim")], 1):
                wt.cell(row=i, column=j, value=v).font = F
            for j, a in enumerate(algos, 4):
                st = status_of.get((name, a))
                if st is None:
                    v = ""
                elif only_ok and st != "ok":
                    v = st                      # timeout / error: shown as text
                else:
                    v = (f"=IFERROR(AVERAGEIFS({rng(key)},{rng('input')},$A{i},"
                         f"{rng('algorithm')},{col(j)}$1),\"\")")
                c = wt.cell(row=i, column=j, value=v)
                c.font = F
                c.number_format = fmt
        wt.column_dimensions["A"].width = 32
        for j in range(2, len(algos) + 4):
            wt.column_dimensions[col(j)].width = 16
        wt.freeze_panes = "B2"
        note = ("the code's own time (s), reading excluded" if key == "time_s"
                else "peak memory (MB) measured from outside, also for runs stopped at the limit")
        wt.cell(row=len(names) + 3, column=1, value=note).font = FN

    # ---- atomic save (tolerates the file being open in Excel)
    path = os.path.join(folder, "summary.xlsx")
    tmp = os.path.join(folder, ".summary.xlsx.tmp")
    wb.save(tmp)
    import time
    for _ in range(10):
        try:
            os.replace(tmp, path)
            return path
        except PermissionError:
            time.sleep(0.5)
    alt = os.path.join(folder, "summary_latest.xlsx")
    os.replace(tmp, alt)
    print(f"  note: summary.xlsx is open in another program; saved {alt} instead")
    return alt


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python experiments_to_excel.py RESULTS/<subfolder>")
    try:
        print(build(sys.argv[1]))
    except Exception as e:  # noqa: BLE001
        sys.exit(f"could not write the Excel summary: {type(e).__name__}: {e}")
