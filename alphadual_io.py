"""
alphadual_io.py
===============

Python side of the C++ `alphadual` program: write its input files, run it,
and read back its output (the ordered alpha filtration) and statistics.

File formats (little-endian):

  input  : b"ALPHAIN1" | int64 N | int64 m | int64 has_weights |
           float64 X[N*m] (row-major) | float64 w[N] (only if has_weights)
  output : b"ALPHOUT2" | int64 M | int64 width (= max_dim + 1) |
           int64 counts[width] | int32 simplices[M*width] (-1 padded) |
           int8 dims[M] | float64 values[M]
           (the filtration in entry order; values = squared entry radius)
  stats  : JSON written next to the output (OUTPUT + ".json"), with
           alpha_time_s (read + compute), read_time_s, compute_time_s,
           write_time_s, complex sizes and QP statistics per dimension.
"""

from __future__ import annotations

import json
import os
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# The executable next to this file (alphadual.exe on Windows).
DEFAULT_EXE = os.path.join(HERE, "alphadual")
if os.path.exists(DEFAULT_EXE + ".exe"):
    DEFAULT_EXE += ".exe"


def write_input(path, X, weights=None):
    """Write a point cloud (N x m array) in alphadual's binary input format."""
    X = np.ascontiguousarray(X, dtype="<f8")
    if X.ndim != 2:
        raise ValueError("points must be an N x m array")
    N, m = X.shape
    with open(path, "wb") as fh:
        fh.write(b"ALPHAIN1")
        fh.write(np.array([N, m, 0 if weights is None else 1], dtype="<i8").tobytes())
        fh.write(X.tobytes())
        if weights is not None:
            fh.write(np.ascontiguousarray(weights, dtype="<f8").tobytes())


def read_output(path):
    """Read alphadual's output: dict with 'simplices' (M x width vertex ids,
    -1 padded; width = max_dim + 1), 'dims' (M,), 'values' (M,), 'counts'
    (simplices per dimension) and 'max_dim'. Same layout as res['ordered'] of
    alpha_dualqp.alpha_complex_dualqp. Reads formats 2 and 1 (version 1.0)."""
    with open(path, "rb") as fh:
        buf = fh.read()
    magic = buf[:8]
    if magic == b"ALPHOUT2":
        M, width = (int(v) for v in np.frombuffer(buf, dtype="<i8", count=2, offset=8))
        off = 8 + 16
        counts = np.frombuffer(buf, dtype="<i8", count=width, offset=off).tolist()
        off += 8 * width
    elif magic == b"ALPHOUT1":                      # older files: always width 3
        M, nV, nE, nT = (int(v) for v in np.frombuffer(buf, dtype="<i8", count=4, offset=8))
        width, counts, off = 3, [nV, nE, nT], 8 + 32
    else:
        raise ValueError(f"{path} is not an alphadual output file")
    simp = np.frombuffer(buf, dtype="<i4", count=width * M, offset=off).reshape(M, width)
    off += 4 * width * M
    dims = np.frombuffer(buf, dtype="i1", count=M, offset=off)
    off += M
    vals = np.frombuffer(buf, dtype="<f8", count=M, offset=off)
    return dict(simplices=simp.astype(np.int64), dims=dims.copy(), values=vals.copy(),
                counts=counts, max_dim=width - 1,
                nV=counts[0], nE=counts[1] if width > 1 else 0,
                nT=counts[2] if width > 2 else 0)


def read_stats(path):
    with open(path) as fh:
        return json.load(fh)


def run_alphadual(input_path, output_path, r=None, a1=None, solver="daqp", max_dim=2,
                  exe=DEFAULT_EXE):
    """Run the C++ program on an existing input file; returns its stats dict.
    (Timing is measured INSIDE the program, see alpha_time_s.)"""
    cmd = [exe, input_path, output_path, "--solver", solver, "--max-dim", str(int(max_dim))]
    if a1 is not None:
        cmd += ["--a1", repr(float(a1))]
    else:
        cmd += ["--r", "inf" if np.isinf(r) else repr(float(r))]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"alphadual failed ({proc.returncode}): {proc.stderr.strip()}")
    return read_stats(output_path + ".json")
