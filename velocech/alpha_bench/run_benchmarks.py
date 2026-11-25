#!/usr/bin/env python3

import csv
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import tadasets


# ============================================================
# TaDAsets generators
# ============================================================

def generate_point_cloud(shape, n, intrinsic_dim, ambient_dim, noise, seed):
    """
    Generate a point cloud using TaDAsets.

    Returns:
        numpy.ndarray of shape (n, ambient_dim)
    """

    if shape == "sphere":
        points = tadasets.sphere(
            n=n,
            noise=noise,
            ambient=ambient_dim,
            seed=seed,
        )

    elif shape == "dsphere":
        points = tadasets.dsphere(
            n=n,
            d=intrinsic_dim,
            noise=noise,
            ambient=ambient_dim,
            seed=seed,
        )

    elif shape == "torus":
        points = tadasets.torus(
            n=n,
            noise=noise,
            ambient=ambient_dim,
            seed=seed,
        )

    elif shape == "swiss_roll":
        points = tadasets.swiss_roll(
            n=n,
            noise=noise,
            ambient=ambient_dim,
            seed=seed,
        )

    elif shape == "infty_sign":
        points = tadasets.infty_sign(
            n=n,
            noise=noise,
            ambient=ambient_dim,
            seed=seed,
        )

    elif shape == "eyeglasses":
        points = tadasets.eyeglasses(
            n=n,
            noise=noise,
            ambient=ambient_dim,
            seed=seed,
        )

    else:
        raise ValueError(f"Unknown shape: {shape}")

    points = np.asarray(points)

    if points.ndim != 2:
        raise RuntimeError(
            f"TaDAsets returned an array with shape {points.shape}"
        )

    if points.shape != (n, ambient_dim):
        raise RuntimeError(
            f"Expected point cloud shape {(n, ambient_dim)}, "
            f"but TaDAsets returned {points.shape}"
        )

    return points


# ============================================================
# Input/output
# ============================================================

def write_point_cloud(points, filename):
    """
    Write point cloud in the format expected by alpha_benchmark:

        n dimension
        x11 x12 ...
        x21 x22 ...
        ...
    """

    n, dimension = points.shape

    with open(filename, "w", newline="\n") as f:
        f.write(f"{n} {dimension}\n")

        for point in points:
            f.write(" ".join(f"{x:.17g}" for x in point))
            f.write("\n")


def parse_alpha_output(output):
    """
    Parse:

        full,n,d,runtime,simplex_count

    produced by alpha_benchmark.
    """

    output = output.strip()

    if not output:
        raise RuntimeError("alpha_benchmark produced no output")

    line = output.splitlines()[-1]
    fields = [x.strip() for x in line.split(",")]

    if len(fields) != 5:
        raise RuntimeError(
            f"Unexpected alpha_benchmark output:\n{output}"
        )

    mode = fields[0]
    n = int(fields[1])
    dimension = int(fields[2])
    runtime = float(fields[3])
    simplex_count = int(fields[4])

    return {
        "benchmark_mode": mode,
        "benchmark_n": n,
        "benchmark_dimension": dimension,
        "alpha_time": runtime,
        "simplex_count": simplex_count,
    }


def run_alpha_benchmark(executable, point_cloud):
    """
    Run alpha_benchmark on a point-cloud file.
    """

    result = subprocess.run(
        [str(executable), str(point_cloud)],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            "alpha_benchmark failed.\n\n"
            f"stdout:\n{result.stdout}\n\n"
            f"stderr:\n{result.stderr}"
        )

    return parse_alpha_output(result.stdout)


# ============================================================
# Excel / CSV
# ============================================================

COLUMNS = [
    "shape",
    "intrinsic_dimension",
    "ambient_dimension",
    "num_points",
    "noise",
    "repetition",
    "seed",
    "generation_time",
    "alpha_time",
    "total_time",
    "simplex_count",
    "status",
    "error",
]


def append_result(filename, row):
    """
    Append one result to a CSV file.

    We use CSV internally because it is robust and Excel opens it
    directly. At the end, the script converts it to .xlsx if
    openpyxl is available.
    """

    filename = Path(filename)

    new_file = not filename.exists()

    with open(filename, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS)

        if new_file:
            writer.writeheader()

        writer.writerow(row)


def convert_csv_to_xlsx(csv_file, xlsx_file):
    """
    Convert the accumulated CSV results into an Excel workbook.
    """

    try:
        import openpyxl
    except ImportError:
        print()
        print("WARNING: openpyxl is not installed.")
        print("The raw CSV results were saved to:")
        print(f"  {csv_file}")
        print()
        print("Install it with:")
        print("  python -m pip install openpyxl")
        return

    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.title = "Results"

    with open(csv_file, "r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f)

        for row in reader:
            sheet.append(row)

    # Freeze header
    sheet.freeze_panes = "A2"

    # Autofilter
    sheet.auto_filter.ref = sheet.dimensions

    # Reasonable column widths
    widths = {
        "A": 15,
        "B": 20,
        "C": 20,
        "D": 15,
        "E": 12,
        "F": 12,
        "G": 15,
        "H": 20,
        "I": 20,
        "J": 20,
        "K": 18,
        "L": 12,
        "M": 50,
    }

    for column, width in widths.items():
        sheet.column_dimensions[column].width = width

    workbook.save(xlsx_file)


# ============================================================
# User input
# ============================================================

def ask(prompt, default):
    answer = input(f"{prompt} [{default}]: ").strip()
    return answer if answer else default


def ask_int(prompt, default):
    return int(ask(prompt, default))


def ask_float(prompt, default):
    return float(ask(prompt, default))


def ask_int_list(prompt, default):
    value = ask(prompt, " ".join(map(str, default)))
    return [int(x) for x in value.replace(",", " ").split()]


def choose_shape():
    shapes = [
        "sphere",
        "dsphere",
        "torus",
        "swiss_roll",
        "infty_sign",
        "eyeglasses",
    ]

    print()
    print("Available TaDAsets shapes:")

    for i, shape in enumerate(shapes, 1):
        print(f"  {i}. {shape}")

    choice = ask_int("Shape", 1)

    if choice < 1 or choice > len(shapes):
        raise ValueError("Invalid shape choice")

    return shapes[choice - 1]


def get_parameters():
    print("=" * 60)
    print("Alpha Complex Benchmark")
    print("=" * 60)

    executable = Path(
        ask(
            "Path to alpha_benchmark.exe",
            "build/alpha_benchmark.exe",
        )
    )

    if executable.is_dir():
        raise RuntimeError(
            f"{executable} is a directory, not the benchmark executable.\n"
            "Use for example:\n"
            "  ./build/alpha_benchmark.exe"
        )

    if not executable.exists():
        raise RuntimeError(
            f"Cannot find alpha_benchmark executable:\n"
            f"  {executable}"
        )

    shape = choose_shape()

    intrinsic_dim = ask_int(
        "Intrinsic dimension",
        2,
    )

    noise = ask_float(
        "Noise standard deviation",
        0,
    )

    print()
    print("Benchmark type:")
    print("  1. Increasing number of points")
    print("  2. Increasing ambient dimension")
    print("  3. Both")

    benchmark_type = ask_int(
        "Choice",
        3,
    )

    repetitions = ask_int(
        "Repetitions per configuration",
        3,
    )

    if benchmark_type in (1, 3):
        point_counts = ask_int_list(
            "Number of points",
            [100, 250, 500, 1000, 2500, 5000],
        )
    else:
        point_counts = ask_int_list(
            "Number of points",
            [1000],
        )

    if benchmark_type in (2, 3):
        ambient_dimensions = ask_int_list(
            "Ambient dimensions",
            [2, 3, 4, 5, 6],
        )
    else:
        ambient_dimensions = ask_int_list(
            "Ambient dimensions",
            [intrinsic_dim],
        )

    seed = ask_int(
        "Base random seed",
        12345,
    )

    xlsx_file = Path(
        ask(
            "Excel output file",
            "alpha_benchmark_results.xlsx",
        )
    )

    return {
        "executable": executable,
        "shape": shape,
        "intrinsic_dim": intrinsic_dim,
        "noise": noise,
        "benchmark_type": benchmark_type,
        "repetitions": repetitions,
        "point_counts": point_counts,
        "ambient_dimensions": ambient_dimensions,
        "base_seed": seed,
        "xlsx_file": xlsx_file,
    }


# ============================================================
# Benchmark
# ============================================================

def main():
    # --------------------------------------------------------
    # Verify TaDAsets
    # --------------------------------------------------------

    print()

    try:
        import tadasets
    except ImportError:
        print("ERROR: TaDAsets is not available.")
        print()
        print("Activate your gudhi-env first, then run this script.")
        print()
        print("For example:")
        print("  conda activate gudhi-env")
        print("  python run_benchmarks.py")
        sys.exit(1)

    print(f"TaDAsets: {tadasets.__file__}")

    # --------------------------------------------------------
    # Parameters
    # --------------------------------------------------------

    params = get_parameters()

    executable = params["executable"]
    shape = params["shape"]
    intrinsic_dim = params["intrinsic_dim"]
    noise = params["noise"]
    benchmark_type = params["benchmark_type"]
    repetitions = params["repetitions"]
    point_counts = params["point_counts"]
    ambient_dimensions = params["ambient_dimensions"]
    base_seed = params["base_seed"]
    xlsx_file = params["xlsx_file"]

    # --------------------------------------------------------
    # Determine configurations
    # --------------------------------------------------------

    if benchmark_type == 1:
        # Increasing n, fixed ambient dimension
        configurations = [
            (n, ambient_dimensions[0])
            for n in point_counts
        ]

    elif benchmark_type == 2:
        # Fixed n, increasing ambient dimension
        configurations = [
            (point_counts[0], ambient)
            for ambient in ambient_dimensions
        ]

    elif benchmark_type == 3:
        # Full grid
        configurations = [
            (n, ambient)
            for n in point_counts
            for ambient in ambient_dimensions
        ]

    else:
        raise ValueError("Invalid benchmark type")

    total_runs = len(configurations) * repetitions

    print()
    print("=" * 60)
    print(f"Shape:              {shape}")
    print(f"Intrinsic dimension: {intrinsic_dim}")
    print(f"Point counts:       {point_counts}")
    print(f"Dimensions:         {ambient_dimensions}")
    print(f"Noise:              {noise}")
    print(f"Repetitions:        {repetitions}")
    print(f"Total runs:         {total_runs}")
    print("=" * 60)
    print()

    # --------------------------------------------------------
    # Temporary directory for point clouds
    # --------------------------------------------------------

    with tempfile.TemporaryDirectory(
        prefix="alpha_benchmark_"
    ) as temp_dir:

        temp_dir = Path(temp_dir)

        run_number = 0

        for n, ambient_dim in configurations:

            for repetition in range(1, repetitions + 1):

                run_number += 1

                # Deterministic but different seed for every run
                seed = (
                    base_seed
                    + 1000000 * repetition
                    + 1000 * ambient_dim
                    + n
                )

                print(
                    f"[{run_number}/{total_runs}] "
                    f"n={n}, ambient={ambient_dim}, "
                    f"rep={repetition}, seed={seed}"
                )

                row = {
                    "shape": shape,
                    "intrinsic_dimension": intrinsic_dim,
                    "ambient_dimension": ambient_dim,
                    "num_points": n,
                    "noise": noise,
                    "repetition": repetition,
                    "seed": seed,
                    "generation_time": "",
                    "alpha_time": "",
                    "total_time": "",
                    "simplex_count": "",
                    "status": "OK",
                    "error": "",
                }

                point_cloud = (
                    temp_dir
                    / f"points_n{n}_d{ambient_dim}_rep{repetition}.txt"
                )

                total_start = time.perf_counter()

                try:
                    # ------------------------------------------------
                    # Generate point cloud
                    # ------------------------------------------------

                    start = time.perf_counter()

                    points = generate_point_cloud(
                        shape=shape,
                        n=n,
                        intrinsic_dim=intrinsic_dim,
                        ambient_dim=ambient_dim,
                        noise=noise,
                        seed=seed,
                    )

                    generation_time = time.perf_counter() - start

                    row["generation_time"] = generation_time

                    # ------------------------------------------------
                    # Write point cloud
                    # ------------------------------------------------

                    write_point_cloud(
                        points,
                        point_cloud,
                    )

                    # ------------------------------------------------
                    # Run Alpha Complex benchmark
                    # ------------------------------------------------

                    benchmark = run_alpha_benchmark(
                        executable,
                        point_cloud,
                    )

                    row["alpha_time"] = benchmark["alpha_time"]
                    row["simplex_count"] = benchmark["simplex_count"]

                    row["total_time"] = (
                        time.perf_counter() - total_start
                    )

                    print(
                        f"    generation: "
                        f"{generation_time:.6f} s"
                    )

                    print(
                        f"    alpha:      "
                        f"{benchmark['alpha_time']:.6f} s"
                    )

                    print(
                        f"    simplices:  "
                        f"{benchmark['simplex_count']}"
                    )

                except Exception as e:
                    row["status"] = "ERROR"
                    row["error"] = str(e)

                    row["total_time"] = (
                        time.perf_counter() - total_start
                    )

                    print()
                    print("    ERROR:")
                    print(f"    {e}")

                # ----------------------------------------------------
                # Save result immediately
                # ----------------------------------------------------

                csv_file = xlsx_file.with_suffix(".csv")

                append_result(
                    csv_file,
                    row,
                )

                print()

    # --------------------------------------------------------
    # Convert CSV -> Excel
    # --------------------------------------------------------

    csv_file = xlsx_file.with_suffix(".csv")

    print("=" * 60)
    print("Benchmark complete.")
    print()
    print(f"Raw CSV results:")
    print(f"  {csv_file}")

    convert_csv_to_xlsx(
        csv_file,
        xlsx_file,
    )

    if xlsx_file.exists():
        print()
        print(f"Excel results:")
        print(f"  {xlsx_file}")

    print("=" * 60)


if __name__ == "__main__":
    main()
