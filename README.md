# Alpha filtration by dual active-set QP — benchmark package

Implementation of Algorithm 1 of E. Carlsson & J. Carlsson, *Computing the alpha complex
using dual active set quadratic programming* (arXiv:2310.00536), in **C++** and in
**Python**, with a benchmark on point clouds from
[tadasets](https://github.com/scikit-tda/tadasets), validation against GUDHI, and a
reproduction of the experiments of the paper (Section 4).

Each run computes the **alpha filtration up to `r_max`**: every simplex of dimension
0..d (d chosen by you; 2 = up to triangles) entering at a radius ≤ `r_max`, with its entry
value (squared radius), ordered by entry value (faces before cofaces).
`r_max = inf` gives the complete filtration.

## What is in this folder

| File | Used for |
|---|---|
| `run_benchmark.sh` | **Main benchmark (C++)**: build → generate inputs → timed runs → Excel |
| `make_inputs.py`, `collect_results.py`, `alphadual_io.py` | the Python steps of `run_benchmark.sh` (inputs, results, file I/O) |
| `alphadual.cpp`, `CMakeLists.txt`, `Makefile`, `third_party/` | the C++ program `alphadual` (with DAQP and nanoflann; built on your machine) |
| `run_experiments.py` | the same benchmark with the **Python** implementation (also holds the shared questions and Excel writer) |
| `alpha_dualqp.py` | the Python implementation |
| `summary_to_latex.py` | Excel results → **LaTeX** table |
| `validate_against_gudhi.py` | check against **GUDHI** (Python implementation, or C++ with `--cpp`) |
| `compare_cpp_python.py` | check that C++ and Python give the same filtrations |
| `reproduce_paper_section4.py` | re-run the experiments of **the paper** (Sections 4.1–4.3) |
| `data/` | the paper's data sets (Stanford bunny, cyclooctane, optical patches) |
| `run_all_experiments.sh`, `measure_run.py`, `experiments_to_excel.py` | **all experiments**: every algorithm on every input subfolder, with a time limit, wall time and peak memory measured from outside, one Excel summary per subfolder |
| `alpha_at_half_diameter.py` | alpha filtration at r = diameter/2 (timed: diameter + alpha) |
| `gudhi_alpha/` | GUDHI's alpha complex in any dimension (C++; see its README) |
| `velocech/` | the velocech repository (Cech H1 barcode, `H1_skip.cpp`), unchanged |
| `velocech_build/` | builds `velocech/H1_skip.cpp` for `run_all_experiments.sh` (automatically) |
| `requirements.txt` | Python packages |

All scripts must stay together in this folder (they use each other).

## Setup

```bash
python -m pip install -r requirements.txt   # numpy scipy numba daqp tadasets openpyxl (+ gudhi for validation)
```
The C++ program is compiled automatically by `run_benchmark.sh` — with `make` if available,
otherwise with `cmake` (any compiler CMake finds: Visual C++/Build Tools, MinGW, GCC,
Clang). To build it by hand: `make`, or
`cmake -S . -B build-cmake -DCMAKE_BUILD_TYPE=Release && cmake --build build-cmake --config Release`.
Only source code is shipped.

## 1. Benchmark (C++) — main use

```bash
./run_benchmark.sh                 # asks the questions, creates the folder bench_<date>_<time>/
./run_benchmark.sh myrun           # same, in the folder myrun/
./run_benchmark.sh myrun2 --config myrun/config.json   # re-run a saved configuration
```
Questions: tadasets shape and its parameters, number(s) of points, ambient dimension(s),
`r_max`, maximum simplex dimension d, runs per configuration, seed, solver (`daqp` = the
solver used by the paper's authors, `native` = our own, faster). Lists may be typed as
`500,1000` or ranges `500:5000:500`; every combination is run.

The work folder then contains:

| | |
|---|---|
| `<name>.xlsx` | **results**: Inputs, Runs (one row per run), Summary (mean ± std), Chart |
| `points/` | **the generated point clouds as CSV** + `index.csv` (see section 3) |
| `inputs/`, `outputs/` | binary input / output files of `alphadual` (ordered filtrations) |
| `config.json` | the configuration (re-usable with `--config`) |

**Timing**: measured inside the C++ program, from opening the input file to the finished
ordered filtration (writing excluded). The Excel file has three columns:
`alpha_complex_time_s` = reading + computing, `input_read_time_s` = reading alone (a
negligible part), `compute_time_s` = the computation alone; the Summary sheet gives mean
and standard deviation of each. Nothing is timed inside the computation.

**Crash safety**: results are saved after every run. To continue an interrupted benchmark:
`./run_benchmark.sh myrun --resume` (add `--retry-failed` to re-run failed jobs, e.g. out of
memory). If the Excel file is open in Excel while it is updated, the latest results go to
`<name>_latest.xlsx`.

## 2. Benchmark (Python implementation)

```bash
python run_experiments.py                    # same questions, same Excel report
python run_experiments.py --config <name>_config.json --resume   # after a crash
```
Writes `<name>.xlsx`, `<name>_points/` (point clouds as CSV) and `<name>_config.json`.

## 3. Saved point clouds (to run your own implementation on the same inputs)

Every point cloud generated with tadasets is saved as `<shape>_<points>_<dimension>_<run>.csv`
(e.g. `torus_1000_3_1.csv`, or `2sphere_1000_5_1.csv` for `dsphere` with d = 2; runs with
several `r_max` values share the same cloud): **one point per
line, coordinates separated by commas**, written with 17 significant digits, so the
numbers read back are *exactly* the same doubles. `index.csv` lists every file with its
configuration: shape, shape parameters, number of points, ambient dimension, `r_max`, d
and seed (configuration and run numbers match the rows of the Excel file).

To generate the point clouds **without** running the benchmark:
`python make_inputs.py myrun` (asks the same questions; writes `myrun/points/`).

For a fair speed comparison note what our time includes (reading a binary input file +
computing the ordered filtration up to dimension d) and use the same compiler flags
(`./alphadual --version` shows ours; they are also in the Excel Inputs sheet).

## 4. LaTeX table

```bash
python summary_to_latex.py myrun/<name>.xlsx      # writes <name>_summary.tex
```
Options: `--columns n,ambient,r,runs,time,compute,edges,triangles,dim3` (`compute` =
computation alone, C++ results), `--unit ms`,
`--standalone` (compilable document). Works on partial results too.

## 5. Validation

```bash
python validate_against_gudhi.py          # Python implementation vs GUDHI (35 cases)
python validate_against_gudhi.py --cpp    # C++ program vs GUDHI (needs the program built)
python compare_cpp_python.py              # C++ vs Python (needs ~4 GB RAM)
```
GUDHI computes the full Delaunay triangulation, so it is only usable in low ambient
dimension (≤ 4–5); the tests cover d = 0..5, weighted points and complete filtrations.

## 6. Paper reproduction (Sections 4.1–4.3)

```bash
python reproduce_paper_section4.py --cpp --max-dim 3      # C++, up to tetrahedra (as the paper)
python reproduce_paper_section4.py                        # Python implementation, d = 2
```
Options: `--only 4.1|4.2|4.3`, `--solvers native,daqp`, `--repeats`, `--seeds`,
`--check-gudhi`, `--betti`. Results are also appended to a CSV file run by run.
The bunny is reproduced with r = 0.01 (this gives exactly the paper's counts).

## 7. All experiments (alphadual, half diameter, GUDHI, velocech)

```bash
python make_inputs.py inputs/sphere_noise005      # once per shape / noise: the point clouds
bash run_all_experiments.sh                       # asks the folder (e.g. inputs), time limit, ...
```
Each algorithm runs on every input of every subfolder under the time limit; larger inputs of a
combination (shape, parameters, ambient dimension) are skipped after a failure. velocech is built
automatically from `velocech/` (CMake; with Visual C++ without OpenMP), its inputs are written in
its text format during the untimed preparation, and its own times (distances + H1) are reported as
`time_s`, with `distance_time_s`, `h1_time_s` and `n_h1_bars`. It computes the full Cech H1
barcode (r_max and the maximum dimension do not apply); its bars are on a diameter scale
(= 2 x the alpha radii: on generic inputs they equal GUDHI's alpha H1 bars times 2).

## Notes

* **Windows / Git Bash**: Python must be found by Git Bash (with Anaconda, run
  `conda init bash` in the Anaconda Prompt once, then reopen Git Bash). Start the script
  from an open terminal; on any error it shows a box with the step and the reason.
* **Degenerate inputs**: points lying *exactly* on a sphere (`dsphere`/`sphere` without
  noise) make every set of points a simplex once `r_max` ≥ the sphere's radius; add a little
  noise or use a smaller `r_max`.
* `r_max = inf` (complete filtration) costs roughly n³ and needs the `daqp` solver.
* Third-party code: DAQP v0.10.3 (MIT) and nanoflann 1.14.0 (BSD), licenses in `third_party/`.
