# gudhi_alpha — alpha filtration with GUDHI (any ambient dimension)

`gudhi_alpha` computes the alpha filtration of a point cloud with GUDHI's general alpha
complex (built on CGAL's d-dimensional Delaunay triangulation, so it works in **any**
ambient dimension). It has the same options, file formats and timing protocol as
`alphadual`, so results can be compared run by run (`alphadual_io.read_output` reads its
output). On every test it gives the same simplices as `alphadual` and as GUDHI's Python
interface.

Note: GUDHI always builds the **full Delaunay triangulation** first, whatever `--r` and
`--max-dim` are; its cost grows very quickly with the ambient dimension.

## Files
| | |
|---|---|
| `gudhi_alpha.cpp` | the program |
| `CMakeLists.txt` | build |
| `run_gudhi_folder.sh` | runs every input of a folder, one process per input, with an optional time limit |

Put this folder inside `alphadual_benchmark/` (or anywhere).

## Build

**Automatic:** `run_all_experiments.sh` builds `gudhi_alpha` itself with `g++` (found on the
PATH or in the usual MinGW/MSYS2 folders, or set `MINGW_BIN`), statically on Windows, and
rebuilds it when the source changes. The GUDHI headers are shipped in `third_party/gudhi`
(GUDHI 3.11.0, MIT license); CGAL, Eigen, Boost, GMP and MPFR come from the compiler's
installation. With MSYS2 (UCRT64), install them once in the MSYS2 UCRT64 terminal:
```bash
pacman -S mingw-w64-ucrt-x86_64-cgal mingw-w64-ucrt-x86_64-eigen3 mingw-w64-ucrt-x86_64-boost \
          mingw-w64-ucrt-x86_64-gmp mingw-w64-ucrt-x86_64-mpfr
```
(If something is missing, the script says which library and gives this command.)

**By hand** (alternative):
Dependencies: GUDHI headers (header-only), CGAL ≥ 5 (header-only), Eigen 3, Boost
(headers), GMP and MPFR. GUDHI headers: download a release from
https://github.com/GUDHI/gudhi-devel/releases (e.g. `gudhi.3.11.0.tar.gz`) and extract it;
use its `include` folder.

**Windows (Git Bash, Anaconda, Visual Studio Build Tools, CMake):**
```bash
conda create -n gudhi-cpp -c conda-forge cgal-cpp eigen boost-cpp gmp mpfr
conda activate gudhi-cpp
cmake -S . -B build -DCMAKE_PREFIX_PATH="$CONDA_PREFIX/Library" -DGUDHI_INCLUDE_DIR=/c/path/to/gudhi.3.11.0/include
cmake --build build --config Release
```
Run it from that activated environment (it provides the GMP/MPFR DLLs).
(If conda reports `boost-cpp` unavailable, use `libboost-headers` instead.)

**Linux:** `sudo apt install libcgal-dev libeigen3-dev libboost-dev libgmp-dev libmpfr-dev`,
then the same two `cmake` commands (without `CMAKE_PREFIX_PATH`).

The executable (`gudhi_alpha` / `gudhi_alpha.exe`) is placed in this folder;
`./gudhi_alpha --version` shows the GUDHI/CGAL versions and compiler flags.

## Use
```bash
./gudhi_alpha INPUT OUTPUT.alpha --r 0.4 --max-dim 3
./run_gudhi_folder.sh ../myrun/points --r 0.4 --max-dim 3 --time-limit 600
```
Options: `--r R` (or `inf`, or `half-diameter`: r = diameter/2, the exact diameter being
computed inside the timed part), `--max-dim D` (default 2), `--precision fast|safe|exact`
(GUDHI's three modes, default `safe` as in GUDHI), `--stats FILE.json`,
`--results-csv FILE.csv`.
Inputs: `.csv`/`.txt`/`.xyz` (one point per line, commas or spaces; a first line with only
the number of points is skipped) or `.bin` (alphadual inputs).

`run_gudhi_folder.sh` writes `results.csv` (one line per input: status `ok` /
`time_limit` / `error`, read / compute / total time, peak memory, simplex
counts per dimension) plus each input's `.alpha` output and `.alpha.json` statistics, in
`--out FOLDER` (default `gudhi_results_<date>_<time>`).

## Timing (same as alphadual)
`alpha_time_s` = from opening the input to the finished ordered filtration (writing
excluded) = `read_time_s` + `compute_time_s`. With `--r half-diameter`, `compute_time_s`
includes the diameter (also given alone as `diameter_time_s`).

## Time limit
Nothing runs inside `gudhi_alpha` besides the computation (no timer thread, no checks), so
its timing is not affected. To stop runs that take too long, the folder script's
`--time-limit SECONDS` starts each run under `timeout` (a separate, sleeping process):
the run is stopped at the limit and recorded as `time_limit` in `results.csv`. For a single
input: `timeout 600 ./gudhi_alpha INPUT OUTPUT.alpha --r 0.4` (exit code 124 = stopped).
The peak memory in the statistics is read once, after the timed part.
