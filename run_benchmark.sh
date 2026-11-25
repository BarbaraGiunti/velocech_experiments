#!/usr/bin/env bash
# ===========================================================================
# run_benchmark.sh -- C++ alpha complex benchmark, end to end.
#
#   ./run_benchmark.sh [WORKDIR] [--config CONFIG.json]
#   ./run_benchmark.sh WORKDIR --reuse-inputs [--solver daqp|native]
#   ./run_benchmark.sh WORKDIR --resume [--retry-failed]
#
#   The Excel report is updated after EVERY job, so it is usable even if the
#   benchmark stops early. Each job's results are written by alphadual itself
#   as soon as it finishes. A job that crashes (e.g. out of memory) is
#   recorded as failed and the benchmark continues with the next one.
#
#   --resume       : continue an interrupted benchmark in WORKDIR: finished
#                    and failed jobs are skipped (inputs are reused)
#   --retry-failed : with --resume, run the failed jobs again
#
#   --reuse-inputs : skip step 2 and time alphadual again on the inputs that
#                    already exist in WORKDIR (e.g. to compare solvers or
#                    builds on identical point clouds); --solver overrides
#                    the solver stored in WORKDIR/jobs.tsv.
#
#   1. builds the C++ program alphadual (make, or cmake if make is missing)
#   2. make_inputs.py   : asks the usual questions (or reads --config),
#                         generates every point cloud, writes binary inputs
#   3. alphadual        : run once per job; the alpha complex is timed INSIDE
#                         the program (reading the input + computing the
#                         ordered filtration; writing the output excluded)
#   4. collect_results.py: statistics + outputs -> Excel report
#
# WORKDIR defaults to bench_<date>_<time>. Everything (inputs, outputs, logs,
# Excel file) is written there. Jobs run one at a time, never in parallel,
# so that timings are not disturbed by each other.
# ===========================================================================
set -Eeuo pipefail

# ---- Error reporting --------------------------------------------------------
# With "set -e" any failing command stops the script, and some fail without
# printing anything. So the script records the step it is in, the command
# that failed (ERR trap) and any explanation (die), and repeats them at the
# end, just before exiting -- where they cannot scroll away. In an
# interactive window it then waits for Enter (Git Bash closes the window of
# a script started by double-click as soon as it ends).
STEP="starting"
FAILED_CMD=""
FAILED_LINE=""
DIE_MSG=""
trap 'FAILED_CMD=$BASH_COMMAND; FAILED_LINE=$LINENO' ERR
on_exit() {
    local status=$?
    if [[ $status -ne 0 ]]; then
        {
            echo
            echo "=================================================================="
            echo "run_benchmark.sh stopped with an error (exit code $status)"
            echo "  during step : $STEP"
            if [[ -n "$DIE_MSG" ]]; then
                echo "  reason      : $DIE_MSG"
            elif [[ -n "$FAILED_CMD" ]]; then
                echo "  command     : $FAILED_CMD"
                echo "  at line     : $FAILED_LINE of run_benchmark.sh"
                echo "  (any message printed by that command is just above this box)"
            fi
            echo "  For a full trace run:  bash -x ./run_benchmark.sh 2>&1 | tee trace.log"
            echo "=================================================================="
        } >&2
        if [[ -t 0 ]]; then
            read -r -p "Press Enter to close..." _ || true
        fi
    fi
}
trap on_exit EXIT
die() { DIE_MSG="$*"; echo; echo "ERROR: $*" >&2; exit 1; }

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- Prerequisites ---------------------------------------------------------
STEP="checking Python"
# Python 3: $PYTHON if set, else the first of python3 / python / py that is a
# working Python 3 (on Windows "python3" is often a Microsoft Store stub).
find_python() {
    local c
    for c in ${PYTHON:-} python3 python py; do
        if command -v "$c" >/dev/null 2>&1 &&
           "$c" -c "import sys; sys.exit(0 if sys.version_info[0] == 3 else 1)" >/dev/null 2>&1; then
            echo "$c"; return 0
        fi
    done
    return 1
}
PYTHON="$(find_python)" || die "Python 3 not found (tried: ${PYTHON:-} python3 python py).
Install Python 3 and make sure it is on the PATH, or run:  PYTHON=/path/to/python ./run_benchmark.sh"
export PYTHONUNBUFFERED=1    # show Python's prompts immediately (needed in Git Bash)
STEP="checking the Python packages (numpy, scipy, tadasets, openpyxl) with $PYTHON"
"$PYTHON" -c "import numpy, scipy, tadasets, openpyxl" >/dev/null 2>&1 ||
    die "missing Python packages. Install them with:
    $PYTHON -m pip install -r \"$HERE/requirements.txt\""

WORK=""
EXTRA=()
REUSE=""
RESUME=""
RETRY_FAILED=""
SOLVER_OVERRIDE=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) EXTRA+=(--config "$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"); shift 2 ;;
        --reuse-inputs) REUSE=1; shift ;;
        --resume) RESUME=1; REUSE=1; shift ;;
        --retry-failed) RETRY_FAILED=1; shift ;;
        --solver) SOLVER_OVERRIDE="$2"; shift 2 ;;
        -h|--help) sed -n '2,38p' "$0"; exit 0 ;;
        *) WORK="$1"; shift ;;
    esac
done
WORK="${WORK:-bench_$(date +%Y%m%d_%H%M%S)}"

# 1. build (no-op if up to date). On Windows the program is alphadual.exe.
STEP="checking that all files of the package are present"
missing=()
for f in alphadual.cpp CMakeLists.txt Makefile third_party/daqp/src/daqp.c \
         third_party/daqp/include/api.h third_party/nanoflann/nanoflann.hpp \
         make_inputs.py collect_results.py run_experiments.py alphadual_io.py \
         alpha_dualqp.py requirements.txt; do
    [[ -e "$HERE/$f" ]] || missing+=("$f")
done
[[ ${#missing[@]} -eq 0 ]] || die "files of the package are missing in $HERE:
    ${missing[*]}
Unzip the COMPLETE alphadual_benchmark.zip into one folder and run the script there."
STEP="building the C++ program"
have_cxx() { command -v "${CXX:-g++}" >/dev/null 2>&1 || command -v c++ >/dev/null 2>&1 ||
             command -v clang++ >/dev/null 2>&1; }
exe_path() { if [[ -f "$HERE/alphadual.exe" ]]; then echo "$HERE/alphadual.exe";
             else echo "$HERE/alphadual"; fi; }
# Build tool: make (with the Makefile) if available, otherwise CMake (with
# CMakeLists.txt, using whatever compiler CMake finds: GCC, Clang, MinGW or
# Visual C++). BUILD=make or BUILD=cmake forces one of them.
BUILD="${BUILD:-}"
if [[ -z "$BUILD" ]]; then
    if command -v make >/dev/null 2>&1 && have_cxx; then BUILD=make
    elif command -v cmake >/dev/null 2>&1; then BUILD=cmake
    fi
fi
if [[ "$BUILD" == make ]]; then
    echo "== Building alphadual (make)"
    make -C "$HERE" --no-print-directory || die "the build failed (see the compiler messages above)"
elif [[ "$BUILD" == cmake ]]; then
    echo "== Building alphadual (cmake)"
    BDIR="$HERE/build-cmake"
    CLOG="$HERE/cmake-configure.log"
    BLOG="$HERE/cmake-build.log"
    # Configure with CMake's default generator (or $CMAKE_GENERATOR if set),
    # showing CMake's output and keeping a copy in cmake-configure.log.
    cmake_try() {    # $1 = generator name, "" = default
        rm -rf "$BDIR"                  # a cache from another generator would clash
        if [[ -n "$1" ]]; then
            echo "-- trying the CMake generator \"$1\""
            cmake -S "$HERE" -B "$BDIR" -G "$1" -DCMAKE_BUILD_TYPE=Release 2>&1 | tee "$CLOG"
        else
            cmake -S "$HERE" -B "$BDIR" -DCMAKE_BUILD_TYPE=Release 2>&1 | tee "$CLOG"
        fi
    }
    # The lines of CMake's output that explain an error.
    cmake_errors() { grep -A4 "CMake Error" "$CLOG" 2>/dev/null | head -20 || true; }
    # Toolchain problems (no compiler or build program for the default
    # generator), as opposed to errors in the project itself.
    toolchain_problem() {
        grep -qiE "No CMAKE_(C|CXX)_COMPILER could be found|could not find any instance of Visual Studio|is not able to compile a simple test program|CMAKE_(C|CXX)_COMPILER not set|CMAKE_MAKE_PROGRAM is not set|unable to find a build program|generator doesn't exist|Could not create named generator" "$CLOG" 2>/dev/null
    }
    configured=""
    if [[ -f "$BDIR/CMakeCache.txt" ]] &&
       cmake -S "$HERE" -B "$BDIR" -DCMAKE_BUILD_TYPE=Release >/dev/null 2>&1; then
        configured=1                    # an earlier configuration still works
    elif cmake_try ""; then
        configured=1
    elif toolchain_problem; then
        # Only when the default generator has no compiler or build program:
        # try the other generators available on this system.
        for gen in "Ninja" "MinGW Makefiles" "MSYS Makefiles" "Unix Makefiles"; do
            case "$gen" in
                Ninja)             command -v ninja        >/dev/null 2>&1 || continue ;;
                "MinGW Makefiles") command -v mingw32-make >/dev/null 2>&1 || continue ;;
                *)                 command -v make         >/dev/null 2>&1 || continue ;;
            esac
            if cmake_try "$gen"; then configured=1; break; fi
        done
    fi
    [[ -n "$configured" ]] || die "CMake could not configure the build. CMake said:
$(cmake_errors)
(full CMake output: $CLOG)"
    cmake --build "$BDIR" --config Release 2>&1 | tee "$BLOG" ||
        die "the C++ compilation failed. First error lines:
$(grep -iE "error" "$BLOG" | head -12)
(full compiler output: $BLOG)"
elif [[ -f "$(exe_path)" ]]; then
    echo "== no make or cmake found: using the existing $(basename "$(exe_path)")"
else
    die "cannot build the C++ program: neither 'make' (with g++) nor 'cmake' was found.
Install CMake (https://cmake.org) with a C++ compiler (e.g. Visual Studio
Build Tools, or MinGW-w64), or use WSL / MSYS2 (see README.md)."
fi
EXE="$(exe_path)"
[[ -f "$EXE" ]] || die "the build did not produce $EXE"

# 2. inputs
STEP="generating the inputs (make_inputs.py)"
if [[ -n "$REUSE" ]]; then
    [[ -f "$WORK/jobs.tsv" ]] || { echo "no $WORK/jobs.tsv to reuse" >&2; exit 1; }
    echo "== Reusing the inputs in $WORK"
else
    echo "== Generating inputs in $WORK"
    "$PYTHON" "$HERE/make_inputs.py" "$WORK" "${EXTRA[@]+"${EXTRA[@]}"}"
fi
WORK="$(cd "$WORK" && pwd)"
if [[ -n "$SOLVER_OVERRIDE" ]]; then
    # record the new solver for the jobs and for the Excel report
    awk -v s="$SOLVER_OVERRIDE" 'BEGIN{FS=OFS="\t"} {$5=s; print}' "$WORK/jobs.tsv" > "$WORK/jobs.tmp"
    mv "$WORK/jobs.tmp" "$WORK/jobs.tsv"
    "$PYTHON" - "$WORK/manifest.json" "$SOLVER_OVERRIDE" <<'PY'
import json, sys
p, s = sys.argv[1], sys.argv[2]
m = json.load(open(p)); m["config"]["solver"] = s
json.dump(m, open(p, "w"), indent=2)
PY
fi

STEP="running the jobs"
# 3. one alphadual run per job (jobs.tsv: job, input, output, r, solver, max_dim;
#    job lists written before max_dim existed have 5 columns: d = 2)
total=$(grep -c . "$WORK/jobs.tsv" || true)
[[ "$total" -gt 0 ]] || die "the job list $WORK/jobs.tsv is empty (no runs were configured?)"
echo "== Running $total job(s)"
k=0
while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%$'\r'}"                    # tolerate Windows line endings
    [[ -z "$line" ]] && continue
    IFS=$'\t' read -r job inp out r solver maxdim <<< "$line"
    maxdim="${maxdim:-2}"
    k=$((k + 1))
    if [[ -n "$RESUME" ]]; then
        if [[ -f "$WORK/$out.json" ]]; then
            printf '[%d/%d] %s  already done\n' "$k" "$total" "$job"; continue
        fi
        if [[ -f "$WORK/$out.error" && -z "$RETRY_FAILED" ]]; then
            printf '[%d/%d] %s  failed earlier (use --retry-failed to run it again)\n' \
                "$k" "$total" "$job"; continue
        fi
    fi
    rm -f "$WORK/$out.error" "$WORK/$out.json"
    if "$EXE" "$WORK/$inp" "$WORK/$out" --r "$r" --solver "$solver" --max-dim "$maxdim" \
            2> "$WORK/$out.stderr"; then
        t=$(grep -o '"alpha_time_s": [0-9.eE+-]*' "$WORK/$out.json" | awk '{print $2}')
        printf '[%d/%d] %s  r=%s  d=%s  %s: %.4f s\n' "$k" "$total" "$job" "$r" "$maxdim" "$solver" "$t"
        rm -f "$WORK/$out.stderr"
    else
        # Keep going: the failure is recorded and shown in the Excel file.
        { echo "exit code $?"; cat "$WORK/$out.stderr"; } > "$WORK/$out.error"
        rm -f "$WORK/$out.stderr" "$WORK/$out.json"
        printf '[%d/%d] %s  FAILED (see %s)\n' "$k" "$total" "$job" "$out.error"
    fi
    # Save the Excel report now (outside alphadual: never part of a timing).
    "$PYTHON" "$HERE/collect_results.py" "$WORK" --partial \
        || echo "  (warning: the Excel report was not updated this time, see the reason above;
     no results are lost: they are kept on disk and included in the next update)"
done < "$WORK/jobs.tsv"

# 4. Excel report
STEP="writing the Excel report (collect_results.py)"
echo "== Collecting results"
"$PYTHON" "$HERE/collect_results.py" "$WORK"
