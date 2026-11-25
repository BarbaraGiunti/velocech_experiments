#!/usr/bin/env bash
# ===========================================================================
# run_all_experiments.sh -- run several algorithms on all input point clouds.
#
#   ./run_all_experiments.sh                        # asks for everything
#   ./run_all_experiments.sh --inputs DIR --time-limit 600 [options]
#   ./run_all_experiments.sh --out RESULTS --resume # continue an interrupted run
#
# The input folder DIR contains SUBFOLDERS of point clouds (*.csv *.txt *.xyz
# *.bin; index.csv is skipped). For every subfolder, every algorithm is run on
# every input of that subfolder (algorithm by algorithm), each run under the
# time limit, and an Excel summary is written for the subfolder:
#
#   RESULTS/settings.json                 settings of the whole experiment
#   RESULTS/<subfolder>/summary.xlsx      Runs, Summary, Time by input, Memory by input
#   RESULTS/<subfolder>/runs.jsonl        one record per run (status, wall time, peak memory)
#   RESULTS/<subfolder>/reported.tsv      values reported by the codes themselves
#   RESULTS/<subfolder>/outputs/<algorithm>/<input>.*   outputs, stdout, stderr
#
# For every run, measure_run.py measures from OUTSIDE the program its wall
# time and peak memory (it only waits for the program: nothing runs inside or
# alongside it), and stops it at the time limit. The codes' own times
# (computation, reading excluded) are collected too, as time_s.
# Results are saved after every run; --resume skips runs already done.
#
# Options (asked interactively when missing):
#   --inputs DIR        folder containing the input subfolders
#   --time-limit S      seconds per run (kill switch)
#   --r R               r_max for alphadual and gudhi (a number or inf)
#   --max-dim D         largest simplex dimension (default 2)
#   --solver S          alphadual / half_diameter solver: daqp or native (default daqp)
#   --precision P       gudhi: fast, safe or exact (default safe)
#   --algorithms LIST   comma-separated, default: all of ALL_ALGORITHMS below
#                       (alphadual, half_diameter, gudhi, velocech)
#   --velocech-threads N  OpenMP threads for velocech (default: all cores)
#   --out DIR           results folder (default experiments_<date>_<time>)
#   --resume            continue in --out DIR (settings are read from it)
#   --skip-larger MODE  failure (default): when an algorithm hits the time limit
#                       or fails on an input, it skips the LARGER inputs of the
#                       same group (same shape, shape parameters such as noise,
#                       and ambient dimension); timeout: only after a timeout;
#                       never: run everything. Skipped runs appear in the Excel
#                       file with status "skipped" and the reason. Groups are
#                       read from the index.csv next to the inputs (written by
#                       the benchmark), else from file names
#                       <shape>_<points>_<dim>_<run> within each subfolder.
#                       Inputs run group by group in increasing size.
#
# TO ADD ANOTHER CODE: see "ALGORITHMS" below (two small functions + its name).
# ===========================================================================
set -uo pipefail     # no "set -e": one failing run must not stop the others
set -E               # the ERR trap below also applies inside functions

# ---- Error reporting --------------------------------------------------------
# If the script stops, a box at the very end says where (STEP) and why
# (DIE_MSG, or the last command that failed), cleaned of carriage returns
# (Windows programs' error output contains them; in a terminal they make text
# overwrite itself). The box is also saved to run_all_last_error.txt in this
# folder, so it can be read or sent even after the window is closed.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STEP="starting"
rm -f "$HERE/run_all_last_error.txt" 2>/dev/null     # from an earlier run
DIE_MSG=""
FAILED_CMD=""
FAILED_LINE=""
trap 'FAILED_CMD=$BASH_COMMAND; FAILED_LINE=$LINENO' ERR
on_exit() {          # keep the window open at the end (success or error)
    local status=$?
    # Output was copied to the log by tee: switch back to the terminal and let
    # tee write everything out first (at most 5 s, never hang here).
    if [[ -n "${TEE_PID:-}" ]]; then
        exec 1>&3 2>&4
        local w
        for w in 1 2 3 4 5 6 7 8 9 10; do
            kill -0 "$TEE_PID" 2>/dev/null || break
            sleep 0.5
        done
    fi
    if [[ $status -ne 0 && -z "${STOPPED:-}" ]]; then
        local box
        box="$(
            echo "=================================================================="
            echo "run_all_experiments.sh stopped with an error (exit code $status)"
            echo "  during step: $STEP"
            if [[ -n "$DIE_MSG" ]]; then
                echo "  reason:"
                printf '%s\n' "$DIE_MSG" | sed 's/^/    /'
            elif [[ -n "$FAILED_CMD" ]]; then
                echo "  last failed command (line $FAILED_LINE): $FAILED_CMD"
            fi
            echo "  (this message is also in: $HERE/run_all_last_error.txt)"
            echo "=================================================================="
        )"
        box="$(printf '%s' "$box" | tr -d '\r')"
        printf '\n%s\n' "$box" >&2
        printf '%s\n' "$box" > "$HERE/run_all_last_error.txt" 2>/dev/null
        [[ -n "${OUT:-}" && -d "${OUT:-}" ]] && printf '%s\n' "$box" > "$OUT/run_all_last_error.txt" 2>/dev/null
    fi
    # wait for Enter in an interactive window, unless stopped by Ctrl-C / a signal
    if [[ -t 0 && -z "${STOPPED:-}" ]]; then
        printf '%s' "Press Enter to close..." >&2
        read -r _ || true
    fi
}
trap on_exit EXIT
trap 'STOPPED=1; echo; echo "Stopped (Ctrl-C). Continue later with --out <results folder> --resume"; exit 130' INT
trap 'STOPPED=1; exit 143' TERM
die() {              # stop with a reason (shown now, and again in the final box)
    DIE_MSG="$(printf '%s' "$*" | tr -d '\r')"
    printf '\nERROR: %s\n' "$DIE_MSG" >&2
    exit 1
}
exe() { if [[ -f "$1.exe" ]]; then echo "$1.exe"; else echo "$1"; fi; }   # Windows: .exe

STEP="checking Python"
# Python 3 (python3 / python / py, or $PYTHON)
PY=""
for c in ${PYTHON:-} python3 python py; do
    if command -v "$c" >/dev/null 2>&1 &&
       "$c" -c "import sys; sys.exit(0 if sys.version_info[0] == 3 else 1)" >/dev/null 2>&1; then
        PY="$c"; break
    fi
done
[[ -n "$PY" ]] || die "Python 3 not found (tried python3, python, py)"
"$PY" -c "import numpy, scipy, openpyxl" >/dev/null 2>&1 ||
    die "missing Python packages: $PY -m pip install -r \"$HERE/requirements.txt\""
export PYTHONUNBUFFERED=1

# ===========================================================================
# ALGORITHMS
#
# Each algorithm NAME is defined by:
#   cmd_NAME INPUT PREFIX     sets the array CMD = the command running the code
#                             on ONE input file; the code should write its
#                             outputs to files starting with PREFIX. (Its
#                             stdout/stderr are saved as PREFIX.stdout/.stderr.)
#   report_NAME PREFIX        (optional) prints "key=value" lines read from the
#                             code's own outputs after the run, e.g.
#                               time_s=0.123     the code's own computation time
#                               n_dim0=1000      counts, or any other value
#                               status=error: .. if the code reports a failure
#                             Every key becomes a column of the Excel file.
#   check_NAME                (optional) fails with a message if the code is
#                             not available (checked before starting).
# and its name in ALL_ALGORITHMS (the order is the order of the runs).
# Settings available in the functions: $R $MAXDIM $SOLVER $PRECISION $PY $HERE
# $BIN_DIR (the inputs converted to alphadual's .bin format).
# ===========================================================================
ALL_ALGORITHMS=(alphadual half_diameter gudhi velocech)

ALPHADUAL_EXE="$(exe "$HERE/alphadual")"
GUDHI_EXE_GIVEN="${GUDHI_EXE:-}"             # set by the user: used as is, never rebuilt
GUDHI_DIR="$HERE/gudhi_alpha"
GUDHI_EXE="${GUDHI_EXE:-$(exe "$GUDHI_DIR/gudhi_alpha")}"

# helpers for the report functions: read values from a flat JSON stats file
json_get() { grep -o "\"$2\": [^,}]*" "$1" 2>/dev/null | head -1 | sed 's/^[^:]*: *//; s/"//g'; }
json_counts() {      # "n_simplices": [a, b, c]  ->  n_dim0=a n_dim1=b ...
    grep -o '"n_simplices": \[[^]]*\]' "$1" 2>/dev/null | sed 's/.*\[//; s/\]//' |
        tr ',' '\n' | awk '{gsub(/ /, ""); if ($0 != "") print "n_dim" NR-1 "=" $0}'
}
json_report() {      # time_s from KEY, plus read/compute/total times and counts
    local f="$1" key="$2" v k
    [[ -f "$f" ]] || return 0
    v="$(json_get "$f" "$key")"; [[ -n "$v" ]] && echo "time_s=$v"
    for k in read_time_s compute_time_s alpha_time_s diameter_time_s; do
        v="$(json_get "$f" "$k")"; [[ -n "$v" ]] && echo "$k=$v"
    done
    json_counts "$f"
}

# --- alphadual: the C++ program, at r_max = $R (reads the .bin input) ------
check_alphadual() {
    [[ -f "$ALPHADUAL_EXE" ]] || die "alphadual is not built ($ALPHADUAL_EXE): run ./run_benchmark.sh once, or make / cmake (README.md)"
}
cmd_alphadual() {
    CMD=("$ALPHADUAL_EXE" "$BIN_DIR/$(basename "$1").bin" "$2.alpha" --r "$R"
         --max-dim "$MAXDIM" --solver "$SOLVER" --stats "$2.alpha.json")
}
report_alphadual() { json_report "$1.alpha.json" compute_time_s; }

# --- half_diameter: diameter, then alphadual at r = D/2 ---------------------
# Uses the functions of alpha_at_half_diameter.py directly (reading, exact
# diameter, alphadual at D/2), without that script's own CSV/Excel output, so
# that the run is not slowed down by work unrelated to the algorithm. time_s =
# diameter + alpha construction (reading excluded), measured as in that script
# (the diameter routine is warmed up first).
HALF_DIAMETER_PY='
import json, sys, tempfile
inp, prefix, maxdim, solver, here = sys.argv[1:6]
sys.path.insert(0, here)
import numpy as np
import alpha_at_half_diameter as h
X = h.read_points(inp)                                  # not timed
rng = np.random.default_rng(0)                          # warm-up, not timed
h.diameter(rng.random((200, 3))); h.diameter(rng.random((200, 5)))
with tempfile.TemporaryDirectory() as tmp:
    D, t_diam, t_alpha, counts, deg, failed = h.run_cpp(X, solver, int(maxdim), tmp)
json.dump(dict(time_s=t_diam + t_alpha, diameter_time_s=t_diam, alpha_time_s=t_alpha,
               diameter=D, r=D / 2, n_simplices=counts, failed_QPs=failed),
          open(prefix + ".json", "w"))
'
check_half_diameter() { check_alphadual; [[ -f "$HERE/alpha_at_half_diameter.py" ]] || die "alpha_at_half_diameter.py missing"; }
cmd_half_diameter() {
    CMD=("$PY" -c "$HALF_DIAMETER_PY" "$1" "$2" "$MAXDIM" "$SOLVER" "$HERE")
}
report_half_diameter() {
    json_report "$1.json" time_s
    local v k
    for k in diameter r failed_QPs; do v="$(json_get "$1.json" "$k")"; [[ -n "$v" ]] && echo "$k=$v"; done
}

# --- gudhi: gudhi_alpha (GUDHI's alpha complex), at r_max = $R -------------
# Built automatically (like velocech) with g++ -- found on the PATH or in the
# usual MinGW folders, or MINGW_BIN -- from gudhi_alpha/gudhi_alpha.cpp, the
# GUDHI headers shipped in gudhi_alpha/third_party, and CGAL, Eigen, Boost,
# GMP and MPFR from the compiler's installation (e.g. MSYS2: pacman, see below).
# Statically linked on Windows. Rebuilt when the source changes.
check_gudhi() {
    local src="$GUDHI_DIR/gudhi_alpha.cpp" info="$GUDHI_DIR/build_info.txt" log="$GUDHI_DIR/build.log"
    if [[ -n "$GUDHI_EXE_GIVEN" ]]; then
        [[ -f "$GUDHI_EXE" ]] || die "GUDHI_EXE=$GUDHI_EXE does not exist"
        echo "   gudhi: $GUDHI_EXE (given, not rebuilt)"; gudhi_startup_test; return
    fi
    [[ -f "$src" ]] || die "$src not found (the gudhi_alpha folder belongs in $HERE)"
    [[ -n "$GXX" ]] || GXX="$(find_gxx)"
    if [[ -n "$GXX" ]]; then
        local need="" static=""
        is_windows && static="-static"
        [[ -f "$GUDHI_EXE" ]] || need=1
        grep -q "GNU" "$info" 2>/dev/null || need=1
        [[ -n "$static" ]] && ! grep -q "static" "$info" 2>/dev/null && need=1
        [[ "$src" -nt "$GUDHI_EXE" ]] && need=1
        if [[ -n "$need" ]]; then
            local gpath gprefix env_name pkg d h missing=() incs=(-I"$GUDHI_DIR/third_party/gudhi/include")
            gpath="$(dirname "$GXX"):$PATH"
            gprefix="$(dirname "$(dirname "$GXX")")"           # e.g. /c/msys64/ucrt64
            for d in "${EIGEN_DIR:-}" "$gprefix/include/eigen3" /usr/include/eigen3 \
                     /usr/local/include/eigen3 /opt/homebrew/include/eigen3; do
                [[ -n "$d" && -f "$d/Eigen/Core" ]] && { incs+=(-I"$d"); break; }
            done
            # the libraries must be visible to this compiler
            for h in "CGAL/Epick_d.h:CGAL" "Eigen/Core:Eigen" "boost/version.hpp:Boost" "gmp.h:GMP" "mpfr.h:MPFR"; do
                printf '#include <%s>\n' "${h%%:*}" |
                    PATH="$gpath" "$GXX" -std=c++17 -E -x c++ "${incs[@]}" - > /dev/null 2>&1 || missing+=("${h##*:}")
            done
            if [[ ${#missing[@]} -gt 0 ]]; then
                env_name="$(basename "$gprefix")"
                case "$env_name" in
                    ucrt64) pkg="mingw-w64-ucrt-x86_64" ;; mingw64) pkg="mingw-w64-x86_64" ;;
                    clang64) pkg="mingw-w64-clang-x86_64" ;; *) pkg="" ;;
                esac
                if [[ -n "$pkg" ]]; then
                    die "cannot build gudhi_alpha: these libraries are missing for $GXX: ${missing[*]}.
Install them in the MSYS2 ${env_name^^} terminal with:
    pacman -S $pkg-cgal $pkg-eigen3 $pkg-boost $pkg-gmp $pkg-mpfr
then run this script again."
                else
                    die "cannot build gudhi_alpha: these libraries are missing for $GXX: ${missing[*]}.
Install CGAL, Eigen 3, Boost, GMP and MPFR for this compiler (Linux: sudo apt install
libcgal-dev libeigen3-dev libboost-dev libgmp-dev libmpfr-dev), or set EIGEN_DIR if
Eigen is installed in an unusual place."
                fi
            fi
            local ver flags=(-O3 -std=c++17 -DCGAL_EIGEN3_ENABLED -DNOMINMAX "${incs[@]}")
            local libs=(-lmpfr -lgmp)
            is_windows && libs+=(-lpsapi)
            ver="$(PATH="$gpath" "$GXX" -dumpversion 2>/dev/null)"
            echo "== Compiling gudhi_alpha with $GXX (log: $log; this takes a minute)"
            : > "$log"
            if [[ -n "$static" ]] && PATH="$gpath" "$GXX" "${flags[@]}" "-DGA_BUILD_FLAGS=\"g++_-O3_static\"" \
                   "$src" -o "$GUDHI_EXE" -static "${libs[@]}" >> "$log" 2>&1; then
                echo "gudhi_alpha.cpp | GNU $ver | -O3_static" > "$info"
            elif PATH="$gpath" "$GXX" "${flags[@]}" "-DGA_BUILD_FLAGS=\"g++_-O3\"" \
                   "$src" -o "$GUDHI_EXE" "${libs[@]}" >> "$log" 2>&1; then
                echo "gudhi_alpha.cpp | GNU $ver | -O3" > "$info"
            else
                die "compiling gudhi_alpha with $GXX failed. First compiler messages:
$(tr -d '\r' < "$log" | grep -iE "error|cannot|not found|undefined|no such" | awk '!seen[$0]++' | head -8 | cut -c1-220)
(all compiler messages: $log)"
            fi
        fi
        if is_windows; then          # DLLs next to the program (needed if not static)
            local gdir dll
            gdir="$(dirname "$GXX")"
            for dll in libgmp-10.dll libmpfr-6.dll libstdc++-6.dll libgcc_s_seh-1.dll libwinpthread-1.dll; do
                [[ -f "$gdir/$dll" ]] && cp -f "$gdir/$dll" "$GUDHI_DIR/" 2>/dev/null
            done
        fi
    elif [[ ! -f "$GUDHI_EXE" ]]; then
        die "gudhi_alpha is not built and no g++ was found to build it (searched the PATH and the
usual MinGW folders). Install g++ (e.g. MSYS2), or set MINGW_BIN=/c/path/to/mingw/bin, or build
it as in gudhi_alpha/README.md, or set GUDHI_EXE=/path/to/gudhi_alpha"
    fi
    echo "   gudhi: $(cat "$info" 2>/dev/null) | $("$GUDHI_EXE" --version 2>/dev/null | cut -d'|' -f2-3)"
    gudhi_startup_test
}
gudhi_startup_test() {           # 3 points: stop here, with the reason, if it cannot run
    local t="$GUDHI_DIR/startup_test" rc missing=""
    printf '0,0\n1,0\n0,1\n' > "$t.csv"
    "$GUDHI_EXE" "$t.csv" "$t.alpha" --r inf --max-dim 2 > /dev/null 2> "$t.err"
    rc=$?
    if [[ $rc -ne 0 ]]; then
        command -v ldd >/dev/null 2>&1 &&
            missing="$(ldd "$GUDHI_EXE" 2>/dev/null | grep -i 'not found' | sed 's/^[[:space:]]*/      /')"
        die "the gudhi_alpha program does not start (start-up test, exit code $rc).
Its error output: $(tr '\n' ' ' < "$t.err" | cut -c1-300)
${missing:+DLLs that Windows cannot find:
$missing
}Program: $GUDHI_EXE
Built with: $(cat "$GUDHI_DIR/build_info.txt" 2>/dev/null)
Please send this message."
    fi
    echo "   gudhi start-up test: ok"
}
cmd_gudhi() {
    CMD=("$GUDHI_EXE" "$1" "$2.alpha" --r "$R" --max-dim "$MAXDIM"
         --precision "$PRECISION" --stats "$2.alpha.json")
}
report_gudhi() { json_report "$1.alpha.json" compute_time_s; }

# --- velocech: Cech H1 barcode code (the repository in velocech/) -----------
# Built from velocech/ by velocech_build/CMakeLists.txt (automatically, and
# rebuilt when the code changes). Its program reads text files whose first line
# is the ambient dimension: the preparation step writes each input in that form
# (same exact numbers, not timed) in $TXT_DIR. It computes the full H1 barcode
# (it ignores r_max and the maximum dimension); bars are on a diameter scale
# (2 x radius). time_s = its own "Time for distances" + "Time after computing
# distances" (reading excluded). OpenMP threads: $VELOCECH_THREADS (empty = all).
# Compiled on the machine it runs on: directly with g++ -O3 -fopenmp when g++
# is available -- found on the PATH or in the usual MinGW folders, or set
# MINGW_BIN=/c/path/to/mingw64/bin (VELOCECH_CXX=default: CMake's default compiler).
VELOCECH_SOURCE="${VELOCECH_SOURCE:-H1_skip.cpp}"
VELOCECH_BUILD_DIR="$HERE/velocech_build"
VELOCECH_EXE_GIVEN="${VELOCECH_EXE:-}"        # set by the user: used as is, never rebuilt
VELOCECH_EXE="${VELOCECH_EXE:-$(exe "$VELOCECH_BUILD_DIR/${VELOCECH_SOURCE%.cpp}")}"
# Windows (Git Bash, MSYS2, Cygwin)?
is_windows() {
    [[ "${OS:-}" == Windows_NT ]] && return 0
    case "$(uname -s 2>/dev/null)" in MINGW*|MSYS*|CYGWIN*) return 0 ;; esac
    return 1
}
# Full path of g++: $MINGW_BIN, then the PATH, then (on Windows) the usual
# MinGW install folders -- so it is found even when its folder is not on the
# PATH of this Git Bash window. Empty if there is none.
find_gxx() {
    local d c
    if [[ -n "${MINGW_BIN:-}" ]]; then
        for c in "$MINGW_BIN/g++.exe" "$MINGW_BIN/g++"; do [[ -f "$c" ]] && { echo "$c"; return; }; done
    fi
    command -v g++ 2>/dev/null && return
    is_windows || return 0
    local cands=() saved
    saved="$(shopt -p nullglob)"; shopt -s nullglob
    cands=(/c/msys64/ucrt64/bin /c/msys64/mingw64/bin /d/msys64/ucrt64/bin /d/msys64/mingw64/bin
           /c/mingw64/bin /c/mingw/bin /c/MinGW/bin /c/winlibs*/mingw64/bin /c/WinLibs*/mingw64/bin
           "/c/Program Files/mingw64/bin" "/c/Program Files/mingw-w64/"*/mingw64/bin
           "/c/Program Files (x86)/mingw-w64/"*/mingw64/bin
           /c/ProgramData/chocolatey/lib/mingw/tools/install/mingw64/bin
           /c/ProgramData/mingw64/mingw64/bin "$HOME"/scoop/apps/mingw*/current/bin
           "$HOME"/AppData/Local/Microsoft/WinGet/Packages/*/mingw64/bin
           /c/TDM-GCC-64/bin /c/Strawberry/c/bin)
    eval "$saved"
    for d in "${cands[@]}"; do
        [[ -f "$d/g++.exe" ]] && { echo "$d/g++.exe"; return; }
    done
}
GXX=""
check_velocech() {
    local src="$HERE/velocech/$VELOCECH_SOURCE" info="$VELOCECH_BUILD_DIR/build_info.txt"
    local log="$VELOCECH_BUILD_DIR/build.log"
    GXX="$([[ "${VELOCECH_CXX:-}" != default ]] && find_gxx)"
    if [[ -n "$VELOCECH_EXE_GIVEN" ]]; then
        [[ -f "$VELOCECH_EXE" ]] || die "VELOCECH_EXE=$VELOCECH_EXE does not exist"
        echo "   velocech: $VELOCECH_EXE (given, not rebuilt)"
        velocech_startup_test; return
    fi
    [[ -f "$src" ]] || die "$src not found: put the velocech repository in $HERE/velocech"
    mkdir -p "$VELOCECH_BUILD_DIR"
    if [[ -n "$GXX" ]]; then
        # GCC: compile the single source file directly (as in velocech's notes:
        # g++ -O3 -fopenmp -std=c++17) on this machine; again whenever the source
        # changed, or an earlier build used another compiler or (on Windows) was
        # not static. On Windows the program is linked STATICALLY: the MinGW
        # runtime (libstdc++, libgomp, ...) goes into the .exe, so it needs no
        # DLL from the PATH; if that fails, a normal build is made and the DLLs
        # are copied next to the program (below).
        local need="" f static=""
        is_windows && static="-static"
        [[ -f "$VELOCECH_EXE" ]] || need=1
        grep -q "GNU" "$info" 2>/dev/null || need=1
        [[ -n "$static" ]] && ! grep -q "static" "$info" 2>/dev/null && need=1
        for f in "$src" "$HERE/velocech/"*.h; do [[ "$f" -nt "$VELOCECH_EXE" ]] && need=1; done
        if [[ -n "$need" ]]; then
            local flags=(-O3 -std=c++17 -I"$HERE/velocech") ver
            [[ "$(uname -s)" == Linux ]] && flags=(-I"$VELOCECH_BUILD_DIR/linux_shim" "${flags[@]}")
            # The compiler's own folder goes first on the PATH while compiling:
            # g++.exe runs cc1plus, as and ld, which need the DLLs of that folder
            # (as in MSYS2's own terminal); without it every compile fails when
            # the folder is not on the PATH of this window.
            local gpath="$(dirname "$GXX"):$PATH"
            ver="$(PATH="$gpath" "$GXX" -dumpversion 2>/dev/null)"
            echo "== Compiling velocech with $GXX ($VELOCECH_SOURCE; log: $log)"
            : > "$log"
            if [[ -n "$static" ]] && PATH="$gpath" "$GXX" "${flags[@]}" -fopenmp "$src" -o "$VELOCECH_EXE" -static >> "$log" 2>&1; then
                echo "$VELOCECH_SOURCE | GNU $ver | -O3_OpenMP_static" > "$info"
            elif PATH="$gpath" "$GXX" "${flags[@]}" -fopenmp "$src" -o "$VELOCECH_EXE" >> "$log" 2>&1; then
                echo "$VELOCECH_SOURCE | GNU $ver | -O3_OpenMP" > "$info"
            elif PATH="$gpath" "$GXX" "${flags[@]}" "$src" -o "$VELOCECH_EXE" $static >> "$log" 2>&1; then
                echo "$VELOCECH_SOURCE | GNU $ver | -O3_no_OpenMP${static:+_static}" > "$info"
                echo "   (OpenMP not available: built single-threaded)"
            else
                die "compiling velocech with $GXX failed. First compiler messages:
$(tr -d '\r' < "$log" | grep -iE "error|cannot|not found|undefined|no such" | awk '!seen[$0]++' | head -8 | cut -c1-220)
(all compiler messages: $log)"
            fi
        fi
    elif command -v cmake >/dev/null 2>&1 && [[ -f "$VELOCECH_BUILD_DIR/CMakeLists.txt" ]] &&
         ! grep -q "GNU" "$info" 2>/dev/null; then
        # no g++ at all: CMake with its default compiler (e.g. Visual C++, no OpenMP)
        echo "== Building velocech with CMake ($VELOCECH_SOURCE; log: $log)"
        { cmake -S "$VELOCECH_BUILD_DIR" -B "$VELOCECH_BUILD_DIR/build" -DCMAKE_BUILD_TYPE=Release \
                -DVELOCECH_SOURCE="$VELOCECH_SOURCE" &&
          cmake --build "$VELOCECH_BUILD_DIR/build" --config Release; } > "$log" 2>&1 ||
            die "building velocech failed; the compiler messages are in $log"
    fi
    [[ -f "$VELOCECH_EXE" ]] || die "velocech program not found ($VELOCECH_EXE): no g++ or cmake to build it.
Install g++ (e.g. MinGW), or set MINGW_BIN=/c/path/to/mingw64/bin, or VELOCECH_EXE=/path/to/the/program"
    # Windows: also put the compiler's runtime DLLs NEXT TO the program. Windows
    # looks in the program's own folder first, so these are used rather than
    # missing or other copies on the PATH (e.g. Anaconda's libstdc++ / libgomp).
    if is_windows && [[ -n "$GXX" ]]; then
        local gdir dll copied=()
        gdir="$(dirname "$GXX")"
        for dll in libgomp-1.dll libstdc++-6.dll libgcc_s_seh-1.dll libgcc_s_dw2-1.dll \
                   libwinpthread-1.dll; do
            if [[ -f "$gdir/$dll" ]]; then
                cp -f "$gdir/$dll" "$VELOCECH_BUILD_DIR/" 2>/dev/null && copied+=("$dll")
            fi
        done
        [[ ${#copied[@]} -gt 0 ]] && echo "   runtime DLLs next to the program: ${copied[*]}"
    fi
    echo "   velocech: $(cat "$info" 2>/dev/null)"
    velocech_startup_test
}
# Start-up test on a 3-point cloud: stop here, with the reason, rather than
# failing on every input.
velocech_startup_test() {
    local info="$VELOCECH_BUILD_DIR/build_info.txt"
    local smoke="$VELOCECH_BUILD_DIR/startup_test" rc
    mkdir -p "$VELOCECH_BUILD_DIR"
    printf '2\n0 0\n1 0\n0 1\n' > "$smoke.txt"
    "$VELOCECH_EXE" "$smoke.txt" "$smoke.bars" > /dev/null 2> "$smoke.err"
    rc=$?
    if [[ $rc -ne 0 ]]; then
        local missing=""
        command -v ldd >/dev/null 2>&1 &&
            missing="$(ldd "$VELOCECH_EXE" 2>/dev/null | grep -i 'not found' | sed 's/^[[:space:]]*/      /')"
        die "the velocech program does not start (start-up test, exit code $rc).
Its error output: $(tr '\n' ' ' < "$smoke.err" | cut -c1-300)
${missing:+DLLs that Windows cannot find:
$missing
}Program: $VELOCECH_EXE
Built with: $(cat "$info" 2>/dev/null)
g++: ${GXX:-NOT FOUND (searched the PATH and the usual MinGW folders): set MINGW_BIN=/c/path/to/your/mingw/bin}
Please send this message."
    fi
    echo "   start-up test: ok"
}
cmd_velocech() {
    local in="$TXT_DIR/$(basename "$1").txt"
    if [[ -n "${VELOCECH_THREADS:-}" ]]; then
        CMD=(env "OMP_NUM_THREADS=$VELOCECH_THREADS" "$VELOCECH_EXE" "$in" "$2.bars")
    else
        CMD=("$VELOCECH_EXE" "$in" "$2.bars")
    fi
}
report_velocech() {
    [[ -f "$1.stderr" ]] && awk '{ gsub(/\r/, "") }
        /^Time for distances:/ { d = $4 }
        /^Time after computing distances:/ { h = $5 }
        /^Cofaces skipped:/ { gsub(/,/, "", $3); s = $3; r = $5 }
        END { if (d != "" && h != "") { printf "time_s=%.9f\ndistance_time_s=%s\nh1_time_s=%s\n", d + h, d, h }
              if (s != "") { print "cofaces_skipped=" s; print "cofaces_reduced=" r } }' "$1.stderr"
    [[ -f "$1.bars" ]] && echo "n_h1_bars=$(grep -c . "$1.bars")"
}

# --- TEMPLATE for another code (copy, adapt, and add its name to ALL_ALGORITHMS)
# Example: a program that reads a CSV point cloud and prints "time: <seconds>".
# MYCODE_EXE="/path/to/mycode"
# check_mycode()  { [[ -f "$MYCODE_EXE" ]] || die "mycode not found: $MYCODE_EXE"; }
# cmd_mycode()    { CMD=("$MYCODE_EXE" "$1" --radius "$R" --max-dim "$MAXDIM" --output "$2.out"); }
# report_mycode() { awk '/^time:/ {print "time_s=" $2}' "$1.stdout"; }

# ===========================================================================
# Settings
# ===========================================================================
INPUTS="" TLIMIT="" R="" MAXDIM="" SOLVER="" PRECISION="" ALGOS="" OUT="" RESUME="" SKIP_LARGER=""
VELOCECH_THREADS="${VELOCECH_THREADS:-}" ASKED_THREADS=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --inputs) INPUTS="$2"; shift 2 ;;
        --time-limit) TLIMIT="$2"; shift 2 ;;
        --r) R="$2"; shift 2 ;;
        --max-dim) MAXDIM="$2"; shift 2 ;;
        --solver) SOLVER="$2"; shift 2 ;;
        --precision) PRECISION="$2"; shift 2 ;;
        --algorithms) ALGOS="$2"; shift 2 ;;
        --out) OUT="$2"; shift 2 ;;
        --resume) RESUME=1; shift ;;
        --skip-larger) SKIP_LARGER="$2"; shift 2 ;;
        --velocech-threads) VELOCECH_THREADS="$2"; ASKED_THREADS=1; shift 2 ;;
        -h|--help) sed -n '2,48p' "$0"; exit 0 ;;
        *) die "unknown option $1 (see --help)" ;;
    esac
done

if [[ -n "$RESUME" ]]; then      # take the settings of the interrupted experiment
    [[ -n "$OUT" && -f "$OUT/settings.json" ]] || die "--resume needs --out DIR of an earlier experiment"
    eval "$("$PY" - "$OUT/settings.json" <<'PYEOF'
import json, shlex, sys
s = json.load(open(sys.argv[1]))
for k, v in (("INPUTS", "inputs"), ("TLIMIT", "time_limit_s"), ("R", "r_max"), ("MAXDIM", "max_dim"),
             ("SOLVER", "solver"), ("PRECISION", "gudhi_precision"), ("ALGOS", "algorithms"),
             ("SKIP_LARGER", "skip_larger"), ("VELOCECH_THREADS", "velocech_threads")):
    print(f"{k}={shlex.quote(str(s.get(v, 'failure' if k == 'SKIP_LARGER' else '')))}")
PYEOF
)"
fi

ask() {    # ask VARIABLE "question" default
    local ans=""
    # The question is printed explicitly: "read -p" shows it only when bash sees
    # a terminal, which some Git Bash setups do not report -- the question would
    # then be invisible while the script waits for the answer.
    printf '%s' "$2${3:+ [$3]}: " >&2
    IFS= read -r ans || true
    ans="${ans%$'\r'}"
    printf -v "$1" '%s' "${ans:-$3}"
}
STEP="asking the settings"
[[ -n "$INPUTS" ]]    || ask INPUTS "Folder containing the input subfolders" ""
[[ -n "$TLIMIT" ]]    || ask TLIMIT "Time limit per run (seconds)" "600"
[[ -n "$R" ]]         || ask R "r_max for alphadual and gudhi (a number, or inf)" "inf"
[[ -n "$MAXDIM" ]]    || ask MAXDIM "Maximum simplex dimension" "2"
[[ -n "$SOLVER" ]]    || ask SOLVER "Solver for alphadual / half_diameter (daqp or native)" "daqp"
[[ -n "$PRECISION" ]] || ask PRECISION "GUDHI precision (fast, safe or exact)" "safe"
[[ -n "$ALGOS" ]]     || ask ALGOS "Algorithms to run" "$(IFS=,; echo "${ALL_ALGORITHMS[*]}")"
if [[ ",$ALGOS," == *",velocech,"* && -z "$ASKED_THREADS" && -z "$RESUME" && -z "$VELOCECH_THREADS" ]]; then
    echo "velocech's distance step can use several cores (OpenMP); the other codes use one."
    ask VELOCECH_THREADS "Threads for velocech (Enter = all cores, 1 = single-threaded like the others)" ""
fi
[[ -z "$VELOCECH_THREADS" || "$VELOCECH_THREADS" =~ ^[1-9][0-9]*$ ]] || die "velocech threads must be a positive integer: $VELOCECH_THREADS"
OUT="${OUT:-experiments_$(date +%Y%m%d_%H%M%S)}"
SKIP_LARGER="${SKIP_LARGER:-failure}"

# Accept Windows paths (C:\...\inputs), with or without surrounding quotes.
win_path() {
    local p="$1"
    p="${p#\"}"; p="${p%\"}"; p="${p#\'}"; p="${p%\'}"
    if [[ "$p" =~ ^[A-Za-z]: || "$p" == *\\* ]] && command -v cygpath >/dev/null 2>&1; then
        p="$(cygpath -u "$p")"
    fi
    printf '%s' "$p"
}
INPUTS="$(win_path "$INPUTS")"
OUT="$(win_path "$OUT")"
[[ -d "$INPUTS" ]] || die "input folder not found: $INPUTS"
[[ "$TLIMIT" =~ ^[0-9]+([.][0-9]+)?$ ]] || die "time limit must be a number of seconds: $TLIMIT"
[[ "$R" =~ ^(inf|[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?)$ ]] || die "r_max must be a number or inf: $R"
[[ "$MAXDIM" =~ ^[0-9]+$ ]] || die "maximum dimension must be an integer: $MAXDIM"
[[ "$SOLVER" == daqp || "$SOLVER" == native ]] || die "solver must be daqp or native"
[[ "$PRECISION" =~ ^(fast|safe|exact)$ ]] || die "precision must be fast, safe or exact"
[[ "$SKIP_LARGER" =~ ^(failure|timeout|never)$ ]] || die "--skip-larger must be failure, timeout or never"
STEP="checking the settings and the programs"
IFS=',' read -r -a SELECTED <<< "$ALGOS"
for a in "${SELECTED[@]}"; do
    declare -F "cmd_$a" >/dev/null || die "unknown algorithm '$a' (known: ${ALL_ALGORITHMS[*]})"
    STEP="checking / building the program for '$a'"
    if declare -F "check_$a" >/dev/null; then "check_$a"; fi
done
if [[ "$R" == inf && "$SOLVER" == native ]] && [[ ",$ALGOS," == *",alphadual,"* ]]; then
    die "r_max = inf needs the daqp solver for alphadual"
fi

mkdir -p "$OUT" || die "cannot create $OUT"
OUT="$(cd "$OUT" && pwd)"
INPUTS="$(cd "$INPUTS" && pwd)"
if [[ -z "$RESUME" ]]; then
    "$PY" - "$OUT/settings.json" "$INPUTS" "$TLIMIT" "$R" "$MAXDIM" "$SOLVER" "$PRECISION" "$ALGOS" \
        "$("$ALPHADUAL_EXE" --version 2>/dev/null)" "$("$GUDHI_EXE" --version 2>/dev/null)" \
        "$SKIP_LARGER" "${VELOCECH_THREADS:-all}" "$VELOCECH_SOURCE" \
        "$(cat "$VELOCECH_BUILD_DIR/build_info.txt" 2>/dev/null)" \
        "$(git -C "$HERE/velocech" log -1 --format='%h %ad %s' --date=short 2>/dev/null)" <<'PYEOF'
import json, sys, platform, datetime
a = sys.argv[1:]
json.dump(dict(inputs=a[1], time_limit_s=float(a[2]), r_max=a[3], max_dim=int(a[4]), solver=a[5],
               gudhi_precision=a[6], algorithms=a[7], alphadual_version=a[8], gudhi_version=a[9],
               skip_larger=a[10], velocech_threads=a[11], velocech_source=a[12],
               velocech_build=a[13].strip(), velocech_commit=a[14],
               python=sys.version.split()[0], platform=platform.platform(),
               started=datetime.datetime.now().isoformat(timespec="seconds")),
          open(a[0], "w"), indent=2)
PYEOF
fi

# ===========================================================================
# Runs
# ===========================================================================
# From here on, everything printed is also written to $OUT/run_all.log.
exec 3>&1 4>&2
exec > >(tee -a "$OUT/run_all.log") 2>&1
TEE_PID=$!
echo "=== $(date '+%Y-%m-%d %H:%M:%S')"
echo "Inputs: $INPUTS"
echo "Results: $OUT"
echo "Algorithms: ${SELECTED[*]} | time limit ${TLIMIT} s | r_max $R | max dim $MAXDIM | skip larger after: $SKIP_LARGER"
shopt -s nullglob
subs=("$INPUTS"/*/)
[[ ${#subs[@]} -gt 0 ]] || die "no subfolders in $INPUTS (the inputs must be in subfolders of it)"
NRUNS=0
# FAILED_AT["algorithm|group"] = smallest number of points at which the
# algorithm failed on that group (timeout, or error unless --skip-larger
# timeout); larger inputs of the group are then skipped -- in every subfolder.
declare -A FAILED_AT=()
for subdir in "${subs[@]}"; do
    sub="$(basename "$subdir")"
    files=()
    for f in "$subdir"*; do
        case "${f,,}" in
            *.csv|*.txt|*.xyz|*.bin) [[ "$(basename "$f")" == index.csv ]] || files+=("$f") ;;
        esac
    done
    if [[ ${#files[@]} -eq 0 && -d "${subdir}points" ]]; then
        # a benchmark result folder: its point clouds are in points/
        for f in "${subdir}points/"*.csv; do
            [[ "$(basename "$f")" == index.csv ]] || files+=("$f")
        done
    fi
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "== $sub: no input files (*.csv *.txt *.xyz *.bin) in it, skipped"; continue
    fi
    RES="$OUT/$sub"
    BIN_DIR="$RES/bin"
    TXT_DIR=""                                   # velocech's text inputs (only if needed)
    [[ ",$ALGOS," == *",velocech,"* ]] && TXT_DIR="$RES/velocech_txt" && mkdir -p "$TXT_DIR"
    mkdir -p "$BIN_DIR"
    echo "== $sub: ${#files[@]} input(s)"
    STEP="preparing the inputs of subfolder '$sub'"
    # Preparation (never timed): for each input its number of points, ambient
    # dimension and GROUP (shape + shape parameters such as noise + ambient
    # dimension: inputs of one group differ only in size and run); the inputs
    # converted to alphadual's .bin format; and the run order (order.tsv):
    # group by group, in increasing number of points.
    "$PY" - "$HERE" "$RES" "$BIN_DIR" "$TXT_DIR" "$sub" "${files[@]}" <<'PYEOF' || echo "   (warning: some inputs could not be prepared)"
import csv, json, os, re, sys
here, res, bindir, txtdir, sub, files = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6:]
import numpy as np
sys.path.insert(0, here)
from alpha_at_half_diameter import read_points
import alphadual_io as io
index = {}          # file name -> row of an index.csv next to the inputs
for d in sorted({os.path.dirname(f) for f in files}):
    p = os.path.join(d, "index.csv")
    if os.path.exists(p):
        for row in csv.DictReader(open(p)):
            index[row.get("file", "")] = row
info, order = {}, []
for f in files:
    name = os.path.basename(f)
    stem = os.path.splitext(name)[0]
    m = re.match(r"^(.*)_(\d+)_(\d+)_(\d+)$", stem)       # <shape>_<points>_<dim>_<run>
    run = int(m.group(4)) if m else 0
    if name in index:       # shape and parameters (noise, ...) known exactly
        shape = f"{index[name].get('shape', '')} {index[name].get('shape_params', '')}".strip()
    elif m:                 # from the name only (noise unknown): within this subfolder
        shape = f"{m.group(1)} (subfolder {sub})"
    else:
        shape = f"{stem} (subfolder {sub})"
    try:
        X = read_points(f)
        n, dim = int(X.shape[0]), int(X.shape[1])
        if txtdir:          # velocech's text format: first line = ambient dimension
            t = os.path.join(txtdir, name + ".txt")
            if not os.path.exists(t):
                with open(t, "w", newline="\n") as fh:
                    fh.write(f"{dim}\n")
                    np.savetxt(fh, X, fmt="%.17g")
        b = os.path.join(bindir, name + ".bin")
        if not os.path.exists(b):
            io.write_input(b, X)
        group = f"{shape} | ambient dim {dim}"
        info[name] = dict(n_points=n, ambient_dim=dim, group=group)
    except Exception as e:
        n, group = 0, f"unreadable: {name}"
        info[name] = dict(error=str(e), group=group)
        print(f"   cannot read {name}: {e}")
    order.append((group, n, run, name, f))
json.dump(info, open(os.path.join(res, "inputs.json"), "w"), indent=1)
first = {}                  # groups in order of first appearance
for g, *_ in order:
    first.setdefault(g, len(first))
order.sort(key=lambda t: (first[t[0]], t[1], t[2], t[3]))
with open(os.path.join(res, "order.tsv"), "w", newline="\n") as fh:
    for g, n, run, name, f in order:
        fh.write(f"{f}\t{name}\t{g}\t{n}\n")
PYEOF
    # inputs in that order (path, name, group, number of points)
    files=(); names=(); groups=(); sizes=()
    while IFS=$'\t' read -r f nm g n; do
        files+=("$f"); names+=("$nm"); groups+=("$g"); sizes+=("$n")
    done < "$RES/order.tsv"
    touch "$RES/done.txt"
    for algo in "${SELECTED[@]}"; do
        mkdir -p "$RES/outputs/$algo"
        # On --resume: add the failures already recorded in this subfolder.
        if [[ "$SKIP_LARGER" != never && -f "$RES/runs.jsonl" ]]; then
            while IFS=$'\t' read -r g n; do
                fk="$algo|$g"
                if [[ -z "${FAILED_AT[$fk]:-}" || "$n" -lt "${FAILED_AT[$fk]}" ]]; then FAILED_AT[$fk]="$n"; fi
            done < <("$PY" - "$RES/runs.jsonl" "$RES/order.tsv" "$algo" "$SKIP_LARGER" <<'PYEOF'
import json, sys
runs, order, algo, mode = sys.argv[1:5]
meta = {}
for line in open(order):
    p = line.rstrip("\n").split("\t")
    if len(p) == 4:
        meta[p[1]] = (p[2], p[3])
bad = {"timeout"} if mode == "timeout" else {"timeout", "error"}
for line in open(runs):
    try:
        r = json.loads(line)
    except json.JSONDecodeError:
        continue
    if r.get("algorithm") == algo and r.get("status") in bad and r.get("input") in meta:
        print(*meta[r["input"]], sep="\t")
PYEOF
)
        fi
        for k in "${!files[@]}"; do
            f="${files[$k]}"; name="${names[$k]}"; group="${groups[$k]}"; npts="${sizes[$k]}"
            id="$algo|$name"
            if [[ -n "$RESUME" ]] && grep -Fxq "$id" "$RES/done.txt"; then
                echo "   [$algo] $name: already done"; continue
            fi
            fk="$algo|$group"
            if [[ -n "${FAILED_AT[$fk]:-}" && "$npts" -gt "${FAILED_AT[$fk]}" ]]; then
                msg="skipped: failed at ${FAILED_AT[$fk]} points on the same input type, ambient dimension and parameters"
                printf '{"id": "%s", "algorithm": "%s", "input": "%s", "status": "skipped", "message": "%s"}\n' \
                    "$id" "$algo" "$name" "$msg" >> "$RES/runs.jsonl"
                echo "$id" >> "$RES/done.txt"
                echo "   [$algo] $name: $msg"
                "$PY" "$HERE/experiments_to_excel.py" "$RES" > /dev/null || true
                continue
            fi
            prefix="$RES/outputs/$algo/$name"
            STEP="running '$algo' on '$sub/$name'"
            CMD=()
            "cmd_$algo" "$f" "$prefix"
            "$PY" "$HERE/measure_run.py" --record "$RES/runs.jsonl" --timeout "$TLIMIT" \
                --stdout "$prefix.stdout" --stderr "$prefix.stderr" \
                --set "id=$id" --set "algorithm=$algo" --set "input=$name" -- "${CMD[@]}"
            if declare -F "report_$algo" >/dev/null; then
                "report_$algo" "$prefix" 2>/dev/null | sed "s/^/$id	/" >> "$RES/reported.tsv"
            fi
            echo "$id" >> "$RES/done.txt"
            NRUNS=$((NRUNS + 1))
            last="$(tail -n 1 "$RES/runs.jsonl" | grep -o '"status": "[a-z]*"' | sed 's/.*: "//; s/"//')"
            if [[ "$SKIP_LARGER" != never ]] &&
               [[ "$last" == timeout || ( "$last" == error && "$SKIP_LARGER" == failure ) ]]; then
                if [[ -z "${FAILED_AT[$fk]:-}" || "$npts" -lt "${FAILED_AT[$fk]}" ]]; then
                    FAILED_AT[$fk]="$npts"
                fi
            fi
            # progress line from the record just written
            "$PY" - "$RES/runs.jsonl" "$RES/reported.tsv" "$id" "$prefix.stderr" "$HERE" <<'PYEOF'
import json, os, sys
r = json.loads(open(sys.argv[1]).read().strip().splitlines()[-1])
t, status = "", r["status"]
if os.path.exists(sys.argv[2]):
    for line in open(sys.argv[2]):
        parts = line.rstrip("\n").split("\t", 1)
        if len(parts) == 2 and parts[0] == sys.argv[3]:
            if parts[1].startswith("time_s="):
                t = parts[1][7:]
            elif parts[1].startswith("status=") and status == "ok":
                status = "error (reported by the code)"
mem = r.get("peak_memory_mb")
why = ""
if r["status"] == "error":       # say why: the program's own message, or the exit code
    sys.path.insert(0, sys.argv[5])
    from experiments_to_excel import exit_code_message
    try:
        lines = [l.strip() for l in open(sys.argv[4], errors="replace") if l.strip()]
    except OSError:
        lines = []
    why = (lines[-1][:150] if lines else "") or exit_code_message(r.get("exit_code")) \
        or f"exit code {r.get('exit_code')}"
print(f"   [{r['algorithm']}] {r['input']}: {status}"
      + (f", time_s {float(t):.4f}" if t else "")
      + f", wall {r['wall_time_s']:.3f} s"
      + (f", peak memory {mem:.0f} MB" if mem is not None else "")
      + (f"\n      reason: {why}" if why else ""))
PYEOF
            "$PY" "$HERE/experiments_to_excel.py" "$RES" > /dev/null ||
                echo "   (warning: summary.xlsx not updated this time; results are kept)"
        done
    done
    echo "   summary: $RES/summary.xlsx"
done
echo
if [[ $NRUNS -eq 0 ]]; then
    echo "NOTHING WAS RUN: no input files were found (or all runs were already done)."
    echo "The input folder must contain SUBFOLDERS holding the point clouds (*.csv *.txt *.xyz *.bin)."
else
    echo "Done: $NRUNS run(s). Results in $OUT (one summary.xlsx per subfolder; log: run_all.log)"
fi
