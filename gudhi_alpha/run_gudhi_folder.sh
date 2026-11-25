#!/usr/bin/env bash
# ===========================================================================
# run_gudhi_folder.sh -- run gudhi_alpha on every point cloud of a folder.
#
#   ./run_gudhi_folder.sh INPUT_FOLDER [--out RESULTS_FOLDER] [gudhi_alpha options]
#
# Example (at most 10 minutes per input):
#   ./run_gudhi_folder.sh ../myrun/points --r 0.4 --max-dim 3 --time-limit 600
#
# gudhi_alpha options: --r R|inf|half-diameter (required), --max-dim D,
#   --precision fast|safe|exact.
# --time-limit SECONDS: each run is started under "timeout", a separate
#   process that only sleeps until the limit, then stops the run; nothing
#   extra runs inside gudhi_alpha, so its timing is not affected.
# Inputs: *.csv *.txt *.xyz *.bin in INPUT_FOLDER (index.csv is skipped), in
# name order. Each input runs in its OWN process, so an input stopped by a
# limit (or crashing) does not affect the others.
#
# Results (RESULTS_FOLDER, default gudhi_results_<date>_<time>):
#   results.csv       one line per input: status (ok / time_limit / error),
#                     times, peak memory, simplex counts
#   <input>.alpha     ordered filtration (alphadual output format)
#   <input>.alpha.json   statistics of that run
# ===========================================================================
set -uo pipefail     # no "set -e": keep going when one input fails

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$HERE/gudhi_alpha"
[[ -f "$EXE.exe" ]] && EXE="$EXE.exe"
if [[ ! -f "$EXE" ]]; then
    echo "ERROR: $EXE not found: build it first (see the top of CMakeLists.txt)" >&2
    exit 1
fi
if [[ $# -lt 1 || "$1" == -h || "$1" == --help ]]; then
    sed -n '2,22p' "$0"; exit 0
fi

IN="$1"; shift
[[ -d "$IN" ]] || { echo "ERROR: $IN is not a folder" >&2; exit 1; }
OUT="gudhi_results_$(date +%Y%m%d_%H%M%S)"
ARGS=()
MAXDIM=2
TLIMIT=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --out) OUT="$2"; shift 2 ;;
        --time-limit) TLIMIT="$2"; shift 2 ;;
        --max-dim) MAXDIM="$2"; ARGS+=("$1" "$2"); shift 2 ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
mkdir -p "$OUT"
RESULTS="$OUT/results.csv"
# Header (gudhi_alpha appends one line per successful run; this script adds
# the lines of runs that were stopped or failed).
hdr="file,status,n_points,ambient_dim,r,max_dim,precision,read_time_s,compute_time_s,alpha_time_s,diameter_time_s,peak_memory_mb"
for ((d = 0; d <= MAXDIM; d++)); do hdr+=",n_dim$d"; done
echo "$hdr" > "$RESULTS"
if [[ -n "$TLIMIT" ]] && ! command -v timeout >/dev/null 2>&1; then
    echo "ERROR: --time-limit needs the 'timeout' command (coreutils)" >&2; exit 1
fi

shopt -s nullglob nocaseglob
files=("$IN"/*.csv "$IN"/*.txt "$IN"/*.xyz "$IN"/*.bin)
shopt -u nocaseglob
inputs=()
for f in "${files[@]}"; do
    [[ "$(basename "$f")" == index.csv ]] || inputs+=("$f")
done
IFS=$'\n' inputs=($(printf '%s\n' "${inputs[@]}" | sort)); unset IFS
[[ ${#inputs[@]} -gt 0 ]] || { echo "ERROR: no input files in $IN" >&2; exit 1; }

echo "$("$EXE" --version)"
echo "${#inputs[@]} input(s); results in $OUT"
k=0
for f in "${inputs[@]}"; do
    k=$((k + 1))
    name="$(basename "$f")"
    if [[ -n "$TLIMIT" ]]; then
        # TERM at the limit, KILL 5 s later if it is still running
        timeout --kill-after=5 "$TLIMIT" "$EXE" "$f" "$OUT/$name.alpha" "${ARGS[@]}" \
            --results-csv "$RESULTS" 2> "$OUT/$name.stderr"
    else
        "$EXE" "$f" "$OUT/$name.alpha" "${ARGS[@]}" --results-csv "$RESULTS" 2> "$OUT/$name.stderr"
    fi
    rc=$?
    case $rc in
        0)   t=$(grep -o '"compute_time_s": [0-9.eE+-]*' "$OUT/$name.alpha.json" | awk '{print $2}')
             echo "[$k/${#inputs[@]}] $name: ok, compute ${t} s" ;;
        124|137)
             if [[ -n "$TLIMIT" ]]; then
                 echo "$name,time_limit (${TLIMIT} s)" >> "$RESULTS"
                 echo "[$k/${#inputs[@]}] $name: stopped by the time limit (${TLIMIT} s)"
                 rm -f "$OUT/$name.alpha" "$OUT/$name.alpha.json"
             else
                 echo "$name,error (exit $rc: killed)" >> "$RESULTS"
                 echo "[$k/${#inputs[@]}] $name: killed (exit $rc)"
             fi ;;
        *)   # crashed: record it
             msg="$(tr '\n,' '  ' < "$OUT/$name.stderr" | cut -c1-200)"
             echo "$name,error (exit $rc): $msg" >> "$RESULTS"
             echo "[$k/${#inputs[@]}] $name: ERROR (exit $rc) $msg" ;;
    esac
    [[ -s "$OUT/$name.stderr" ]] || rm -f "$OUT/$name.stderr"
done
echo "Done: $RESULTS"
