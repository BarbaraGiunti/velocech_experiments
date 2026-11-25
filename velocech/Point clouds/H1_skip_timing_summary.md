# How long H1_skip takes on these point clouds (rough guide)

Written 2026-11-13. These are quick, strategic measurements (about 30 minutes of runs) plus educated extrapolation. They are meant for planning a proper benchmark, not as final numbers.

## Setup

- Code: `H1_skip.cpp` with the EXH change (fully reduced stored columns), built with `g++ -O3 -fopenmp -std=c++17`.
- Machine: laptop, Intel i7-1355U, 16 GB RAM. Every run was pinned to the 4 fast cores (P-cores). Only the distance step uses several threads; the main computation runs on one core.
- Input: the `*_run1.bin` files (plus two `run2` files) read directly. Note: `H1_skip.cpp` itself only reads text files at the moment. To read the `.bin` files, change `load_ascii` to `load_point_cloud` in `main` (that function is already in `inputs.h` and detects the format automatically).
- "Distance step" = computing and sorting all pairwise distances. "Main step" = everything after that. Peak memory is the most RAM the process used at any time.

## Measurements

| Cloud | Points | Ambient dim | Distance step | Main step | Peak memory | H1 bars |
|---|---|---|---|---|---|---|
| 2-sphere | 2000 | 3 | 0.1 s | 0.2 s | 87 MB | 3,782 |
| 2-sphere | 5000 | 3 | 0.8-0.9 s | 1.6-2.7 s (same file twice; run2 file: 1.5 s) | 514 MB | 10,773 |
| 2-sphere | 10000 | 3 | 4.1 s | 15.9 s | 2.04 GB | 23,556 |
| 2-sphere | 10000 | 20 | 4.5 s | 19.0 s | 2.04 GB | 23,556 |
| 2-sphere | 2000 | 100 | 0.2 s | 0.9 s | 90 MB | 3,782 |
| 2-sphere | 5000 | 100 | 1.0 s | 7.8 s | 522 MB | 10,773 |
| torus | 2000 | 3 | 0.1 s | 0.2 s | 87 MB | 2,768 |
| torus | 5000 | 3 | 0.8 s | 2.2 s | 514 MB | 7,814 |
| torus | 10000 | 3 | 4.0 s | 14.0 s | 2.04 GB | 17,880 |
| torus | 5000 | 20 | 0.9 s | 4.4 s | 516 MB | 7,814 |
| torus | 5000 | 100 | 1.0 s | 8.8 s | 522 MB | 7,814 |
| torus | 10000 | 100 | 4.7 s | 47.6 s | 2.04 GB | 17,880 |
| 19-sphere | 1000 | 20 | 0.05 s | 5.6 s | 48 MB | 324,650 |
| 19-sphere | 1500 | 20 | 0.08 s | 16.4 s (run2: 15.7 s) | 95 MB | 658,441 |
| 19-sphere | 2000 | 20 | 0.14 s | 36.9 s | 161 MB | 1,092,838 |
| 19-sphere | 2000 | 50 | 0.12 s | 63.0 s | 161 MB | 1,092,838 |
| 19-sphere | 2000 | 100 | 0.16 s | 121.4 s | 164 MB | 1,092,838 |
| 19-sphere | 3000 | 20 | 0.29 s | 113 s | 356 MB | 2,187,036 |
| 19-sphere | 4000 | 20 | 0.60 s | 312 s | 609 MB | 3,512,513 |
| 19-sphere | 5000 | 20 | 0.79 s | 576 s (9.6 min) | 949 MB | 5,095,757 |

## What the numbers say

**1. Two very different kinds of input.** The 2-sphere and the torus (low intrinsic dimension) are easy: a few thousand bars, and even 10,000 points take well under a minute. The 19-sphere is hard: it has an enormous number of short-lived bars (about 1.1 million at 2000 points, 5.1 million at 5000), and the time grows much faster with the number of points.

**2. How time grows with the number of points.**

- 2-sphere and torus: doubling the points makes the main step roughly 6 times slower (time grows like n^2.5 to n^2.7). The distance step grows a bit more slowly (about 4.5 times per doubling).
- 19-sphere: doubling the points makes the main step roughly 7 to 8 times slower (from 1000 to 5000 points the time grows like n^2.9; the step from 3000 to 4000 looked steeper, but 4000 to 5000 was back on the n^2.9 curve, so that was probably noise). The number of bars grows like n^1.6 to n^1.7.

**3. The ambient dimension matters, even though the shape is the same.** The clouds in higher ambient dimension are the same shapes embedded differently (the bar counts are identical across ambient dimensions, so this is a useful sanity check in the real test: the bars should agree up to rounding). But every distance computation costs more in more coordinates, and the main step computes many distances:

- 2-sphere / torus at 5000-10000 points: R^20 is 1.2 to 2 times slower than R^3, R^100 is 3 to 3.5 times slower.
- 19-sphere at 2000 points: R^50 is 1.7 times slower than R^20, R^100 is 3.3 times slower.

**4. Memory is not a problem for these files.**

- 2-sphere / torus: about 20 bytes x n^2, i.e. 2.04 GB at 10,000 points (it would be about 8 GB at 20,000).
- 19-sphere: about 38-40 bytes x n^2 (extra columns for all the bars), but n is smaller in practice: 0.6 GB at 4000 points, 0.95 GB at 5000 points, so roughly 4 GB at 10,000 points.
- The ambient dimension barely changes memory.

**5. Run-to-run variation is large on this laptop.** The very same file (2-sphere, 5000 points, run1) took 2.7 s the first time and 1.6 s the second time; the run2 file took 1.5 s. So the variation comes mostly from the laptop (background load, power management), not from the random sample: different samples of the same size behave very similarly (the 19-sphere at 1500 points: 16.4 s for run1, 15.7 s for run2). Expect +-20% (sometimes more) for single runs on a laptop.

## Estimates for the full grid (on this laptop, one job at a time)

Main step + distance step, per single run:

| Cloud | n = 2000 | n = 5000 | n = 10000 |
|---|---|---|---|
| 2-sphere / torus, R^3 | < 1 s | 3-4 s | 20 s |
| 2-sphere / torus, R^20 | 1 s | 5 s | 25-35 s |
| 2-sphere / torus, R^100 | 1-2 s | 9-10 s | 50-60 s |
| 19-sphere, R^20 | 40 s | 10 min (measured) | 1-1.5 hours |
| 19-sphere, R^50 | 65 s | 15-20 min | 2-2.5 hours |
| 19-sphere, R^100 | 2 min | 30-35 min | 4-5 hours |

The 19-sphere estimates assume the main step keeps growing like n^2.9, and that the slowdown from the ambient dimension (x1.7 for R^50, x3.3 for R^100, measured at 2000 points) stays the same for larger n. Both are guesses; the 10,000-point numbers could easily be off by a factor of 1.5 in either direction.

Totals (all 13 sizes from 500 to 10,000 points):

- 2-sphere or torus: about 10-12 minutes per run for all 7 ambient dimensions together, so about 2 hours for both shapes and all 5 runs.
- 19-sphere: about 3.7 hours per run in R^20, 6 hours in R^50 and 12 hours in R^100, so about 22 hours per run and about 4-5 days for all 5 runs, one job at a time. Most of it is the largest sizes: the sizes 8000-10,000 alone are about 70% of the total.

## Suggestions for the proper test

- **Run several jobs in parallel**, one per core. The main step uses only one core, so a machine with many cores finishes the grid many times faster. Memory allows this easily for the 2-sphere/torus (at most 2 GB per job) and for the 19-sphere (at most roughly 4 GB per job at 10,000 points). Set `OMP_NUM_THREADS` to 1 or 2 when running many jobs at once, so the distance steps don't fight over cores.
- **Put a time limit on the 19-sphere jobs**, or cap n (for example at 5000-6000 points), unless you have days of machine time. A faster desktop CPU might be 1.5-2 times faster per core than this laptop, not more.
- **Disk space for the bar files:** about 40 bytes per bar. The 19-sphere files get large: 42 MB at 2000 points, 198 MB at 5000 points, and probably 0.6-0.8 GB per file at 10,000 points (about 16 million bars). All 19-sphere bar files together would be roughly 30 GB. The 2-sphere/torus bar files are under 1 MB. If you only need timings, consider not writing the 19-sphere bars, or writing only bars above some length.
- **Pin jobs to fast cores or use a machine without efficiency cores**, and repeat a few runs to see the noise.
- Use the same compiler flags (`-O3`; avoid `-march=native`, which changes the last digits of the results).
