// ===========================================================================
// alphadual.cpp
//
// Alpha filtration up to a chosen dimension D (default 2: vertices, edges,
// triangles) with the
// dual active-set QP method of
//
//     E. Carlsson and J. Carlsson, "Computing the alpha complex using dual
//     active set quadratic programming", arXiv:2310.00536 (2023), Algorithm 1.
//
// C++ port of alpha_dualqp.py (same algorithm, tolerances and outcome codes;
// see that file's READING GUIDE for the mathematics). Two dual QP solvers:
//   --solver daqp   : DAQP (Arnstrom, Bemporad & Axehill 2022), ref. [2] of
//                     the paper, called through its C API (vendored v0.10.3)
//   --solver native : C++ port of our own dual active-set solver
//
// Usage:
//   alphadual INPUT.bin OUTPUT.bin --r R [--max-dim D] [--solver daqp|native]
//             [--stats S.json]           (D = largest simplex dimension, default 2)
//   alphadual INPUT.bin OUTPUT.bin --a1 A1 ...      (weighted: power level a1)
//   R may be "inf" (complete filtration; requires --solver daqp).
//
// Input file (little-endian, written by alphadual_io.py):
//   char[8] "ALPHAIN1" | int64 N | int64 m | int64 has_weights |
//   double X[N*m] (row-major) | double w[N] (only if has_weights)
// Output file (format 2):
//   char[8] "ALPHOUT2" | int64 M | int64 width (= D+1) | int64 counts[width] |
//   int32 simplices[M*width] (vertex ids, -1 padded) | int8 dims[M] |
//   double values[M] (squared entry radius)
//   -- the filtration in entry order: by value, then dimension (faces before
//      cofaces), then vertex ids.
// Stats file (JSON): timings, sizes and QP statistics.
//
// Timing: ONE clock interval, "alpha_time_s", from just before the input file
// is opened to just after the ordered filtration is built (reading included,
// writing excluded). One extra clock reading at the end of the read gives
// "read_time_s", so you can check that reading is a small part; nothing is
// timed inside the computation itself.
// ===========================================================================

// Windows: keep <windows.h> (if any header pulls it in) from defining the
// min/max macros, which would break std::min / std::max.
#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "api.h"          // DAQP (its headers have their own extern "C" guards)
#include "utils.h"        // DAQP: daqp_update_ldp
#include "nanoflann.hpp"  // KD-tree for the Cech graph

using std::vector;
using i64 = int64_t;
using i32 = int32_t;

// QP outcome codes (same values as alpha_dualqp.py)
enum : int8_t { IN_COMPLEX = 0, OUT_BOUND = 1, OUT_INFEASIBLE = 2, QP_FAILED = 3 };

// Tolerances (same values as alpha_dualqp.py)
constexpr double DEFAULT_PIV_TOL = 1e-10;   // native: singular working set
constexpr double DEFAULT_FEAS_TOL = 1e-10;  // native: primal feasibility
constexpr int DEFAULT_ITER_FACTOR = 10;     // native: max_iter = 10 n + 50
constexpr double DAQP_PRIMAL_TOL = 1e-10;   // DAQP, on the rescaled problem
constexpr double GRAM_FACTOR_TOL = 1e-15;   // pivoted Cholesky stopping rule
constexpr int DAQP_EQ = 5;                  // DAQP sense: active + immutable

// Phase timers for development only: compiled in with -DALPHADUAL_PROFILE
// (make profile). The normal build has no timers inside the computation.
#ifdef ALPHADUAL_PROFILE
#define PROF(label) prof_mark(label)
static void prof_mark(const char* label) {
    static auto last = std::chrono::steady_clock::now();
    const auto t = std::chrono::steady_clock::now();
    std::fprintf(stderr, "  [profile] %-28s %9.4f s\n", label,
                 std::chrono::duration<double>(t - last).count());
    last = t;
}
#else
#define PROF(label) do {} while (0)
#endif

[[noreturn]] static void die(const std::string& msg) {
    std::fprintf(stderr, "alphadual: %s\n", msg.c_str());
    std::exit(1);
}

// ---------------------------------------------------------------------------
// Input
// ---------------------------------------------------------------------------
struct Input {
    i64 N = 0, m = 0;
    vector<double> X;   // N x m, row-major
    vector<double> w;   // weights p(x) (zeros if unweighted)
};

static Input read_input(const char* path) {
    FILE* f = std::fopen(path, "rb");
    if (!f) die(std::string("cannot open input ") + path);
    char magic[8];
    i64 hdr[3];
    if (std::fread(magic, 1, 8, f) != 8 || std::memcmp(magic, "ALPHAIN1", 8) != 0)
        die("input is not an ALPHAIN1 file");
    if (std::fread(hdr, sizeof(i64), 3, f) != 3) die("truncated input header");
    Input in;
    in.N = hdr[0];
    in.m = hdr[1];
    const bool has_w = hdr[2] != 0;
    if (in.N < 0 || in.m <= 0) die("bad N or m in input");
    in.X.resize(size_t(in.N * in.m));
    if (std::fread(in.X.data(), sizeof(double), in.X.size(), f) != in.X.size())
        die("truncated point data");
    in.w.assign(size_t(in.N), 0.0);
    if (has_w && std::fread(in.w.data(), sizeof(double), in.w.size(), f) != in.w.size())
        die("truncated weights");
    std::fclose(f);
    return in;
}

// ---------------------------------------------------------------------------
// Graphs in CSR form: neighbours of x are idx[ptr[x] .. ptr[x+1]), sorted.
// ---------------------------------------------------------------------------
struct Graph {
    vector<i64> ptr;
    vector<i32> idx;
};

static Graph csr_from_lists(vector<vector<i32>>& lists) {
    Graph g;
    const size_t N = lists.size();
    g.ptr.assign(N + 1, 0);
    for (size_t i = 0; i < N; ++i) g.ptr[i + 1] = g.ptr[i] + i64(lists[i].size());
    g.idx.resize(size_t(g.ptr[N]));
    for (size_t i = 0; i < N; ++i) {
        std::sort(lists[i].begin(), lists[i].end());
        std::copy(lists[i].begin(), lists[i].end(), g.idx.begin() + g.ptr[i]);
        vector<i32>().swap(lists[i]);  // free memory as we go
    }
    return g;
}

// nanoflann adaptor over the live points (rows of X listed in ids)
struct Cloud {
    const double* X;
    i64 m;
    const vector<i32>* ids;
    size_t kdtree_get_point_count() const { return ids->size(); }
    double kdtree_get_pt(size_t i, size_t d) const { return X[size_t((*ids)[i]) * m + d]; }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};
using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, Cloud>, Cloud, -1, uint32_t>;

// Line 1 of Algorithm 1: Cech graph of the weighted ball cover at level a1
// (ball radius sqrt(a1 + p(x)); x ~ y iff the balls meet). a1 = inf: all pairs.
static bool use_kdtree(i64 K, i64 m);   // KD-tree or brute force (see below)

static Graph cech_graph(const Input& in, double a1, vector<char>& alive) {
    const i64 N = in.N, m = in.m;
    alive.assign(size_t(N), 0);
    vector<vector<i32>> nb(static_cast<size_t>(N));
    if (!std::isfinite(a1)) {
        for (i64 i = 0; i < N; ++i) {
            alive[i] = 1;
            nb[i].reserve(size_t(N - 1));
            for (i64 j = 0; j < N; ++j)
                if (j != i) nb[i].push_back(i32(j));
        }
        return csr_from_lists(nb);
    }
    vector<double> R(size_t(N), 0.0);
    vector<i32> ids;
    double Rmax = 0.0;
    for (i64 i = 0; i < N; ++i) {
        const double r2 = a1 + in.w[i];
        if (r2 >= 0) {
            alive[i] = 1;
            R[i] = std::sqrt(r2);
            Rmax = std::max(Rmax, R[i]);
            ids.push_back(i32(i));
        }
    }
    const double slack = 1.0 + 1e-9;  // extra candidates are harmless
    bool weighted = false;
    for (i32 i : ids) weighted |= (R[i] != Rmax);
    if (ids.size() >= 2 && use_kdtree(i64(ids.size()), m)) {
        Cloud cloud{in.X.data(), m, &ids};
        KDTree tree(size_t(m), cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
        tree.buildIndex();
        const double rad = 2.0 * Rmax * slack;
        std::vector<nanoflann::ResultItem<uint32_t, double>> res;
        nanoflann::SearchParameters sp;
        sp.sorted = false;
        for (size_t a = 0; a < ids.size(); ++a) {
            const i32 i = ids[a];
            res.clear();
            (void)tree.radiusSearch(&in.X[size_t(i) * m], rad * rad, res, sp);
            for (const auto& it : res) {
                const i32 j = ids[it.first];
                if (j == i) continue;
                if (weighted) {
                    const double t = (R[i] + R[j]) * slack;
                    if (it.second > t * t) continue;
                }
                nb[i].push_back(j);
            }
        }
    } else if (ids.size() >= 2) {
        // Brute force, vectorised: coordinates stored transposed (XT[d][a]),
        // so for a fixed point a the squared distances to all b > a are
        // accumulated one coordinate at a time over contiguous memory.
        const size_t K = ids.size();
        vector<double> XT(size_t(m) * K), acc(K);
        for (size_t a = 0; a < K; ++a)
            for (i64 d = 0; d < m; ++d) XT[size_t(d) * K + a] = in.X[size_t(ids[a]) * m + d];
        for (size_t a = 0; a + 1 < K; ++a) {
            const size_t len = K - a - 1;
            double* ac = &acc[a + 1];
            std::fill(ac, ac + len, 0.0);
            for (i64 d = 0; d < m; ++d) {
                const double xa = XT[size_t(d) * K + a];
                const double* col = &XT[size_t(d) * K + a + 1];
                for (size_t b = 0; b < len; ++b) {
                    const double t = col[b] - xa;
                    ac[b] += t * t;
                }
            }
            const i32 ia = ids[a];
            for (size_t b = a + 1; b < K; ++b) {
                const i32 ib = ids[b];
                const double t = (R[ia] + R[ib]) * slack;
                if (acc[b] <= t * t) {
                    nb[ia].push_back(ib);
                    nb[ib].push_back(ia);
                }
            }
        }
    }
    return csr_from_lists(nb);
}

// Neighbour search strategy. ALPHADUAL_FORCE_NN=kdtree|brute overrides it
// (used only to benchmark the two).
// Rule chosen from measurements (bunny, tori in 3D/50D/200D, uniform cubes
// in 6D/10D, cyclooctane in 24D): the KD-tree wins for large point sets in
// low/moderate dimension; brute force wins for small sets and high dimension
// (KD-trees degrade there), unless N is so large that its N^2 term dominates.
static bool use_kdtree(i64 K, i64 m) {
    if (const char* f = std::getenv("ALPHADUAL_FORCE_NN")) return std::string(f) == "kdtree";
    return (K >= 2000 && m <= 32) || (K >= 20000 && m <= 128);
}

// ---------------------------------------------------------------------------
// Per-vertex local data (lines 10-13): neighbours x_i, D = x_i - x, the Gram
// matrix B = D D^T (the only place the ambient dimension enters), U and c1.
// ---------------------------------------------------------------------------
struct Local {
    i64 n = 0;
    vector<double> D, B, U;   // D stored TRANSPOSED: D[d * n + i] = (x_i - x)_d
};

// B[i][j] = sum_d D[d][i] * D[d][j], accumulated row by row as m contiguous
// "axpy" updates (vectorised, sequential memory access). Each entry is the
// same sum in the same order d = 0..m-1, so B is exactly symmetric.
static void build_local(const Input& in, const Graph& g, i64 x, Local& L) {
    const i64 m = in.m, s = g.ptr[x], n = g.ptr[x + 1] - s;
    L.n = n;
    L.D.resize(size_t(n * m));
    L.B.resize(size_t(n * n));
    L.U.resize(size_t(n));
    const double* px = &in.X[size_t(x) * m];
    for (i64 i = 0; i < n; ++i) {
        const double* pi = &in.X[size_t(g.idx[s + i]) * m];
        for (i64 d = 0; d < m; ++d) L.D[d * n + i] = pi[d] - px[d];
    }
    for (i64 i = 0; i < n; ++i) {                              // line 11
        double* Bi = &L.B[i * n];
        const double di0 = L.D[i];
        const double* D0 = &L.D[0];
        for (i64 j = 0; j < n; ++j) Bi[j] = di0 * D0[j];
        for (i64 d = 1; d < m; ++d) {
            const double did = L.D[d * n + i];
            const double* Dd = &L.D[d * n];
            for (i64 j = 0; j < n; ++j) Bi[j] += did * Dd[j];
        }
    }
    for (i64 i = 0; i < n; ++i)                                // line 12
        L.U[i] = 0.5 * (in.w[g.idx[s + i]] - in.w[x] - L.B[i * n + i]);
}

// ---------------------------------------------------------------------------
// Native solver (port of dual_qp in alpha_dualqp.py)
// max -1/2 l^T B l + U^T l  s.t. l_i >= 0 for i not in J
// ---------------------------------------------------------------------------
struct NativeWork {
    i64 n = 0;
    vector<double> lam, L, b1, b2, g;
    vector<i64> W;
    vector<char> inW;
    vector<int8_t> kind;  // 0 inequality, 1 equality, 2 redundant equality
    // Grow-only buffers, reused across vertices. No zero-filling: dual_qp
    // resets lam/inW/kind itself and writes every L/b1/b2/W entry before
    // reading it (zero-filling the (n+1)^2 matrix L was a large waste).
    void resize(i64 n_) {
        n = n_;
        if (i64(lam.size()) < n) {
            lam.resize(size_t(n));
            g.resize(size_t(n));
            inW.resize(size_t(n));
            kind.resize(size_t(n));
        }
        if (i64(b1.size()) < n + 1) {
            b1.resize(size_t(n + 1));
            b2.resize(size_t(n + 1));
            W.resize(size_t(n + 1));
        }
        if (i64(L.size()) < (n + 1) * (n + 1)) L.resize(size_t((n + 1) * (n + 1)));
    }
};

// Cholesky of B[W,W]; returns -1, or the position of the first zero pivot.
static i64 chol(const double* B, i64 n, const i64* W, i64 w, double* L, i64 ld, double piv_tol) {
    for (i64 i = 0; i < w; ++i) {
        const i64 bi = W[i];
        for (i64 j = 0; j < i; ++j) {
            double s = B[bi * n + W[j]];
            for (i64 t = 0; t < j; ++t) s -= L[i * ld + t] * L[j * ld + t];
            L[i * ld + j] = s / L[j * ld + j];
        }
        double d = B[bi * n + bi];
        for (i64 t = 0; t < i; ++t) d -= L[i * ld + t] * L[i * ld + t];
        if (d <= piv_tol * B[bi * n + bi]) return i;
        L[i * ld + i] = std::sqrt(d);
    }
    return -1;
}

static void chol_solve(const double* L, i64 ld, i64 q, const double* b, double* z) {
    for (i64 i = 0; i < q; ++i) {
        double s = b[i];
        for (i64 t = 0; t < i; ++t) s -= L[i * ld + t] * z[t];
        z[i] = s / L[i * ld + i];
    }
    for (i64 i = q - 1; i >= 0; --i) {
        double s = z[i];
        for (i64 t = i + 1; t < q; ++t) s -= L[t * ld + i] * z[t];
        z[i] = s / L[i * ld + i];
    }
}

static i64 remove_at(i64* W, i64 w, i64 r) {
    for (i64 i = r; i < w - 1; ++i) W[i] = W[i + 1];
    return w - 1;
}

static int8_t dual_qp(const Local& Lc, const i32* J, int nJ, double c1, NativeWork& wk,
                      i64 max_iter, double piv_tol, double feas_tol, double& c_out,
                      i64& it_out) {
    const i64 n = Lc.n, ld = n + 1;
    const double* B = Lc.B.data();
    const double* U = Lc.U.data();
    double* lam = wk.lam.data();
    double* L = wk.L.data();
    double* b1 = wk.b1.data();
    double* b2 = wk.b2.data();
    i64* W = wk.W.data();
    char* inW = wk.inW.data();
    int8_t* kind = wk.kind.data();
    for (i64 i = 0; i < n; ++i) { lam[i] = 0.0; inW[i] = 0; kind[i] = 0; }
    i64 w = 0;
    for (int t = 0; t < nJ; ++t) { W[w++] = J[t]; inW[J[t]] = 1; kind[J[t]] = 1; }
    it_out = 0;
    c_out = 0.0;
    // Equality constraints: drop redundant ones, detect inconsistent ones.
    while (w > 0) {
        const i64 q = chol(B, n, W, w, L, ld, piv_tol);
        if (q < 0) break;
        for (i64 i = 0; i < q; ++i) b1[i] = B[W[i] * n + W[q]];
        chol_solve(L, ld, q, b1, b2);
        double s = U[W[q]], sabs = std::fabs(U[W[q]]);
        for (i64 i = 0; i < q; ++i) { s -= U[W[i]] * b2[i]; sabs += std::fabs(U[W[i]] * b2[i]); }
        if (std::fabs(s) > 1e-8 * sabs) { c_out = INFINITY; return OUT_INFEASIBLE; }
        kind[W[q]] = 2;
        inW[W[q]] = 0;
        w = remove_at(W, w, q);
    }
    i64 it = 0;
    while (it < max_iter) {
        ++it;
        it_out = it;
        const i64 q = chol(B, n, W, w, L, ld, piv_tol);
        if (q < 0) {
            for (i64 i = 0; i < w; ++i) b1[i] = U[W[i]];
            chol_solve(L, ld, w, b1, b2);  // b2 = lambda*_W
            double alpha = 1.0;
            i64 r = -1;
            for (i64 i = 0; i < w; ++i)
                if (kind[W[i]] == 0 && b2[i] < 0.0) {
                    const double li = lam[W[i]], a = li / (li - b2[i]);
                    if (a < alpha) { alpha = a; r = i; }
                }
            if (r < 0) {
                double c = 0.0;
                for (i64 i = 0; i < w; ++i) { lam[W[i]] = b2[i]; c += U[W[i]] * b2[i]; }
                c *= 0.5;
                if (c > c1) { c_out = c; return OUT_BOUND; }   // early exit
                // Gradient g = B lambda - U for all constraints, built from
                // ROWS of B (B is symmetric): contiguous, vectorised updates,
                // each g[k] accumulated in the same order i = 0..w-1.
                double* gv = wk.g.data();
                for (i64 k = 0; k < n; ++k) gv[k] = -U[k];
                for (i64 i = 0; i < w; ++i) {
                    const double li = b2[i];
                    const double* Bi = &B[W[i] * n];
                    for (i64 k = 0; k < n; ++k) gv[k] += Bi[k] * li;
                }
                double gmin = -feas_tol;
                i64 kmin = -1;
                for (i64 k = 0; k < n; ++k) {
                    if (inW[k] || kind[k] == 2) continue;
                    if (gv[k] < gmin) { gmin = gv[k]; kmin = k; }
                }
                if (kmin < 0) { c_out = c; return IN_COMPLEX; }
                W[w++] = kmin;
                inW[kmin] = 1;
            } else {
                for (i64 i = 0; i < w; ++i) {
                    const double li = lam[W[i]];
                    lam[W[i]] = li + alpha * (b2[i] - li);
                }
                lam[W[r]] = 0.0;
                inW[W[r]] = 0;
                w = remove_at(W, w, r);
            }
        } else {
            // Singular working set: step along a null direction p.
            for (i64 i = 0; i < q; ++i) b1[i] = B[W[i] * n + W[q]];
            chol_solve(L, ld, q, b1, b2);
            double gp = 0.0, pmax = 1.0;
            for (i64 i = 0; i <= q; ++i) {
                const double pi = i < q ? -b2[i] : 1.0;
                pmax = std::max(pmax, std::fabs(pi));
                const i64 bi = W[i];
                double gj = -U[bi];
                for (i64 t = 0; t < w; ++t) gj += B[bi * n + W[t]] * lam[W[t]];
                gp += pi * gj;
            }
            if (gp == 0.0) return QP_FAILED;
            const double sgn = gp < 0.0 ? 1.0 : -1.0;
            double tmax = INFINITY;
            i64 r = -1;
            for (i64 i = 0; i <= q; ++i) {
                const double pi = sgn * (i < q ? -b2[i] : 1.0);
                if (kind[W[i]] == 0 && pi < -1e-12 * pmax) {
                    const double t = lam[W[i]] / (-pi);
                    if (t < tmax) { tmax = t; r = i; }
                }
            }
            if (r < 0) { c_out = INFINITY; return OUT_INFEASIBLE; }   // dual unbounded
            for (i64 i = 0; i <= q; ++i) lam[W[i]] += tmax * sgn * (i < q ? -b2[i] : 1.0);
            lam[W[r]] = 0.0;
            inW[W[r]] = 0;
            w = remove_at(W, w, r);
        }
    }
    return QP_FAILED;
}

// ---------------------------------------------------------------------------
// DAQP: pivoted Cholesky B/s2 = G G^T (dot products only), then DAQP solves
// min 1/2|w|^2 s.t. G w <= -U/s2, whose dual is the paper's problem (10).
// ---------------------------------------------------------------------------
static i64 gram_factor(const vector<double>& B, i64 n, double scale, double tol, vector<double>& G) {
    vector<double> d(static_cast<size_t>(n));
    double dmax = 0.0;
    for (i64 i = 0; i < n; ++i) { d[i] = B[i * n + i] / scale; dmax = std::max(dmax, d[i]); }
    vector<char> used(size_t(n), 0);
    vector<double> Gf(size_t(n * std::max<i64>(n, 1)), 0.0);   // n x n, column k
    const i64 ldg = std::max<i64>(n, 1);
    i64 k = 0;
    while (k < n) {
        i64 p = -1;
        double best = tol * dmax;
        for (i64 i = 0; i < n; ++i)
            if (!used[i] && d[i] > best) { best = d[i]; p = i; }
        if (p < 0) break;
        used[p] = 1;
        const double piv = std::sqrt(d[p]);
        for (i64 i = 0; i < n; ++i) {
            if (used[i] && i != p) continue;
            double s = B[i * n + p] / scale;
            for (i64 t = 0; t < k; ++t) s -= Gf[i * ldg + t] * Gf[p * ldg + t];
            Gf[i * ldg + k] = s / piv;
        }
        for (i64 i = 0; i < n; ++i)
            if (!used[i]) d[i] -= Gf[i * ldg + k] * Gf[i * ldg + k];
        ++k;
    }
    const i64 kk = std::max<i64>(k, 1);
    G.assign(size_t(n * kk), 0.0);   // compact n x kk, row-major (DAQP's A)
    for (i64 i = 0; i < n; ++i)
        for (i64 t = 0; t < kk && t < ldg; ++t) G[i * kk + t] = Gf[i * ldg + t];
    return kk;
}

// ---------------------------------------------------------------------------
// One dimension of the main loop (lines 9-20) for all candidates.
// cand_ptr[x]..cand_ptr[x+1]: candidates owned by x; candJ: k positions each.
// ---------------------------------------------------------------------------
struct DimResult {
    vector<int8_t> status;
    vector<double> value;
    vector<i32> iters;
};

static DimResult process_dimension(const Input& in, const Graph& g, double a1,
                                   const vector<i64>& cand_ptr, const vector<i32>& candJ,
                                   int k, bool use_daqp) {
    const i64 N = in.N;
    const i64 C = cand_ptr[N];
    DimResult R;
    R.status.assign(size_t(C), QP_FAILED);
    R.value.assign(size_t(C), INFINITY);
    R.iters.assign(size_t(C), 0);
    Local Lc;
    NativeWork wk;
    vector<double> G, H, f, bu, bl, xbuf, lambuf;
    vector<int> sense;
    for (i64 x = 0; x < N; ++x) {
        const i64 s0 = cand_ptr[x], s1 = cand_ptr[x + 1];
        if (s0 == s1) continue;
        const double c1 = 0.5 * (a1 + in.w[x]);                  // line 13
        build_local(in, g, x, Lc);                               // lines 10-12
        const i64 n = Lc.n;
        if (!use_daqp) {
            double scale = 1e-300;
            for (i64 i = 0; i < n; ++i) scale = std::max(scale, Lc.B[i * n + i]);
            const double feas_tol = DEFAULT_FEAS_TOL * scale;
            const i64 max_iter = DEFAULT_ITER_FACTOR * n + 50;
            wk.resize(n);
            for (i64 c = s0; c < s1; ++c) {                      // line 14
                double cval;
                i64 its;
                const int8_t st = dual_qp(Lc, &candJ[size_t(c) * std::max(k, 1)], k, c1, wk,
                                          max_iter, DEFAULT_PIV_TOL, feas_tol, cval, its);
                R.status[c] = st;
                R.iters[c] = i32(its);
                if (st == IN_COMPLEX) R.value[c] = 2.0 * cval - in.w[x];   // line 19
            }
            continue;
        }
        // ---- DAQP ----
        if (n == 0) {   // isolated vertex: lambda = 0, c* = 0
            for (i64 c = s0; c < s1; ++c) {
                R.status[c] = c1 >= 0 ? IN_COMPLEX : OUT_BOUND;
                R.value[c] = -in.w[x];
            }
            continue;
        }
        double s2 = 0.0;
        for (i64 i = 0; i < n; ++i) s2 = std::max(s2, Lc.B[i * n + i]);
        if (s2 <= 0.0) s2 = 1.0;
        const i64 kk = gram_factor(Lc.B, n, s2, GRAM_FACTOR_TOL, G);
        H.assign(size_t(kk * kk), 0.0);
        for (i64 i = 0; i < kk; ++i) H[i * kk + i] = 1.0;
        f.assign(size_t(kk), 0.0);
        bu.resize(size_t(n));
        for (i64 i = 0; i < n; ++i) bu[i] = -Lc.U[i] / s2;
        bl.assign(size_t(n), -DAQP_INF);
        sense.assign(size_t(n), 0);
        xbuf.assign(size_t(kk), 0.0);
        lambuf.assign(size_t(n), 0.0);
        DAQPProblem qp;
        std::memset(&qp, 0, sizeof(qp));
        qp.n = int(kk);
        qp.m = int(n);
        qp.ms = 0;
        qp.H = H.data();
        qp.f = f.data();
        qp.A = G.data();
        qp.bupper = bu.data();
        qp.blower = bl.data();
        qp.sense = sense.data();
        // Same sequence as DAQP's Python bindings (Model.setup / settings /
        // update(sense=...) / solve), which is what alpha_dualqp.py validated.
        DAQPWorkspace* work = static_cast<DAQPWorkspace*>(std::calloc(1, sizeof(DAQPWorkspace)));
        if (!work) die("out of memory (DAQP workspace)");
        allocate_daqp_settings(work);
        double setup_time = 0.0;
        if (setup_daqp_main(&qp, work, &setup_time, 0) < 0) {
            std::free(work);   // leave the candidates marked QP_FAILED
            continue;
        }
        work->settings->fval_bound = std::min(c1 / s2, 1e30);
        work->settings->primal_tol = DAQP_PRIMAL_TOL;
        DAQPResult res;
        std::memset(&res, 0, sizeof(res));
        res.x = xbuf.data();
        res.lam = lambuf.data();
        for (i64 c = s0; c < s1; ++c) {                          // line 14
            const i32* J = &candJ[size_t(c) * std::max(k, 1)];
            for (int t = 0; t < k; ++t) sense[J[t]] = DAQP_EQ;
            daqp_update_ldp(DAQP_UPDATE_sense, work, &qp);
            daqp_solve(&res, work);                              // line 16
            for (int t = 0; t < k; ++t) sense[J[t]] = 0;
            R.iters[c] = res.iter;
            const double cval = res.fval * s2;
            const int flag = res.exitflag;
            if (flag >= 1 && cval <= c1) {                       // line 17
                R.status[c] = IN_COMPLEX;
                R.value[c] = 2.0 * cval - in.w[x];               // line 19
            } else {
                R.status[c] = (flag == 1 || flag == 2 || flag == -1) ? OUT_BOUND : QP_FAILED;
            }
        }
        free_daqp_workspace(work);
        free_daqp_ldp(work);
        std::free(work);
    }
    return R;
}

// ---------------------------------------------------------------------------
// Accepted simplices of one dimension k ("level"): rows of k+1 sorted vertex
// ids, rows sorted lexicographically, so the rows owned by vertex x (x = the
// smallest vertex) are rows ptr[x] .. ptr[x+1]-1.
// ---------------------------------------------------------------------------
struct Level {
    int k = 0;
    vector<i32> rows;     // n x (k+1)
    vector<double> vals;  // filtration values (squared radii)
    vector<i64> faces;    // k >= 1: n x (k+1), row in level k-1 of the face omitting vertex i
    vector<i64> ptr;      // N + 1
    i64 n() const { return i64(vals.size()); }
};

static void level_ptr(Level& L, i64 N) {
    const int w = L.k + 1;
    L.ptr.assign(size_t(N + 1), 0);
    for (i64 r = 0; r < L.n(); ++r) ++L.ptr[L.rows[size_t(r) * w] + 1];
    for (i64 x = 0; x < N; ++x) L.ptr[x + 1] += L.ptr[x];
}

// Binary search for the row `key` (w ids) among rows lo..hi-1 of S.
static i64 find_row(const Level& S, i64 lo, i64 hi, const i32* key, int w) {
    while (lo < hi) {
        const i64 mid = (lo + hi) / 2;
        const i32* r = &S.rows[size_t(mid) * w];
        int c = 0;
        for (int t = 0; t < w; ++t)
            if (r[t] != key[t]) { c = r[t] < key[t] ? -1 : 1; break; }
        if (c == 0) return mid;
        if (c < 0) lo = mid + 1; else hi = mid;
    }
    return -1;
}

// Candidate k-simplices (k >= 2): Lazy_{k-1}(X)_k, the sets of k+1 points all
// of whose (k-1)-faces are accepted. Each accepted (k-1)-simplex sigma =
// (x, y1, .., y_{k-1}) of S is extended by every accepted-edge neighbour
// z > y_{k-1} of y_{k-1}, and kept if the k faces "drop a vertex of sigma,
// add z" are in S too (binary searches). Candidates come out sorted
// lexicographically and grouped by owner x. Same rule as alpha_dualqp.py.
static void extend_candidates(i64 N, int k, const Level& S, const Graph& a,
                              const vector<i32>& aid, const Graph& c,
                              vector<i64>& cand_ptr, vector<i32>& candJ, vector<i32>& cand,
                              vector<i64>& faces) {
    cand_ptr.assign(size_t(N + 1), 0);
    candJ.clear(); cand.clear(); faces.clear();
    if (k == 2) {
        // Fast path for triangles (same candidates, same order): for owner x,
        // mark[y] = row of edge (x, y); aid[t] = row of the edge stored at
        // a.idx[t]. Every face lookup is then O(1) instead of a binary search.
        vector<i32> mark(static_cast<size_t>(N), -1);
        for (i64 x = 0; x < N; ++x) {
            const i64 s0 = S.ptr[x], s1 = S.ptr[x + 1];
            if (s0 == s1) continue;
            for (i64 e = s0; e < s1; ++e) mark[S.rows[2 * e + 1]] = i32(e);
            const i32* cn = c.idx.data() + c.ptr[x];
            const i32* ce = c.idx.data() + c.ptr[x + 1];
            for (i64 e = s0; e < s1; ++e) {
                const i32 y = S.rows[2 * e + 1];
                const i32 py = i32(std::lower_bound(cn, ce, y) - cn);
                for (i64 t = a.ptr[y]; t < a.ptr[y + 1]; ++t) {
                    const i32 z = a.idx[t];
                    if (z <= y || mark[z] < 0) continue;
                    cand.push_back(i32(x)); cand.push_back(y); cand.push_back(z);
                    candJ.push_back(py);
                    candJ.push_back(i32(std::lower_bound(cn, ce, z) - cn));
                    faces.push_back(aid[t]);      // drop x: edge (y, z)
                    faces.push_back(mark[z]);     // drop y: edge (x, z)
                    faces.push_back(e);           // drop z: edge (x, y)
                    ++cand_ptr[x + 1];
                }
            }
            for (i64 e = s0; e < s1; ++e) mark[S.rows[2 * e + 1]] = -1;
        }
        for (i64 x = 0; x < N; ++x) cand_ptr[x + 1] += cand_ptr[x];
        return;
    }
    vector<i32> key(static_cast<size_t>(k));
    vector<i64> fidx(static_cast<size_t>(k + 1));
    for (i64 row = 0; row < S.n(); ++row) {
        const i32* sg = &S.rows[size_t(row) * k];
        const i32 x = sg[0], last = sg[k - 1];
        for (i64 t = a.ptr[last]; t < a.ptr[last + 1]; ++t) {
            const i32 z = a.idx[t];
            if (z <= last) continue;
            bool ok = true;
            for (int i = 0; i < k && ok; ++i) {        // face: drop sg[i], add z
                int q = 0;
                for (int u = 0; u < k; ++u)
                    if (u != i) key[q++] = sg[u];
                key[k - 1] = z;                        // z is largest: still sorted
                const i64 j = find_row(S, S.ptr[key[0]], S.ptr[key[0] + 1], key.data(), k);
                if (j < 0) ok = false; else fidx[i] = j;
            }
            if (!ok) continue;
            fidx[k] = row;                             // dropping z gives sigma
            const i32* cn = c.idx.data() + c.ptr[x];
            const i32* ce = c.idx.data() + c.ptr[x + 1];
            for (int u = 0; u < k; ++u) cand.push_back(sg[u]);
            cand.push_back(z);
            for (int u = 1; u < k; ++u) candJ.push_back(i32(std::lower_bound(cn, ce, sg[u]) - cn));
            candJ.push_back(i32(std::lower_bound(cn, ce, z) - cn));
            for (int u = 0; u <= k; ++u) faces.push_back(fidx[u]);
            ++cand_ptr[x + 1];
        }
    }
    for (i64 x = 0; x < N; ++x) cand_ptr[x + 1] += cand_ptr[x];
}

// ---------------------------------------------------------------------------
struct Output {
    int width = 1;               // max_dim + 1
    vector<i64> counts;          // simplices per dimension
    vector<i32> simp;            // M x width, -1 padded
    vector<int8_t> dims;         // M
    vector<double> vals;         // M
};

static double now_s() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

// Output format 2: "ALPHOUT2" | int64 M | int64 width | int64 counts[width] |
// int32 simplices[M*width] | int8 dims[M] | double values[M]
static void write_output(const char* path, const Output& o) {
    FILE* f = std::fopen(path, "wb");
    if (!f) die(std::string("cannot write ") + path);
    const i64 hdr[2] = {i64(o.dims.size()), i64(o.width)};
    std::fwrite("ALPHOUT2", 1, 8, f);
    std::fwrite(hdr, sizeof(i64), 2, f);
    std::fwrite(o.counts.data(), sizeof(i64), o.counts.size(), f);
    std::fwrite(o.simp.data(), sizeof(i32), o.simp.size(), f);
    std::fwrite(o.dims.data(), 1, o.dims.size(), f);
    std::fwrite(o.vals.data(), sizeof(double), o.vals.size(), f);
    if (std::fclose(f) != 0) die(std::string("error writing ") + path);
}

static void dim_stats_json(FILE* f, int k, const DimResult& r, bool daqp, bool last) {
    i64 acc = 0, bnd = 0, inf = 0, fail = 0, its = 0;
    for (size_t i = 0; i < r.status.size(); ++i) {
        acc += r.status[i] == IN_COMPLEX;
        bnd += r.status[i] == OUT_BOUND;
        inf += r.status[i] == OUT_INFEASIBLE;
        fail += r.status[i] == QP_FAILED;
        its += r.iters[i];
    }
    std::fprintf(f, "  \"dim%d\": {\"candidates\": %lld, \"accepted\": %lld, \"rejected\": %lld, ",
                 k, (long long)r.status.size(), (long long)acc, (long long)(bnd + inf));
    if (daqp)
        std::fprintf(f, "\"early_terminated\": null, \"infeasible\": null, ");
    else
        std::fprintf(f, "\"early_terminated\": %lld, \"infeasible\": %lld, ",
                     (long long)bnd, (long long)inf);
    std::fprintf(f, "\"failed\": %lld, \"iterations\": %lld}%s\n", (long long)fail,
                 (long long)its, last ? "" : ",");
}

#ifndef ALPHADUAL_BUILD_FLAGS
#define ALPHADUAL_BUILD_FLAGS "unknown"
#endif
// Compiler name: __VERSION__ exists in GCC/Clang only; MSVC has _MSC_FULL_VER.
#define ALPHADUAL_STR2(x) #x
#define ALPHADUAL_STR(x) ALPHADUAL_STR2(x)
#if defined(__VERSION__)
#define ALPHADUAL_COMPILER __VERSION__
#elif defined(_MSC_FULL_VER)
#define ALPHADUAL_COMPILER "MSVC " ALPHADUAL_STR(_MSC_FULL_VER)
#else
#define ALPHADUAL_COMPILER "unknown"
#endif
static const char* VERSION_STRING =
    "alphadual 1.1 (Carlsson & Carlsson Algorithm 1, any max dimension) | DAQP v0.10.3"
    " | nanoflann 1.14.0 | compiler " ALPHADUAL_COMPILER " | flags " ALPHADUAL_BUILD_FLAGS;

int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--version") {
        std::printf("%s\n", VERSION_STRING);
        return 0;
    }
    if (argc < 3) {
        std::fprintf(stderr,
                     "usage: alphadual INPUT OUTPUT (--r R | --a1 A1) [--max-dim D] "
                     "[--solver daqp|native] [--stats FILE.json]\n       alphadual --version\n");
        return 2;
    }
    const char* in_path = argv[1];
    const char* out_path = argv[2];
    double a1 = NAN;
    int max_dim = 2;
    std::string solver = "daqp", stats_path = std::string(out_path) + ".json";
    for (int i = 3; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() -> const char* {
            if (i + 1 >= argc) die("missing value after " + a);
            return argv[++i];
        };
        if (a == "--r") { const double r = std::strtod(next(), nullptr); a1 = r * r; }
        else if (a == "--a1") a1 = std::strtod(next(), nullptr);
        else if (a == "--max-dim") max_dim = std::atoi(next());
        else if (a == "--solver") solver = next();
        else if (a == "--stats") stats_path = next();
        else die("unknown option " + a);
    }
    if (std::isnan(a1)) die("give --r or --a1");
    if (max_dim < 0) die("--max-dim must be >= 0");
    if (solver != "daqp" && solver != "native") die("solver must be daqp or native");
    const bool use_daqp = solver == "daqp";
    if (!std::isfinite(a1) && !use_daqp)
        die("the complete filtration (r = inf) needs --solver daqp");

    // ===================== timed region: read + compute =====================
    const double t0 = now_s();
    PROF("start");
    const Input in = read_input(in_path);
    const double t_read = now_s();
    PROF("read input");
    const i64 N = in.N;

    vector<char> alive;
    const Graph g = cech_graph(in, a1, alive);                   // line 1
    PROF("cech graph");

    vector<Level> lev(static_cast<size_t>(max_dim + 1));
    vector<DimResult> res(static_cast<size_t>(max_dim + 1));
    for (int k = 0; k <= max_dim; ++k) lev[k].k = k;

    // k = 0: vertices of G
    {
        vector<i64> cptr(size_t(N + 1), 0);
        for (i64 x = 0; x < N; ++x) cptr[x + 1] = cptr[x] + (alive[x] ? 1 : 0);
        vector<i32> candJ0(size_t(cptr[N]), 0);
        res[0] = process_dimension(in, g, a1, cptr, candJ0, 0, use_daqp);
        for (i64 x = 0; x < N; ++x)
            if (alive[x] && res[0].status[cptr[x]] == IN_COMPLEX) {
                lev[0].rows.push_back(i32(x));
                lev[0].vals.push_back(res[0].value[cptr[x]]);
            }
        level_ptr(lev[0], N);
        PROF("dim 0 QPs");
    }
    vector<i64> vpos(size_t(N), -1);   // row of vertex v in level 0 (-1: none)
    for (i64 i = 0; i < lev[0].n(); ++i) vpos[lev[0].rows[i]] = i;

    // k = 1: edges of G, owned by their smaller vertex
    Graph ag;                          // graph of accepted edges (for k >= 2)
    vector<i32> aid;                   // aid[t] = row in level 1 of edge at ag.idx[t]
    if (max_dim >= 1) {
        vector<i64> eptr(size_t(N + 1), 0);
        vector<i32> eJ, ey;
        for (i64 x = 0; x < N; ++x) {
            for (i64 t = g.ptr[x]; t < g.ptr[x + 1]; ++t)
                if (g.idx[t] > x) { eJ.push_back(i32(t - g.ptr[x])); ey.push_back(g.idx[t]); }
            eptr[x + 1] = i64(eJ.size());
        }
        PROF("edge candidates");
        res[1] = process_dimension(in, g, a1, eptr, eJ, 1, use_daqp);
        PROF("dim 1 QPs");
        Level& E = lev[1];
        for (i64 x = 0; x < N; ++x)
            for (i64 c = eptr[x]; c < eptr[x + 1]; ++c)
                if (res[1].status[c] == IN_COMPLEX) {
                    E.rows.push_back(i32(x));
                    E.rows.push_back(ey[c]);
                    E.vals.push_back(res[1].value[c]);
                    E.faces.push_back(vpos[ey[c]]);    // face omitting x
                    E.faces.push_back(vpos[x]);        // face omitting y
                }
        level_ptr(E, N);
        // Accepted-edge graph in CSR form by counting sort; lists come out
        // sorted because E is sorted by (x, y).
        ag.ptr.assign(size_t(N + 1), 0);
        for (i64 e = 0; e < E.n(); ++e) { ++ag.ptr[E.rows[2 * e] + 1]; ++ag.ptr[E.rows[2 * e + 1] + 1]; }
        for (i64 x = 0; x < N; ++x) ag.ptr[x + 1] += ag.ptr[x];
        ag.idx.resize(size_t(ag.ptr[N]));
        aid.resize(size_t(ag.ptr[N]));
        vector<i64> fill(ag.ptr.begin(), ag.ptr.end() - 1);
        for (i64 e = 0; e < E.n(); ++e) {
            const i64 p = fill[E.rows[2 * e + 1]]++;
            ag.idx[p] = E.rows[2 * e]; aid[p] = i32(e);
        }
        for (i64 e = 0; e < E.n(); ++e) {
            const i64 p = fill[E.rows[2 * e]]++;
            ag.idx[p] = E.rows[2 * e + 1]; aid[p] = i32(e);
        }
    }

    // k >= 2: Lazy_{k-1}(X)_k from the accepted (k-1)-simplices
    for (int k = 2; k <= max_dim; ++k) {
        vector<i64> cptr, cfaces;
        vector<i32> cJ, cand;
        extend_candidates(N, k, lev[k - 1], ag, aid, g, cptr, cJ, cand, cfaces);
        PROF("candidates (k >= 2)");
        res[k] = process_dimension(in, g, a1, cptr, cJ, k, use_daqp);
        PROF("QPs (k >= 2)");
        Level& L = lev[k];
        for (i64 c = 0; c < cptr[N]; ++c)
            if (res[k].status[c] == IN_COMPLEX) {
                L.rows.insert(L.rows.end(), &cand[size_t(c) * (k + 1)], &cand[size_t(c) * (k + 1)] + k + 1);
                L.faces.insert(L.faces.end(), &cfaces[size_t(c) * (k + 1)],
                               &cfaces[size_t(c) * (k + 1)] + k + 1);
                L.vals.push_back(res[k].value[c]);
            }
        level_ptr(L, N);
    }

    // Monotonicity guard (a face never enters after its cofaces), dimension by
    // dimension upwards so that corrections propagate; then the global order.
    for (int k = 1; k <= max_dim; ++k)
        for (i64 i = 0; i < lev[k].n(); ++i)
            for (int u = 0; u <= k; ++u)
                lev[k].vals[i] = std::max(lev[k].vals[i],
                                          lev[k - 1].vals[lev[k].faces[size_t(i) * (k + 1) + u]]);
    Output out;
    out.width = max_dim + 1;
    i64 M = 0;
    for (int k = 0; k <= max_dim; ++k) { out.counts.push_back(lev[k].n()); M += lev[k].n(); }
    // Sort compact records (value, dim, location of the vertex ids); vertex
    // ids are only read on ties. Same total order as numpy.lexsort in Python.
    struct Rec { double v; int8_t d; i64 row; };
    vector<Rec> recs;
    recs.reserve(size_t(M));
    for (int k = 0; k <= max_dim; ++k)
        for (i64 i = 0; i < lev[k].n(); ++i) recs.push_back(Rec{lev[k].vals[i], int8_t(k), i});
    std::sort(recs.begin(), recs.end(), [&](const Rec& p, const Rec& q) {
        if (p.v != q.v) return p.v < q.v;
        if (p.d != q.d) return p.d < q.d;
        const i32* a = &lev[p.d].rows[size_t(p.row) * (p.d + 1)];
        const i32* b = &lev[q.d].rows[size_t(q.row) * (q.d + 1)];
        for (int t = 0; t <= p.d; ++t)
            if (a[t] != b[t]) return a[t] < b[t];
        return false;
    });
    out.simp.assign(size_t(M) * out.width, -1);
    out.dims.resize(size_t(M));
    out.vals.resize(size_t(M));
    for (i64 i = 0; i < M; ++i) {
        const Rec& r = recs[i];
        const i32* s = &lev[r.d].rows[size_t(r.row) * (r.d + 1)];
        std::copy(s, s + r.d + 1, &out.simp[size_t(i) * out.width]);
        out.dims[i] = r.d;
        out.vals[i] = r.v;
    }
    PROF("guard + sort + output arrays");
    const double t1 = now_s();
    // ====================== end of timed region =============================

    write_output(out_path, out);
    const double t2 = now_s();

    i64 maxdeg = 0;
    for (i64 x = 0; x < N; ++x) maxdeg = std::max(maxdeg, g.ptr[x + 1] - g.ptr[x]);
    FILE* f = std::fopen(stats_path.c_str(), "w");
    if (!f) die("cannot write " + stats_path);
    std::fprintf(f, "{\n  \"implementation\": \"cpp\",\n  \"solver\": \"%s\",\n", solver.c_str());
    std::fprintf(f, "  \"N\": %lld,\n  \"m\": %lld,\n  \"max_dim\": %d,\n", (long long)N,
                 (long long)in.m, max_dim);
    if (std::isfinite(a1)) std::fprintf(f, "  \"a1\": %.17g,\n", a1);
    else std::fprintf(f, "  \"a1\": \"inf\",\n");
    std::fprintf(f, "  \"alpha_time_s\": %.9f,\n  \"read_time_s\": %.9f,\n", t1 - t0, t_read - t0);
    std::fprintf(f, "  \"compute_time_s\": %.9f,\n  \"write_time_s\": %.9f,\n", t1 - t_read, t2 - t1);
    std::fprintf(f, "  \"n_simplices\": [");
    for (int k = 0; k <= max_dim; ++k) std::fprintf(f, "%s%lld", k ? ", " : "", (long long)lev[k].n());
    std::fprintf(f, "],\n");
    auto cnt = [&](int k) { return (long long)(k <= max_dim ? lev[k].n() : 0); };
    std::fprintf(f, "  \"n_vertices\": %lld,\n  \"n_edges\": %lld,\n  \"n_triangles\": %lld,\n",
                 cnt(0), cnt(1), cnt(2));
    // largest filtration value (last of the ordered filtration): lets the
    // Python side report the largest entry radius without reading the output
    if (out.vals.empty()) std::fprintf(f, "  \"max_value\": null,\n");
    else std::fprintf(f, "  \"max_value\": %.17g,\n", out.vals.back());
    std::fprintf(f, "  \"cech_edges\": %lld,\n  \"max_cech_degree\": %lld,\n",
                 (long long)(g.idx.size() / 2), (long long)maxdeg);
    for (int k = 0; k <= max_dim; ++k) dim_stats_json(f, k, res[k], use_daqp, k == max_dim);
    std::fprintf(f, "}\n");
    std::fclose(f);
    return 0;
}
