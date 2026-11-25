"""
alpha_dualqp.py
===============

Alpha complex up to dimension 2 (vertices, edges, triangles) computed with the
dual active-set quadratic-programming approach of

    E. Carlsson and J. Carlsson,
    "Computing the alpha complex using dual active set quadratic programming",
    arXiv:2310.00536 (2023).

This is a Python/Numba re-implementation of Algorithm 1 of the paper,
restricted to d <= 2:

  1. Build the Cech graph G of the (weighted) ball cover at level a1.
  2. For k = 0, 1, 2:
       * candidate k-simplices: vertices of G (k=0), edges of G (k=1),
         3-cliques of the alpha edges found so far, i.e. Lazy_1(X)_2 (k=2);
       * for every vertex x: build B = ((x_i-x).(x_j-x)) and
         U = 1/2 (p(x_i) - p(x) - |x_i-x|^2) over the Cech neighbours x_i of x,
         and c1 = (a1 + p(x)) / 2 (lines 10-13 of Algorithm 1);
       * for every candidate sigma = [x, x_j1, ..., x_jk] with x its smallest
         vertex, solve the dual QP (10) with the equality set J = {j1..jk},
         stopping early as soon as a dual-feasible value exceeds c1.
         sigma is in the complex iff c* <= c1; its filtration value is
         w(sigma) = 2c* - p(x) and its witness is x - sum_j lambda_j (x_j - x).

Two interchangeable solvers for the dual QP (10) (argument `solver`):

* "daqp" (default): DAQP, the dual active-set solver of Arnstrom, Bemporad &
  Axehill (IEEE TAC 67(8), 2022) -- reference [2] of the paper, the solver the
  authors used -- through its official Python bindings (pip install daqp).
  DAQP takes a primal problem min 1/2|w|^2 s.t. G w <= -U, and works on its
  dual, whose Hessian is G G^T. Choosing G as a (pivoted) Cholesky factor of
  B = G G^T makes DAQP's internal dual exactly problem (10): the solver sees
  only dot products, never the ambient dimension. Equality constraints J use
  DAQP's equality flag, and the early termination of the paper is DAQP's
  `fval_bound` setting (= c1). Each local problem is rescaled to unit size so
  that DAQP's tolerances are relative. One DAQP workspace is set up per
  vertex x (lines 10-13) and only the equality flags change per simplex.
  Note: DAQP reports "bound exceeded" and "infeasible" with the same exit
  code, so these are counted together as rejected simplices.

* "native": a Numba re-implementation of the same dual active-set scheme
  (iterates always dual feasible, so the dual objective is a valid lower bound
  and early termination is exact; singular working sets handled by stepping
  along a null-space direction, which also detects primal infeasibility).
  It avoids the per-QP Python->C call overhead of the DAQP bindings and
  recomputes small factorizations instead of DAQP's recursive LDL^T updates.

Both give identical complexes (checked by validate_against_gudhi.py).

Conventions: unweighted alpha complex Alpha(S, r) = Alpha(S, p=0, a1=r^2);
filtration values are squared radii (same convention as GUDHI).


READING GUIDE (the math behind the code, in the paper's notation)
-----------------------------------------------------------------
Power distance of a point y to a site x with weight p(x):
    pi_x(y) = |y - x|^2 - p(x).            (p = 0: squared Euclidean distance)

A set sigma = {x, x_j1, ..., x_jk} is a simplex of the alpha complex at level
a1 iff some point y is "equally close" (in power distance) to all points of
sigma, no other point is closer, and that common value is <= a1. The smallest
such value is the filtration value (entry "time") of sigma: for p = 0 it is
the squared radius of the smallest empty sphere through the points of sigma.
Finding it is the quadratic program (13) of the paper:

    minimize  pi_x(y)
    subject to pi_x(y) =  pi_z(y)  for z in sigma \\ {x}      (equalities, J)
               pi_x(y) <= pi_z(y)  for every other point z  (inequalities)

Expanding the squares, every constraint is LINEAR in y:
    2 (z - x).y <= |z|^2 - |x|^2 - p(z) + p(x)   i.e.  A_i y <= V_i.
So (13) is the primal QP (9). Its Lagrangian dual is (10):

    maximize   -1/2 lambda^T B lambda + U^T lambda
    subject to lambda_i >= 0 for the inequality constraints (i not in J)

with  B_ij = (x_i - x).(x_j - x)  (a Gram matrix: only dot products!) and
      U_i  = 1/2 (p(x_i) - p(x) - |x_i - x|^2).

Why the dual? Weak duality: ANY feasible lambda gives a LOWER bound on the
primal optimum. The primal optimum is c* = (w(sigma) + p(x)) / 2, so as soon
as an iterate's dual value exceeds c1 = (a1 + p(x)) / 2 we know w(sigma) > a1
and sigma can be thrown away immediately ("early termination"). lambda = 0 is
always feasible, so no start-up work is needed. Most candidates are rejected
after only a handful of iterations, touching only a few nearby points (the
"active set"), e.g. the two points r, s of Figure 3 in the paper.

At the optimum, the KKT conditions give the witness (the point y realising
the minimum):   y* = x - sum_i lambda*_i (x_i - x).

Index bookkeeping used throughout:
  * Points are numbered 0..N-1 (rows of P).
  * Graphs are stored in CSR form (compressed sparse rows): the neighbours of
    vertex x are idx[ptr[x]:ptr[x+1]], sorted increasingly.
  * Every candidate simplex is "owned" by its smallest vertex x. Candidates
    are grouped by owner: those of x are rows cand_ptr[x]:cand_ptr[x+1] of the
    candidate arrays. This is what lets B, U be built ONCE per vertex and
    reused for all simplices owned by it (the paper's main-loop trick).
  * cand_J stores the OTHER vertices of a candidate not by point number but
    by POSITION in x's neighbour list (= the row/column index in B and U),
    because that is how the QP addresses its constraints.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix      # sparse matrices, used to build graphs
from scipy.spatial import cKDTree        # fast neighbour search in low dimension

# DAQP (the dual active-set solver of ref. [2]) is optional: without it, only
# the "native" solver is available.
try:
    import daqp
    HAVE_DAQP = True
except ImportError:
    daqp = None
    HAVE_DAQP = False

# Numba compiles the decorated functions (@njit) to machine code the first
# time they are called. Without Numba, the same code runs as plain Python:
# correct, but 100x+ slower. The fallback below defines a do-nothing @njit.
try:
    from numba import njit
    HAVE_NUMBA = True
except ImportError:  # pure-Python fallback (works, but is very slow)
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        # Supports both "@njit" and "@njit(cache=True)" usages.
        if len(args) == 1 and callable(args[0]) and not kwargs:
            return args[0]

        def deco(f):
            return f
        return deco


# QP outcome codes (one per candidate simplex)
IN_COMPLEX = 0       # optimum found and c* <= c1
OUT_BOUND = 1        # early termination: a dual-feasible value exceeded c1
OUT_INFEASIBLE = 2   # primal infeasible (dual unbounded): empty Voronoi face
QP_FAILED = 3        # iteration limit / numerical breakdown

# Tolerances and settings. "relative" = multiplied by the local length scale
# (the largest |x_i - x|^2 among the neighbours of x).
DEFAULT_PIV_TOL = 1e-10    # relative pivot tolerance for singular working sets
DEFAULT_FEAS_TOL = 1e-10   # relative primal-feasibility tolerance
DEFAULT_ITER_FACTOR = 10   # max_iter = DEFAULT_ITER_FACTOR * n_neighbours + 50
DEFAULT_SOLVER = "daqp"
SOLVERS = ("daqp", "native")
DAQP_PRIMAL_TOL = 1e-10    # DAQP primal tolerance on the rescaled problem
GRAM_FACTOR_TOL = 1e-15    # pivoted Cholesky stops when residual diag <= tol*max


# ---------------------------------------------------------------------------
# Small dense linear algebra on the working set (indices W[0:w] of B)
#
# The native solver repeatedly needs to solve small linear systems
# B[W, W] z = b, where W is the current "working set" (the constraints
# treated as equalities). W is tiny (a handful of indices), so a plain
# Cholesky factorisation B[W,W] = L L^T, recomputed from scratch, is cheap.
# ---------------------------------------------------------------------------
@njit(cache=True)
def _chol(B, W, w, L, piv_tol):
    """Cholesky of B[W,W]. Returns -1 on success, otherwise the position q of
    the first pivot that is (numerically) zero; rows 0..q-1 of L are valid."""
    # Standard row-by-row Cholesky: row i of L is computed from rows 0..i-1.
    for i in range(w):
        bi = W[i]
        # Off-diagonal entries L[i, j], j < i.
        for j in range(i):
            s = B[bi, W[j]]
            for t in range(j):
                s -= L[i, t] * L[j, t]
            L[i, j] = s / L[j, j]
        # Diagonal entry: d = B_ii - sum_t L[i,t]^2 = squared length of the
        # part of constraint i that is NOT explained by constraints 0..i-1.
        d = B[bi, bi]
        for t in range(i):
            d -= L[i, t] * L[i, t]
        # If that remainder is (relatively) zero, constraint i is linearly
        # dependent on the previous ones: B[W,W] is singular. Report where.
        if d <= piv_tol * B[bi, bi]:
            return i
        L[i, i] = np.sqrt(d)
    return -1


@njit(cache=True)
def _chol_solve(L, q, b, z):
    """Solve (L L^T) z = b using the leading q x q block of L."""
    # Forward substitution: solve L y = b (y stored in z).
    for i in range(q):
        s = b[i]
        for t in range(i):
            s -= L[i, t] * z[t]
        z[i] = s / L[i, i]
    # Back substitution: solve L^T z = y.
    for i in range(q - 1, -1, -1):
        s = z[i]
        for t in range(i + 1, q):
            s -= L[t, i] * z[t]
        z[i] = s / L[i, i]


@njit(cache=True)
def _remove(W, w, r):
    # Delete position r from the list W[0:w] (shift the tail left by one) and
    # return the new length.
    for i in range(r, w - 1):
        W[i] = W[i + 1]
    return w - 1


# ---------------------------------------------------------------------------
# Dual active-set QP:  max  -1/2 l^T B l + U^T l   s.t. l_i >= 0 (i not in J)
#
# This is the "native" solver. Internally it MINIMISES f(l) = 1/2 l^T B l
# - U^T l (the negative of the dual objective). It is a primal active-set
# method applied to the dual, as in DAQP:
#   * W ("working set") = indices whose l_i is free to be nonzero; all other
#     l_i are held at 0. Equality constraints (J) are always in W.
#   * Each iteration solves the small system B[W,W] l_W = U_W (the best l
#     using only the constraints in W), then either
#       - moves towards that solution, dropping an index whose l_i would turn
#         negative (it was not really active), or
#       - if the solution is already feasible: checks whether some primal
#         constraint outside W is violated; if so, adds the most violated one
#         to W, otherwise we are optimal.
#   * The iterate l is always dual feasible, so the dual value is always a
#     valid lower bound on the answer -> early termination when it beats c1.
#
# Interpretation of the gradient g = B l - U: for constraint k,
#   g_k = V_k - A_k y   with y = x - sum l_i (x_i - x),
# i.e. g_k < 0 means the current witness y violates constraint k (point x_k
# is strictly closer to y than x is).
# ---------------------------------------------------------------------------
@njit(cache=True)
def dual_qp(B, U, J, nJ, c1, lam, W, inW, kind, L, b1, b2,
            max_iter, piv_tol, feas_tol):
    """Returns (status, c, iterations); on IN_COMPLEX, lam holds lambda*.

    kind[i]: 0 = inequality, 1 = equality (in J), 2 = redundant equality."""
    n = U.shape[0]
    # Reset the work arrays (they are reused across all candidates of x).
    for i in range(n):
        lam[i] = 0.0
        inW[i] = False
        kind[i] = 0
    # The equality constraints (the other vertices of sigma) start in W.
    w = 0
    for t in range(nJ):
        j = J[t]
        W[w] = j
        w += 1
        inW[j] = True
        kind[j] = 1

    # Equality constraints: drop redundant ones, detect inconsistent ones.
    # If the equalities themselves are linearly dependent (e.g. 3 collinear
    # points), then either they contradict each other -> no point is
    # equidistant from all of sigma -> sigma is not a simplex, or one of them
    # is implied by the others -> it can simply be dropped.
    while w > 0:
        q = _chol(B, W, w, L, piv_tol)
        if q < 0:
            break                              # independent: nothing to do
        # Express constraint W[q] as a combination of W[0..q-1] (null vector
        # p = (-b2, 1) of B[W,W]) and test whether the right-hand sides agree.
        for i in range(q):
            b1[i] = B[W[i], W[q]]
        _chol_solve(L, q, b1, b2)
        s = U[W[q]]
        sabs = abs(U[W[q]])
        for i in range(q):
            s -= U[W[i]] * b2[i]
            sabs += abs(U[W[i]] * b2[i])
        if abs(s) > 1e-8 * sabs:
            return OUT_INFEASIBLE, np.inf, 0   # contradictory equalities
        kind[W[q]] = 2                         # redundant: drop it for good
        inW[W[q]] = False
        w = _remove(W, w, q)

    it = 0
    while it < max_iter:
        it += 1
        q = _chol(B, W, w, L, piv_tol)
        if q < 0:
            # ---- Case 1: B[W,W] is invertible -------------------------------
            # Equality-constrained subproblem on the working set.
            # Its solution lambda*_W solves B[W,W] lambda*_W = U[W].
            for i in range(w):
                b1[i] = U[W[i]]
            _chol_solve(L, w, b1, b2)          # b2 = lambda*_W
            # Would some inequality multiplier become negative? Then we can
            # only move a fraction alpha of the way (until the first one hits
            # zero), remembering which one (r) "blocks".
            alpha = 1.0
            r = -1
            for i in range(w):
                if kind[W[i]] == 0 and b2[i] < 0.0:
                    li = lam[W[i]]
                    a = li / (li - b2[i])
                    if a < alpha:
                        alpha = a
                        r = i
            if r < 0:
                # Full step is feasible: accept lambda*_W.
                c = 0.0
                for i in range(w):
                    lam[W[i]] = b2[i]
                    c += U[W[i]] * b2[i]
                # At the subproblem optimum, lambda^T B lambda = U^T lambda, so
                # the dual value -1/2 l^T B l + U^T l equals 1/2 U^T lambda.
                c *= 0.5                       # dual value at subproblem optimum
                if c > c1:                     # weak duality: early exit
                    return OUT_BOUND, c, it
                # Find the most violated primal constraint among those not in
                # W (most negative g_k = (B lambda)_k - U_k).
                gmin = -feas_tol
                kmin = -1
                for k in range(n):             # most violated primal constraint
                    if inW[k] or kind[k] == 2:
                        continue
                    g = -U[k]
                    for i in range(w):
                        g += B[k, W[i]] * b2[i]
                    if g < gmin:
                        gmin = g
                        kmin = k
                if kmin < 0:
                    # Nothing violated: lambda is optimal, c is the optimum.
                    return IN_COMPLEX, c, it
                # Otherwise make that constraint active and iterate again.
                W[w] = kmin
                w += 1
                inW[kmin] = True
            else:
                # Partial step: move a fraction alpha towards lambda*_W, then
                # drop the blocking constraint r from the working set.
                for i in range(w):
                    li = lam[W[i]]
                    lam[W[i]] = li + alpha * (b2[i] - li)
                lam[W[r]] = 0.0
                inW[W[r]] = False
                w = _remove(W, w, r)
        else:
            # ---- Case 2: B[W,W] is singular ---------------------------------
            # The constraint just added (at position q) is a linear
            # combination of the others: more active constraints than the
            # local affine dimension. The objective is then LINEAR along a
            # null direction p of B[W,W], so we can move along p until some
            # multiplier hits zero (and drop it). If no multiplier ever
            # blocks, the dual is unbounded: the primal is infeasible, i.e.
            # the Voronoi face of sigma is empty.
            # Singular working set: move along a null direction p of B[W,W]
            # (p = 1 on W[q], -B^{-1} B[.,W[q]] on the rest), linear in l.
            for i in range(q):
                b1[i] = B[W[i], W[q]]
            _chol_solve(L, q, b1, b2)
            # gp = directional derivative of f along p (g^T p).
            gp = 0.0
            pmax = 1.0
            for i in range(q + 1):
                pi = -b2[i] if i < q else 1.0
                if abs(pi) > pmax:
                    pmax = abs(pi)
                bi = W[i]
                gj = -U[bi]
                for t in range(w):
                    gj += B[bi, W[t]] * lam[W[t]]
                gp += pi * gj
            if gp == 0.0:
                return QP_FAILED, 0.0, it      # no descent direction: give up
            # Orient p so that f decreases (dual objective increases).
            sgn = 1.0 if gp < 0.0 else -1.0
            # Longest step t keeping all inequality multipliers >= 0.
            tmax = np.inf
            r = -1
            for i in range(q + 1):
                pi = sgn * (-b2[i] if i < q else 1.0)
                if kind[W[i]] == 0 and pi < -1e-12 * pmax:
                    t = lam[W[i]] / (-pi)
                    if t < tmax:
                        tmax = t
                        r = i
            if r < 0:                          # dual unbounded
                return OUT_INFEASIBLE, np.inf, it
            # Take the step and drop the multiplier that reached zero.
            for i in range(q + 1):
                pi = sgn * (-b2[i] if i < q else 1.0)
                lam[W[i]] += tmax * pi
            lam[W[r]] = 0.0
            inW[W[r]] = False
            w = _remove(W, w, r)
    return QP_FAILED, 0.0, it                  # iteration limit reached


# ---------------------------------------------------------------------------
# One dimension of the main loop of Algorithm 1 (native solver)
#
# Loops over all vertices x, builds B, U, c1 for x once (lines 10-13), then
# solves one dual QP per candidate simplex owned by x (lines 14-20).
# Results are written into the out_* arrays, one entry per candidate.
# ---------------------------------------------------------------------------
@njit(cache=True)
def _process_dimension(P, pw, a1, c_ptr, c_idx, cand_ptr, cand_J, k,
                       want_witness, iter_factor, piv_tol, feas_tol_rel,
                       out_status, out_w, out_iter, out_wit):
    N = P.shape[0]          # number of points
    m = P.shape[1]          # ambient dimension
    J = np.empty(max(k, 1), dtype=np.int64)   # equality indices of one sigma
    for x in range(N):
        # Candidates owned by x (x is their smallest vertex).
        s0 = cand_ptr[x]
        s1 = cand_ptr[x + 1]
        if s0 == s1:
            continue                           # nothing to test for this x
        # Cech neighbours of x: the only points that can matter (line 10).
        nb = c_idx[c_ptr[x]:c_ptr[x + 1]]
        n = nb.shape[0]
        # D = matrix of difference vectors x_i - x (one row per neighbour).
        D = np.empty((n, m))
        for i in range(n):
            for d in range(m):
                D[i, d] = P[nb[i], d] - P[x, d]
        # Gram matrix: the ONLY place where the ambient dimension m enters.
        B = np.dot(D, D.T)                                   # line 11
        U = np.empty(n)
        scale = 1e-300
        for i in range(n):
            U[i] = 0.5 * (pw[nb[i]] - pw[x] - B[i, i])       # line 12
            if B[i, i] > scale:
                scale = B[i, i]                # largest squared distance
        # Upper bound on the dual value: c* > c1  <=>  w(sigma) > a1.
        c1 = 0.5 * (a1 + pw[x])                              # line 13
        feas_tol = feas_tol_rel * scale
        max_iter = iter_factor * n + 50
        # Work arrays for dual_qp, allocated once per x and reused.
        lam = np.zeros(n)
        W = np.empty(n + 1, dtype=np.int64)
        inW = np.zeros(n, dtype=np.bool_)
        kind = np.zeros(n, dtype=np.int64)
        L = np.zeros((n + 1, n + 1))
        b1 = np.empty(n + 1)
        b2 = np.empty(n + 1)
        for c in range(s0, s1):                              # line 14
            # Positions (in x's neighbour list) of the other vertices of sigma.
            for t in range(k):
                J[t] = cand_J[c, t]
            st, cval, its = dual_qp(B, U, J, k, c1, lam, W, inW, kind, L,
                                    b1, b2, max_iter, piv_tol, feas_tol)
            out_status[c] = st
            out_iter[c] = its
            if st == IN_COMPLEX:
                # Filtration value = squared radius (unweighted case).
                out_w[c] = 2.0 * cval - pw[x]                # line 19
                if want_witness:                             # line 20
                    # KKT: y* = x - sum_i lambda_i (x_i - x).
                    for d in range(m):
                        v = P[x, d]
                        for i in range(n):
                            v -= lam[i] * D[i, d]
                        out_wit[c, d] = v
            else:
                out_w[c] = np.inf                  # not in the complex


@njit(cache=True)
def _find_row(S, lo, hi, key, k):
    """Binary search for the row `key` (k vertex ids) among rows lo..hi-1 of S,
    whose rows are sorted lexicographically. Returns its index, or -1."""
    while lo < hi:
        mid = (lo + hi) // 2
        c = 0
        for t in range(k):
            if S[mid, t] != key[t]:
                c = -1 if S[mid, t] < key[t] else 1
                break
        if c == 0:
            return mid
        if c < 0:
            lo = mid + 1
        else:
            hi = mid
    return -1


@njit(cache=True)
def _extend_candidates(N, k, S, S_ptr, a_ptr, a_idx, c_ptr, c_idx):
    """Candidate k-simplices (k >= 2): Lazy_{k-1}(X)_k, i.e. the sets of k+1
    points ALL of whose (k-1)-faces are accepted (k-1)-simplices.

    S     : accepted (k-1)-simplices, one per row (k sorted vertex ids), rows
            sorted lexicographically, so rows owned by vertex x (x = smallest
            vertex) are S[S_ptr[x]:S_ptr[x+1]].
    a_ptr, a_idx : graph of accepted edges (CSR, sorted neighbour lists).
    c_ptr, c_idx : Cech graph (CSR), to locate vertices in B and U.

    Every candidate is a face sigma = (x, y1, .., y_{k-1}) in S extended by a
    vertex z > y_{k-1} adjacent to y_{k-1}; it is kept if the k other faces
    (drop one vertex of sigma, add z) are in S too. Candidates therefore come
    out sorted lexicographically and grouped by owner x.

    Returns cand_ptr (N+1), cand_J (C x k: positions of the k other vertices
    in x's Cech neighbour list), cand (C x (k+1) vertex ids) and faces
    (C x (k+1): faces[c, i] = row of S of the face that omits vertex i)."""
    M = S.shape[0]
    key = np.empty(k, dtype=np.int64)
    fidx = np.empty(k + 1, dtype=np.int64)
    counts = np.zeros(N, dtype=np.int64)
    cand_ptr = np.zeros(N + 1, dtype=np.int64)
    cand_J = np.empty((0, k), dtype=np.int64)
    cand = np.empty((0, k + 1), dtype=np.int64)
    faces = np.empty((0, k + 1), dtype=np.int64)
    for pass_ in range(2):                 # pass 0: count, pass 1: fill
        if pass_ == 1:
            for x in range(N):
                cand_ptr[x + 1] = cand_ptr[x] + counts[x]
            C = cand_ptr[N]
            cand_J = np.empty((C, k), dtype=np.int64)
            cand = np.empty((C, k + 1), dtype=np.int64)
            faces = np.empty((C, k + 1), dtype=np.int64)
        o = 0
        for row in range(M):
            x = S[row, 0]
            last = S[row, k - 1]
            for t in range(a_ptr[last], a_ptr[last + 1]):
                z = a_idx[t]
                if z <= last:
                    continue
                ok = True
                for i in range(k):         # face: drop vertex i of sigma, add z
                    q = 0
                    for u in range(k):
                        if u != i:
                            key[q] = S[row, u]
                            q += 1
                    key[k - 1] = z         # z is the largest: key stays sorted
                    f0 = key[0]
                    j = _find_row(S, S_ptr[f0], S_ptr[f0 + 1], key, k)
                    if j < 0:
                        ok = False
                        break
                    fidx[i] = j
                if not ok:
                    continue
                fidx[k] = row              # the face dropping z is sigma itself
                if pass_ == 0:
                    counts[x] += 1
                    continue
                cn = c_idx[c_ptr[x]:c_ptr[x + 1]]   # x's Cech neighbours
                for u in range(k):
                    cand[o, u] = S[row, u]
                cand[o, k] = z
                for u in range(1, k):
                    cand_J[o, u - 1] = np.searchsorted(cn, S[row, u])
                cand_J[o, k - 1] = np.searchsorted(cn, z)
                for u in range(k + 1):
                    faces[o, u] = fidx[u]
                o += 1
    return cand_ptr, cand_J, cand, faces


# ---------------------------------------------------------------------------
# Cech graph (line 1 of Algorithm 1)
#
# Two points are Cech-neighbours when their balls (radius r, or sqrt(a1+p)
# in the weighted case) overlap, i.e. |x - y| <= 2r. This graph serves two
# purposes (Section 3.3 of the paper):
#   * every alpha edge is a Cech edge, so Cech edges are the edge candidates;
#   * only Cech neighbours of x can rule out a simplex owned by x, so each QP
#     only needs those points as constraints (exactly the same answer).
# ---------------------------------------------------------------------------
def cech_graph(P, pw, a1, kdtree_max_dim=64):
    """1-skeleton of the Cech complex of the weighted ball cover at level a1:
    ball U_x has radius sqrt(a1 + p(x)) (empty if negative); x~y iff the two
    balls meet. Returns (alive mask, indptr, indices) of a symmetric CSR graph
    with sorted neighbour lists."""
    N, m = P.shape
    if not np.isfinite(a1):     # all radii: every pair of balls meets
        # Complete filtration: every pair of points is a neighbour pair.
        alive = np.ones(N, dtype=bool)
        iu, ju = np.triu_indices(N, k=1)       # all pairs i < j
        A = csr_matrix((np.ones(2 * len(iu), dtype=np.int8),
                        (np.concatenate([iu, ju]), np.concatenate([ju, iu]))),
                       shape=(N, N))
        A.sort_indices()
        return alive, A.indptr.astype(np.int64), A.indices.astype(np.int64)
    # Ball radius of each point; a point whose ball is empty (a1 + p < 0) is
    # "dead" and takes no part (only possible in the weighted case).
    R2 = a1 + pw
    alive = R2 >= 0
    R = np.sqrt(np.maximum(R2, 0.0))
    ids = np.flatnonzero(alive)                # numbers of the live points
    Q = P[ids]
    Ra = R[ids]
    slack = 1.0 + 1e-9          # extra candidates are harmless (QP rejects)
    if len(ids) < 2:
        rows = np.empty(0, np.int64)
        cols = np.empty(0, np.int64)
    elif m <= kdtree_max_dim:
        # Low/moderate dimension: KD-tree returns all pairs within distance
        # 2 * (largest radius) directly.
        tree = cKDTree(Q)
        pairs = tree.query_pairs(2.0 * Ra.max() * slack, output_type="ndarray")
        if len(pairs) and np.ptp(Ra) > 0:
            # Weighted case (radii differ): keep pairs with |x-y| <= R_x + R_y.
            d = np.linalg.norm(Q[pairs[:, 0]] - Q[pairs[:, 1]], axis=1)
            pairs = pairs[d <= (Ra[pairs[:, 0]] + Ra[pairs[:, 1]]) * slack]
        rows, cols = pairs[:, 0], pairs[:, 1]
    else:                        # brute force, chunked (high ambient dim)
        # KD-trees degrade in high dimension, so compute all squared distances
        # |a-b|^2 = |a|^2 + |b|^2 - 2 a.b, a block of rows at a time (to keep
        # memory bounded at ~1e6 entries per block).
        sq = np.einsum("ij,ij->i", Q, Q)       # squared norms |q|^2
        nq = len(ids)
        chunk = max(1, int(1e6 // max(nq, 1)))
        rl, cl = [], []
        for s in range(0, nq, chunk):
            e = min(nq, s + chunk)
            d2 = sq[s:e, None] + sq[None, :] - 2.0 * (Q[s:e] @ Q.T)
            thr = ((Ra[s:e, None] + Ra[None, :]) * slack) ** 2
            ii, jj = np.nonzero(d2 <= thr)
            ii = ii + s                        # block-local -> global row
            keep = jj > ii                     # each pair once, no self-pairs
            rl.append(ii[keep])
            cl.append(jj[keep])
        rows = np.concatenate(rl)
        cols = np.concatenate(cl)
    # Back to original point numbers, then build a symmetric CSR adjacency
    # (each pair stored in both directions) with sorted neighbour lists.
    r_orig = ids[rows]
    c_orig = ids[cols]
    A = csr_matrix((np.ones(2 * len(r_orig), dtype=np.int8),
                    (np.concatenate([r_orig, c_orig]),
                     np.concatenate([c_orig, r_orig]))), shape=(N, N))
    A.sort_indices()
    return alive, A.indptr.astype(np.int64), A.indices.astype(np.int64)


def _csr_from_edges(N, E):
    # Symmetric CSR adjacency (sorted neighbour lists) of an edge list E (e x 2).
    A = csr_matrix((np.ones(2 * len(E), dtype=np.int8),
                    (np.concatenate([E[:, 0], E[:, 1]]),
                     np.concatenate([E[:, 1], E[:, 0]]))), shape=(N, N))
    A.sort_indices()
    return A.indptr.astype(np.int64), A.indices.astype(np.int64)


# ---------------------------------------------------------------------------
# DAQP backend
#
# DAQP solves PRIMAL problems   min 1/2 w^T H w + f^T w  s.t.  A w <= b
# by working on their dual. We give it
#     H = identity,  f = 0,  A = G,  b = -U   with  G G^T = B,
# whose dual is   max -1/2 l^T (G G^T) l + U^T l = the paper's dual (10).
# So DAQP solves exactly problem (10), and since G is computed from B alone,
# the solver only ever sees dot products (never the ambient dimension).
# ---------------------------------------------------------------------------
@njit(cache=True)
def _gram_factor(B, tol):
    """Pivoted Cholesky B ~= G G^T with G of size n x k, k = numerical rank
    (at least 1). Uses only the Gram matrix B (dot products)."""
    # Greedy algorithm: at each step pick the row with the largest remaining
    # "unexplained" squared length d[i], use it as a new pivot column, and
    # update the others. Stop when every remaining d[i] is (relatively)
    # zero. For data lying on a k-dimensional affine piece (e.g. a torus
    # embedded in R^1000) this stops after k columns, keeping DAQP small.
    n = B.shape[0]
    d = np.empty(n)
    dmax = 0.0
    for i in range(n):
        d[i] = B[i, i]
        if d[i] > dmax:
            dmax = d[i]
    used = np.zeros(n, dtype=np.bool_)        # rows already used as pivots
    G = np.zeros((n, max(n, 1)))
    k = 0                                      # number of columns built so far
    while k < n:
        # Choose the pivot: largest remaining diagonal.
        p = -1
        best = tol * dmax
        for i in range(n):
            if not used[i] and d[i] > best:
                best = d[i]
                p = i
        if p < 0:
            break                              # rest is numerically zero
        used[p] = True
        piv = np.sqrt(d[p])
        # New column k of G (for the pivot and all not-yet-used rows).
        for i in range(n):
            if used[i] and i != p:
                continue
            s = B[i, p]
            for t in range(k):
                s -= G[i, t] * G[p, t]
            G[i, k] = s / piv
        # Update the remaining squared lengths.
        for i in range(n):
            if not used[i]:
                d[i] -= G[i, k] * G[i, k]
        k += 1
    return np.ascontiguousarray(G[:, :max(k, 1)])


def _process_dimension_daqp(P, pw, a1, c_ptr, c_idx, cand_ptr, cand_J, k,
                            want_witness, out_status, out_w, out_iter, out_wit):
    """Same loop as _process_dimension, with each dual QP solved by DAQP."""
    N = P.shape[0]
    EQ = 5                                   # DAQP: active + immutable
    for x in range(N):
        # Candidates owned by x; skip x if there are none.
        s0, s1 = cand_ptr[x], cand_ptr[x + 1]
        if s0 == s1:
            continue
        nb = c_idx[c_ptr[x]:c_ptr[x + 1]]      # Cech neighbours (line 10)
        n = nb.shape[0]
        c1 = 0.5 * (a1 + pw[x])                               # line 13
        if n == 0:                       # isolated vertex: lambda = 0, c* = 0
            # No constraints at all: the minimum of pi_x is at x itself.
            out_status[s0:s1] = IN_COMPLEX if c1 >= 0 else OUT_BOUND
            out_w[s0:s1] = -pw[x]
            out_iter[s0:s1] = 0
            if want_witness:
                out_wit[s0:s1] = P[x]
            continue
        D = P[nb] - P[x]                       # rows: x_i - x
        B = D @ D.T                                           # line 11
        U = 0.5 * (pw[nb] - pw[x] - np.diag(B))               # line 12
        # Rescale the local problem so that the largest |x_i - x|^2 is 1:
        # DAQP's tolerances are absolute, so this makes them relative.
        # (All dual values scale by 1/s2; the optimal lambda is unchanged.)
        s2 = float(B.diagonal().max())
        if s2 <= 0.0:                    # all neighbours coincide with x
            s2 = 1.0
        G = _gram_factor(B / s2, GRAM_FACTOR_TOL)             # B/s2 = G G^T
        kk = G.shape[1]
        bupper = -U / s2
        # sense[i] = constraint type for DAQP: 0 = inequality, EQ = equality.
        sense = np.zeros(n, dtype=np.int32)
        # One DAQP workspace per vertex x, reused for all its candidates.
        model = daqp.Model()
        model.setup(np.eye(kk), np.zeros(kk), G, bupper,
                    np.full(n, -1e30), sense)          # lower bounds: -inf
        # fval_bound = the paper's early termination: DAQP stops as soon as a
        # dual iterate proves the optimum exceeds c1 (rescaled).
        model.settings = {"fval_bound": min(c1 / s2, 1e30),
                          "primal_tol": DAQP_PRIMAL_TOL}
        update, solve = model.update, model.solve
        Jlist = cand_J[s0:s1, :k].tolist()      # plain ints: cheaper per QP
        st_l, w_l, it_l = [], [], []
        px = float(pw[x])
        for jj, J in enumerate(Jlist):                        # line 14
            # Mark the other vertices of sigma as equality constraints ...
            for j in J:
                sense[j] = EQ
            update(sense=sense)
            _, fval, flag, info = solve()                     # line 16
            # ... and reset them for the next candidate.
            for j in J:
                sense[j] = 0
            it_l.append(info["iterations"])
            cval = fval * s2                   # undo the rescaling
            # flag >= 1: optimum found; -1: bound exceeded OR infeasible.
            if flag >= 1 and cval <= c1:                      # line 17
                st_l.append(IN_COMPLEX)
                w_l.append(2.0 * cval - px)                   # line 19
                if want_witness:                              # line 20
                    out_wit[s0 + jj] = P[x] - info["lam"] @ D
            else:
                st_l.append(OUT_BOUND if flag in (1, 2, -1) else QP_FAILED)
                w_l.append(np.inf)
        # Write this vertex's results into the output arrays in one go.
        out_status[s0:s1] = st_l
        out_w[s0:s1] = w_l
        out_iter[s0:s1] = it_l


def _run_dim(P, pw, a1, c_ptr, c_idx, cand_ptr, cand_J, k, want_witness,
             iter_factor, piv_tol, feas_tol, solver):
    # Allocate the per-candidate output arrays for dimension k, then dispatch
    # to the chosen solver. C = number of candidate k-simplices.
    C = int(cand_ptr[-1])
    status = np.empty(C, dtype=np.int64)       # outcome code per candidate
    wval = np.empty(C)                         # filtration value (or inf)
    iters = np.empty(C, dtype=np.int64)        # active-set iterations
    wit = np.empty((C if want_witness else 0, P.shape[1]))   # witness points
    if C and solver == "daqp":
        _process_dimension_daqp(P, pw, float(a1), c_ptr, c_idx, cand_ptr,
                                np.ascontiguousarray(cand_J, dtype=np.int64), k,
                                want_witness, status, wval, iters, wit)
    elif C:
        _process_dimension(P, pw, float(a1), c_ptr, c_idx, cand_ptr,
                           np.ascontiguousarray(cand_J, dtype=np.int64), k,
                           want_witness, iter_factor, piv_tol, feas_tol,
                           status, wval, iters, wit)
    return status, wval, iters, wit


def _qp_stats(status, iters):
    """'rejected' = ruled out by the early-termination bound or infeasible
    (DAQP does not distinguish the two; the native solver does, see
    'early_terminated' / 'infeasible', which are None for DAQP)."""
    return dict(candidates=int(len(status)),
                accepted=int(np.sum(status == IN_COMPLEX)),
                rejected=int(np.sum((status == OUT_BOUND) | (status == OUT_INFEASIBLE))),
                early_terminated=int(np.sum(status == OUT_BOUND)),
                infeasible=int(np.sum(status == OUT_INFEASIBLE)),
                failed=int(np.sum(status == QP_FAILED)),
                iterations=int(iters.sum()) if len(iters) else 0)


def compute_stats(res):
    """Diagnostics (Cech graph size, QP outcome counts, iterations), computed
    from the raw arrays kept by alpha_complex_dualqp. Call it AFTER stopping
    any timer: it is bookkeeping, not part of building the complex."""
    raw = res["_raw"]
    deg = raw["cech_degree"]                   # number of Cech neighbours
    out = dict(cech_edges=int(deg.sum() // 2),  # each edge counted twice
               max_cech_degree=int(deg.max()) if len(deg) else 0)
    for k in range(res["max_dim"] + 1):
        out[f"dim{k}"] = _qp_stats(raw["status"][k], raw["iters"][k])
        if res["solver"] == "daqp":
            # DAQP cannot tell these two apart (same exit code).
            out[f"dim{k}"]["early_terminated"] = None
            out[f"dim{k}"]["infeasible"] = None
    return out


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
def alpha_complex_dualqp(points, r=None, a1=None, weights=None,
                         compute_witness=False, solver=DEFAULT_SOLVER, max_dim=2,
                         piv_tol=DEFAULT_PIV_TOL, feas_tol=DEFAULT_FEAS_TOL,
                         iter_factor=DEFAULT_ITER_FACTOR):
    """max_dim-skeleton of Alpha(S, p, a1) via Algorithm 1 of Carlsson & Carlsson
    (default max_dim = 2: vertices, edges, triangles).

    Contains no timers and no statistics bookkeeping, so that it can be timed
    as a single block from the input point cloud to the output simplices.

    Parameters
    ----------
    points : (N, m) array
    r : float, radius of the unweighted alpha complex Alpha(S, r) (a1 = r^2)
    a1 : float, maximum power (use instead of r for the weighted complex)
    weights : (N,) array p(x), default 0
    compute_witness : also return the witness map Phi (line 20)
    solver : "daqp" (DAQP, ref. [2] of the paper; default) or "native"
    max_dim : largest simplex dimension d (>= 0) to compute

    r = inf (or a1 = inf) gives the complete alpha filtration (all radii), i.e.
    the whole Delaunay 2-skeleton; the Cech graph is then all pairs of points,
    so cost grows roughly cubically with N.

    Returns a dict with keys
      'ordered'    : the filtration in entry order, a dict with
                     'simplices' (M x (max_dim+1) vertex indices, -1 padded),
                     'dims' (M,) and 'values' (M,) = squared entry radius,
                     sorted by (value, dimension, vertices): every face
                     precedes its cofaces;
      'simplices'  : dict dim k -> (n_k x (k+1)) vertex ids, rows sorted;
      'filtration' : dict dim k -> squared-radius values (same row order);
      'vertices', 'edges', 'triangles' (shortcuts for dims 0, 1, 2; empty
      arrays if max_dim is smaller), 'max_dim',
      'witness' (dict or None), and '_raw' (for compute_stats()).
    """
    # ---- Input checks ------------------------------------------------------
    if solver not in SOLVERS:
        raise ValueError(f"solver must be one of {SOLVERS}")
    if solver == "daqp" and not HAVE_DAQP:
        raise ImportError("solver='daqp' needs the DAQP package: pip install daqp")
    P = np.ascontiguousarray(points, dtype=np.float64)
    N = P.shape[0]
    # The paper works with the power level a1; for the unweighted complex
    # a1 = r^2 (all filtration values are squared radii).
    if a1 is None:
        if r is None:
            raise ValueError("give r (unweighted) or a1")
        a1 = float(r) ** 2
    max_dim = int(max_dim)
    if max_dim < 0:
        raise ValueError("max_dim must be >= 0")
    if not np.isfinite(a1) and solver != "daqp":
        raise ValueError(
            "the complete filtration (r = inf) needs solver='daqp': with all points "
            "as constraints, nearly unbounded empty Voronoi faces on the convex hull "
            "defeat the native solver's infeasibility detection (validated: DAQP "
            "matches GUDHI, native does not)")
    pw = (np.zeros(N) if weights is None
          else np.ascontiguousarray(weights, dtype=np.float64))
    # Per-dimension results: filtration values, witnesses, QP status codes
    # and iteration counts (dicts keyed by dimension 0, 1, 2).
    filt, wit, sts, its = {}, {}, {}, {}

    # Cech graph (line 1)
    alive, c_ptr, c_idx = cech_graph(P, pw, a1)
    deg = np.diff(c_ptr)                       # neighbour count per point

    # k = 0: vertices of G
    # One candidate per live point, owned by itself, with no equalities.
    # (Unweighted: every point is a vertex, entering at radius 0. Weighted:
    # a point's power cell can be empty, so the QP really decides.)
    cand_ptr = np.concatenate([[0], np.cumsum(alive)]).astype(np.int64)
    cand_J = np.zeros((int(cand_ptr[-1]), 1), dtype=np.int64)   # unused (k=0)
    sts[0], wv, its[0], wt = _run_dim(P, pw, a1, c_ptr, c_idx, cand_ptr, cand_J,
                                      0, compute_witness, iter_factor, piv_tol,
                                      feas_tol, solver)
    ok = sts[0] == IN_COMPLEX
    V = np.flatnonzero(alive)[ok]              # accepted vertices
    filt[0] = wv[ok]
    wit[0] = wt[ok] if compute_witness else None

    # Accepted simplices per dimension k: rows of k+1 sorted vertex ids
    # (rows sorted lexicographically), their values, and for k >= 1 the
    # indices of their faces in level k-1 (used by the rounding guard below).
    simp = {0: V[:, None].astype(np.int64)}
    faces_of = {}
    # vpos[v] = row of vertex v among the accepted vertices (-1 if none).
    vpos = np.full(N, -1, dtype=np.int64)
    vpos[V] = np.arange(len(V))

    # k = 1: edges of G, grouped by their smaller vertex
    # Read the Cech graph entry by entry: entry number t is the edge
    # (row[t], c_idx[t]), and pos[t] is the position of c_idx[t] in row[t]'s
    # neighbour list (= its index in B and U). Keep each edge once (x < y).
    if max_dim >= 1:
        row = np.repeat(np.arange(N, dtype=np.int64), deg)
        pos = np.arange(len(c_idx), dtype=np.int64) - c_ptr[row]
        mask = c_idx > row
        ex, ey, epos = row[mask], c_idx[mask], pos[mask]
        # Candidates are already sorted by owner ex; cand_ptr = where each
        # owner's block starts (prefix sums of the number of edges per owner).
        cand_ptr = np.concatenate([[0], np.cumsum(np.bincount(ex, minlength=N))]
                                  ).astype(np.int64)
        sts[1], wv, its[1], wt = _run_dim(P, pw, a1, c_ptr, c_idx, cand_ptr,
                                          epos[:, None], 1, compute_witness,
                                          iter_factor, piv_tol, feas_tol, solver)
        ok = sts[1] == IN_COMPLEX
        E = np.stack([ex[ok], ey[ok]], axis=1)   # accepted edges, sorted (x, y)
        simp[1] = E
        filt[1] = wv[ok]
        wit[1] = wt[ok] if compute_witness else None
        faces_of[1] = np.stack([vpos[E[:, 1]], vpos[E[:, 0]]], axis=1)  # drop x, drop y
        a_ptr, a_idx = _csr_from_edges(N, E)     # graph of accepted edges

    # k >= 2: Lazy_{k-1}(X)_k, built from the accepted (k-1)-simplices
    for k in range(2, max_dim + 1):
        S = simp[k - 1]
        S_ptr = np.concatenate([[0], np.cumsum(np.bincount(S[:, 0], minlength=N))]
                               ).astype(np.int64) if len(S) else np.zeros(N + 1, np.int64)
        cand_ptr, cand_J, cand, faces = _extend_candidates(
            N, k, np.ascontiguousarray(S), S_ptr, a_ptr, a_idx, c_ptr, c_idx)
        sts[k], wv, its[k], wt = _run_dim(P, pw, a1, c_ptr, c_idx, cand_ptr, cand_J,
                                          k, compute_witness, iter_factor, piv_tol,
                                          feas_tol, solver)
        ok = sts[k] == IN_COMPLEX
        simp[k] = cand[ok]
        filt[k] = wv[ok]
        faces_of[k] = faces[ok]
        wit[k] = wt[ok] if compute_witness else None

    # Dimensions that were not reached (no candidates left) are empty.
    for k in range(max_dim + 1):
        if k not in simp:
            simp[k] = np.empty((0, k + 1), dtype=np.int64)
            filt[k] = np.empty(0)
            sts[k] = np.empty(0, dtype=np.int64)
            its[k] = np.empty(0, dtype=np.int64)
            wit[k] = np.empty((0, P.shape[1])) if compute_witness else None

    # Filtration order. A face never enters after its cofaces (V_sigma grows
    # when sigma shrinks); since faces and cofaces come from separate QPs,
    # rounding can break this by ~1e-15, so values are raised to the max over
    # their faces (as GUDHI does), dimension by dimension upwards so that the
    # correction propagates. Then sort by (value, dimension, vertices).
    for k in range(1, max_dim + 1):
        if len(simp[k]):
            filt[k] = np.maximum(filt[k], filt[k - 1][faces_of[k]].max(axis=1))
    width = max_dim + 1
    counts = [len(simp[k]) for k in range(width)]
    table = np.full((sum(counts), width), -1, dtype=np.int64)
    o = 0
    for k in range(width):
        table[o:o + counts[k], :k + 1] = simp[k]
        o += counts[k]
    dims = np.repeat(np.arange(width, dtype=np.int8), counts)
    vals = np.concatenate([filt[k] for k in range(width)])
    # np.lexsort sorts by its LAST key first: value, then dimension (faces
    # before cofaces on ties), then vertex ids (deterministic order).
    keys = tuple(table[:, c] for c in range(width - 1, -1, -1)) + (dims, vals)
    order = np.lexsort(keys)
    ordered = dict(simplices=table[order], dims=dims[order], values=vals[order])

    empty = lambda k: np.empty((0, k + 1), dtype=np.int64)  # noqa: E731
    return dict(vertices=V, edges=simp.get(1, empty(1)),
                triangles=simp[2] if max_dim >= 2 else empty(2),
                simplices=simp, filtration=filt, ordered=ordered, max_dim=max_dim,
                witness=wit if compute_witness else None, a1=a1, solver=solver,
                _raw=dict(cech_degree=deg, status=sts, iters=its))


def warmup(solver=DEFAULT_SOLVER):
    """Trigger Numba compilation (and DAQP loading) so that it is not counted
    in timings."""
    # Run every code path once on a tiny random cloud: both solvers, with and
    # without witnesses, and a 23-dimensional copy so that the brute-force
    # Cech-graph branch and the BLAS routines are initialised too.
    rng = np.random.default_rng(0)
    X = rng.standard_normal((60, 3))
    solvers = SOLVERS if HAVE_DAQP else ("native",)
    for sv in dict.fromkeys((solver,) + solvers):   # unique, chosen one first
        alpha_complex_dualqp(X, r=0.8, compute_witness=True, solver=sv)
        alpha_complex_dualqp(X, r=0.8, compute_witness=False, solver=sv)
        Y = np.hstack([X, np.zeros((60, 20))])   # brute-force Cech path + BLAS
        alpha_complex_dualqp(Y, r=0.8, compute_witness=False, solver=sv)
