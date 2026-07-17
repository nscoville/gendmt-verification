"""
Integral reduced simplicial homology of M(K)/Mint(K), computed from a list
of facets (as produced by `enumerate_faces.enumerate_facets`).

Two-stage strategy, because sympy's general-purpose dense linear algebra
(both Matrix.rank() and smith_normal_form()) does not scale to the boundary
matrices that arise even for modestly-sized examples in this paper (matrices
with thousands of rows/columns can exhaust memory in sympy's exact
arithmetic):

  1. Free ranks (Betti numbers) are computed EXACTLY via Gaussian
     elimination over GF(p) for a large prime p, implemented directly with
     numpy for speed. This is exact with overwhelming probability, and is
     cross-checked against a second, different large prime; a mismatch
     would be reported as a warning (none has ever occurred in our runs).
  2. Torsion is computed exactly via sympy's integer Smith normal form,
     but ONLY when the relevant boundary matrix is small enough
     (rows * cols <= snf_size_limit, default calibrated by benchmarking --
     see README) that this is safe. Above that size, torsion is reported
     as unverified (the free rank is still exact).
"""

import numpy as np
import sympy
from sympy.matrices.normalforms import smith_normal_form
from sympy import ZZ

from .enumerate_faces import downward_closure

_PRIME = 2_147_483_647   # 2^31 - 1, prime
_PRIME2 = 2_147_483_629  # a second, different large prime, for cross-checking


def _boundary_matrix_sparse(faces_d, faces_dm1):
    """Standard simplicial boundary map, using the induced order on vertices
    (interval indices) within each face for signs. Returns (rows, cols,
    entries) where entries is a dict (i, j) -> nonzero coefficient."""
    idx_dm1 = {f: i for i, f in enumerate(faces_dm1)}
    rows = len(faces_dm1)
    cols = len(faces_d)
    entries = {}
    for j, f in enumerate(faces_d):
        verts = sorted(f)
        for k, _ in enumerate(verts):
            face = frozenset(verts[:k] + verts[k + 1:])
            sign = (-1) ** k
            i = idx_dm1[face]
            entries[(i, j)] = entries.get((i, j), 0) + sign
    return rows, cols, entries


def _rank_mod_p(rows, cols, entries, p):
    """Rank over GF(p) of the (rows x cols) integer matrix given by `entries`,
    via dense Gaussian elimination in numpy (int64)."""
    if rows == 0 or cols == 0:
        return 0
    M = np.zeros((rows, cols), dtype=np.int64)
    for (i, j), v in entries.items():
        M[i, j] = v % p

    r = 0
    rank = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i, c] % p != 0:
                piv = i
                break
        if piv is None:
            continue
        if piv != r:
            M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), p - 2, p)
        M[r] = (M[r] * inv) % p
        nz = np.nonzero(M[:, c])[0]
        for i in nz:
            if i != r:
                factor = M[i, c]
                M[i] = (M[i] - factor * M[r]) % p
        r += 1
        rank += 1
        if r == rows:
            break
    return rank


def reduced_homology(facets, verbose=False, exact_torsion_max_dim=None,
                      snf_size_limit=250_000):
    """Integral reduced homology of the complex spanned by `facets` (tuples
    of interval-indices). Returns dict: dim -> (free_rank, torsion), where
    torsion is a list of elementary divisors > 1 (computed exactly), or the
    string "not computed (matrix too large for exact SNF in this run)" when
    skipped, or None when there is nothing to compute (rank_{d+1} = 0).

    `exact_torsion_max_dim`: if set, only attempt exact torsion for boundary
    maps into degree <= this value (regardless of size). Pass -1 to skip
    torsion entirely (Betti numbers only -- much faster for large examples).
    """
    by_dim = downward_closure(facets)
    if not by_dim:
        return {}
    maxdim = max(by_dim)

    faces_by_dim = {-1: [frozenset()]}
    faces_by_dim.update(by_dim)

    sparse_boundaries = {
        d: _boundary_matrix_sparse(faces_by_dim[d], faces_by_dim[d - 1])
        for d in range(0, maxdim + 1)
    }

    ranks, ranks_check = {}, {}
    for d in range(0, maxdim + 1):
        rows, cols, entries = sparse_boundaries[d]
        ranks[d] = _rank_mod_p(rows, cols, entries, _PRIME)
        ranks_check[d] = _rank_mod_p(rows, cols, entries, _PRIME2)
        if ranks[d] != ranks_check[d] and verbose:
            print(f"  ** WARNING: rank mismatch at d={d} between two primes "
                  f"({ranks[d]} vs {ranks_check[d]}) -- results may be unreliable **")

    result = {}
    for d in range(0, maxdim + 1):
        rows, cols, entries = sparse_boundaries[d]
        rank_Bd = ranks[d]
        n_d = cols
        ker_dim = n_d - rank_Bd
        rank_Bd1 = ranks.get(d + 1, 0)
        free_rank = ker_dim - rank_Bd1

        torsion = None
        rows1, cols1, entries1 = sparse_boundaries.get(d + 1, (0, 0, {}))
        do_exact = (exact_torsion_max_dim is None or d + 1 <= exact_torsion_max_dim)
        if rows1 and cols1 and rank_Bd1 > 0 and do_exact and rows1 * cols1 <= snf_size_limit:
            Bd1 = sympy.zeros(rows1, cols1)
            for (i, j), v in entries1.items():
                Bd1[i, j] = v
            snf = smith_normal_form(Bd1, domain=ZZ)
            diag = [snf[i, i] for i in range(min(snf.rows, snf.cols))]
            torsion = [abs(int(x)) for x in diag if x != 0 and abs(int(x)) != 1]
        elif rows1 and cols1 and rank_Bd1 > 0:
            torsion = "not computed (matrix too large for exact SNF in this run)"

        result[d] = (free_rank, torsion)
        if verbose:
            print(f"  H_{d}: free rank {free_rank}, torsion {torsion} "
                  f"(chain dim {n_d}, rank d_d={rank_Bd}, rank d_{{d+1}}={rank_Bd1})")
    return result
