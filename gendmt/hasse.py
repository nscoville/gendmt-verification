"""
The generalized Hasse diagram and its acyclicity test.

Implements Definition "defn: generalized Hasse": for a chosen collection of
pairwise-disjoint non-singular intervals (with every remaining simplex of K
regarded as its own singular interval), draw an arrow tau -> sigma for every
codimension-one pair tau < sigma of K:
  - upward (tau -> sigma) if tau and sigma lie in the SAME chosen interval;
  - downward (sigma -> tau) otherwise.

By Theorem "thm: main equivalence" (the paper's central result), the induced
generalized discrete vector field is a generalized gradient vector field if
and only if this directed graph has no directed cycle. That is the
computational criterion used everywhere in this package: membership of a set
of intervals as a face of Mint(K) (or M(K), restricting to rank-1 intervals)
is decided entirely by this acyclicity check.
"""


def is_acyclic(K, chosen_intervals):
    """True iff the generalized Hasse diagram of K, with `chosen_intervals`
    (a list of pairwise-disjoint Interval objects) as the non-singular
    intervals and every other simplex singular, has no directed cycle.

    Cost is linear in |K.covers| (i.e. in the fixed complex K), independent
    of which intervals are chosen -- this is what makes repeated calls
    inside the search of `enumerate_faces.py` cheap.
    """
    owner = [-1] * K.n
    for idx, iv in enumerate(chosen_intervals):
        mm = iv.mask
        i = 0
        while mm:
            if mm & 1:
                owner[i] = idx
            mm >>= 1
            i += 1

    succ = [[] for _ in range(K.n)]
    for tau, sigma in K.covers:
        if owner[tau] != -1 and owner[tau] == owner[sigma]:
            succ[tau].append(sigma)   # upward arrow: same interval
        else:
            succ[sigma].append(tau)   # downward arrow: different intervals

    WHITE, GRAY, BLACK = 0, 1, 2
    color = [WHITE] * K.n

    def dfs(u):
        color[u] = GRAY
        for w in succ[u]:
            if color[w] == GRAY:
                return True
            if color[w] == WHITE and dfs(w):
                return True
        color[u] = BLACK
        return False

    for u in range(K.n):
        if color[u] == WHITE and dfs(u):
            return False
    return True
