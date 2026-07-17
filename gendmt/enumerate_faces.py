"""
Enumeration of the facets of M(K) and Mint(K), and their downward closure.

A face of Mint(K) (Definition "defn: gen morse complex") is a set of
pairwise-disjoint non-singular intervals whose induced generalized discrete
vector field is acyclic (see `hasse.is_acyclic`). M(K) is the sub-poset
spanned by rank-one intervals only (Proposition "prop: subcomplex"). Faces
of M(K)/Mint(K) are automatically downward closed (Lemma "lem: downward
closed"), so it suffices to enumerate the FACETS (maximal faces); the full
face lattice needed for homology is then just the union of the power sets
of the facets.
"""

from itertools import combinations
from .hasse import is_acyclic


def enumerate_facets(K, intervals, rank_filter=None):
    """Enumerate the facets (maximal faces) of the complex whose vertices are
    `intervals` (optionally restricted to those of rank in `rank_filter`,
    e.g. {1} for M(K)) and whose faces are pairwise-disjoint acyclic subsets.

    Returns a list of tuples of interval-indices (indices into `intervals`),
    one per facet.

    Implementation note (a bug we hit and fixed): the search below always
    extends the currently-chosen set with a STRICTLY LARGER interval-index,
    which visits every subset exactly once and is efficient. But that means
    a naive "no further extension found by only trying larger indices"
    check is NOT a correct maximality test: a set can be extendable only by
    an interval of SMALLER index, reachable from a different branch of the
    search that picks that smaller-indexed interval first. So maximality of
    a visited set is verified separately, by checking ALL remaining
    candidates regardless of index -- see the two loops in `rec` below.
    """
    if rank_filter is not None:
        idxs = [i for i, iv in enumerate(intervals) if iv.rank in rank_filter]
    else:
        idxs = list(range(len(intervals)))
    idxs.sort()
    m = len(idxs)
    masks = [intervals[i].mask for i in idxs]

    facets = set()

    def rec(pos, used_mask, chosen_ivs, chosen_positions):
        for j in range(pos, m):
            if used_mask & masks[j]:
                continue
            candidate_iv = intervals[idxs[j]]
            new_chosen = chosen_ivs + [candidate_iv]
            if is_acyclic(K, new_chosen):
                rec(j + 1, used_mask | masks[j], new_chosen, chosen_positions + [idxs[j]])

        if chosen_positions:
            truly_maximal = True
            for j2 in range(m):
                if used_mask & masks[j2]:
                    continue
                candidate_iv2 = intervals[idxs[j2]]
                if is_acyclic(K, chosen_ivs + [candidate_iv2]):
                    truly_maximal = False
                    break
            if truly_maximal:
                facets.add(tuple(sorted(chosen_positions)))

    rec(0, 0, [], [])
    return list(facets)


def downward_closure(facets):
    """facets: list of tuples of interval-indices. Returns dict
    dim -> sorted list of faces (frozensets of interval-indices) of that
    dimension (dimension = |face| - 1, empty face excluded)."""
    all_faces = set()
    for f in facets:
        fs = frozenset(f)
        items = sorted(fs)
        for r in range(1, len(fs) + 1):
            for c in combinations(items, r):
                all_faces.add(frozenset(c))
    by_dim = {}
    for f in all_faces:
        by_dim.setdefault(len(f) - 1, []).append(f)
    for d in by_dim:
        by_dim[d].sort(key=lambda s: sorted(s))
    return by_dim
