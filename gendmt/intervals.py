"""
Intervals of the face poset of K.

Corresponds to Definition "defn: interval": an interval [alpha, beta] with
alpha subseteq beta both faces of K, consisting of every gamma with
alpha <= gamma <= beta. Rank = dim(beta) - dim(alpha) (equivalently
|beta| - |alpha|). Non-singular means alpha proper subset beta (rank >= 1).
"""

from itertools import combinations


class Interval:
    """An interval [alpha, beta] of the face poset of K.

    `mask` is a bitmask (python int) over face ids: bit i is set iff the
    face with that id lies in [alpha, beta]. This makes disjointness checks
    (used throughout `enumerate_faces.py`) a single AND operation.
    """

    __slots__ = ("alpha", "beta", "mask", "rank")

    def __init__(self, alpha, beta, mask, rank):
        self.alpha = alpha
        self.beta = beta
        self.mask = mask
        self.rank = rank

    def __repr__(self):
        return f"[{sorted(self.alpha)},{sorted(self.beta)}]"


def all_nonsingular_intervals(K):
    """Every non-singular interval [alpha, beta] of K's face poset
    (Definition "defn: interval", restricted to alpha proper subset beta,
    matching the vertex set of Mint(K) in Definition "defn: gen morse complex").
    """
    intervals = []
    for alpha in K.faces:
        for beta in K.faces:
            if alpha < beta:  # proper subset ==> non-singular
                diff = beta - alpha
                mask = 0
                for r in range(0, len(diff) + 1):
                    for c in combinations(sorted(diff), r):
                        gamma = alpha | frozenset(c)
                        mask |= (1 << K.index[gamma])
                rank = len(beta) - len(alpha)
                intervals.append(Interval(alpha, beta, mask, rank))
    return intervals
