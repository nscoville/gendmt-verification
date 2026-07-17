"""
Finite abstract simplicial complexes, given by a list of maximal faces.

Corresponds to the ambient complex K throughout the paper.
"""

from itertools import combinations


class SimplicialComplex:
    """A finite abstract simplicial complex K, given by its maximal faces.

    Attributes
    ----------
    faces : list of frozenset
        All nonempty faces of K (the "simplices" of the paper), sorted by
        (dimension, vertex set).
    index : dict frozenset -> int
        Maps a face to its position in `faces`.
    n : int
        Number of nonempty faces of K, i.e. |K| in the paper's notation.
    covers : list of (int, int)
        All codimension-one ("covering") pairs (tau_id, sigma_id) with
        tau a facet of sigma in the sense of the (renamed) Definition
        "codimension-one face", i.e. tau \\lessdot sigma.
    """

    def __init__(self, maximal_faces):
        maximal_faces = [frozenset(f) for f in maximal_faces]
        faces = set()
        for mf in maximal_faces:
            items = sorted(mf)
            for r in range(1, len(mf) + 1):
                for c in combinations(items, r):
                    faces.add(frozenset(c))
        self.faces = sorted(faces, key=lambda f: (len(f), sorted(f)))
        self.index = {f: i for i, f in enumerate(self.faces)}
        self.n = len(self.faces)

        self.covers = []
        for sigma in self.faces:
            for v in sigma:
                tau = sigma - {v}
                if tau:  # exclude the empty face
                    self.covers.append((self.index[tau], self.index[sigma]))

    def face_id(self, s):
        return self.index[frozenset(s)]

    def __repr__(self):
        return f"SimplicialComplex(n={self.n} faces)"
