"""
Regression tests: assert the table values already known/verified in the
paper, plus the previously-missing row. These pin down correct behavior of
`enumerate_facets` and `reduced_homology` -- in particular they would catch
a reintroduction of the facet-maximality bug described in
`gendmt/enumerate_faces.py`, where a naive "increasing index only" extension
check wrongly reports non-maximal sets as facets (e.g. it originally
reported M(full triangle) as having 24 facets instead of the correct 9).
"""

import pytest

from gendmt import SimplicialComplex, analyze


def betti(homology, dim):
    """Extract the free rank at a given dimension from a homology dict,
    treating a missing entry as rank 0."""
    entry = homology.get(dim)
    return entry[0] if entry is not None else 0


CASES = [
    # (name, maximal_faces, dim_M, facets_M, top_M, betti_M,
    #        dim_Mint, facets_Mint, top_Mint, betti_Mint)
    (
        "full triangle", [[0, 1, 2]],
        2, 9, 9, {1: 4},
        2, 15, 9, {1: 7},
    ),
    (
        "triangle + pendant edge", [[0, 1, 2], [2, 3]],
        3, 14, 12, {2: 2},
        3, 22, 12, {2: 2},
    ),
    (
        "triangle + two pendant edges", [[0, 1, 2], [2, 3], [1, 4]],
        4, 23, 15, {},
        4, 34, 15, {2: 1},
    ),
    (
        "two triangles, shared edge", [[0, 1, 2], [0, 1, 3]],
        4, 54, 32, {2: 2, 3: 5},
        4, 102, 32, {2: 19, 3: 5},
    ),
]


@pytest.mark.parametrize(
    "name,faces,dim_M,facets_M,top_M,betti_M,dim_Mint,facets_Mint,top_Mint,betti_Mint",
    CASES,
    ids=[c[0] for c in CASES],
)
def test_table_row(name, faces, dim_M, facets_M, top_M, betti_M,
                    dim_Mint, facets_Mint, top_Mint, betti_Mint):
    K = SimplicialComplex(faces)
    res = analyze(K, name=name, compute_homology=True, verbose=False)

    assert res["dim_M"] == dim_M
    assert res["num_facets_M"] == facets_M
    assert res["num_top_facets_M"] == top_M
    assert res["dim_Mint"] == dim_Mint
    assert res["num_facets_Mint"] == facets_Mint
    assert res["num_top_facets_Mint"] == top_Mint

    for dim in range(dim_M + 1):
        assert betti(res["homology_M"], dim) == betti_M.get(dim, 0), \
            f"H_{dim}(M({name})) mismatch"
    for dim in range(dim_Mint + 1):
        assert betti(res["homology_Mint"], dim) == betti_Mint.get(dim, 0), \
            f"H_{dim}(Mint({name})) mismatch"

    # torsion-freeness is asserted in the paper for every computed example;
    # check it here wherever it was computed exactly (a list, not the
    # "not computed" string / None placeholder).
    for h in (res["homology_M"], res["homology_Mint"]):
        for dim, (rank, torsion) in h.items():
            if isinstance(torsion, list):
                assert torsion == [], f"unexpected torsion at dim {dim}: {torsion}"


def test_facet_maximality_regression():
    """The specific bug case: a single interval [01,012] pairing an edge
    with the top face of the full triangle must NOT appear as a facet of
    M(triangle) by itself, since it can always be extended (e.g. by pairing
    vertex 0 with edge 02)."""
    K = SimplicialComplex([[0, 1, 2]])
    res = analyze(K, name="regression check", compute_homology=False, verbose=False)
    singleton_facets = [f for f in res["facets_M"] if len(f) == 1]
    assert singleton_facets == [], (
        "Found a size-1 facet of M(triangle); the maximality check has "
        "regressed to only considering later-indexed extensions."
    )
