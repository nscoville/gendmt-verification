"""Top-level convenience function tying the pipeline together for a single K."""

from .intervals import all_nonsingular_intervals
from .enumerate_faces import enumerate_facets
from .homology import reduced_homology


def analyze(K, name="K", compute_homology=True, verbose=True, **homology_kwargs):
    """Compute M(K) and Mint(K) facet/dimension data for K, and optionally
    their integral reduced homology.

    Returns a dict with keys:
      intervals, facets_M, facets_Mint,
      dim_M, dim_Mint, num_facets_M, num_facets_Mint,
      num_top_facets_M, num_top_facets_Mint,
      and (if compute_homology) homology_M, homology_Mint.
    """
    intervals = all_nonsingular_intervals(K)
    rank1 = [i for i, iv in enumerate(intervals) if iv.rank == 1]

    if verbose:
        print(f"=== {name} ===")
        print(f"  |K| (nonempty faces) = {K.n}")
        print(f"  # non-singular intervals total (Mint vertices) = {len(intervals)}")
        print(f"  # rank-1 intervals (M vertices) = {len(rank1)}")

    facets_M = enumerate_facets(K, intervals, rank_filter={1})
    facets_Mint = enumerate_facets(K, intervals, rank_filter=None)

    dim_M = max(len(f) for f in facets_M) - 1
    dim_Mint = max(len(f) for f in facets_Mint) - 1
    top_M = [f for f in facets_M if len(f) - 1 == dim_M]
    top_Mint = [f for f in facets_Mint if len(f) - 1 == dim_Mint]

    out = {
        "intervals": intervals,
        "facets_M": facets_M,
        "facets_Mint": facets_Mint,
        "dim_M": dim_M,
        "dim_Mint": dim_Mint,
        "num_facets_M": len(facets_M),
        "num_facets_Mint": len(facets_Mint),
        "num_top_facets_M": len(top_M),
        "num_top_facets_Mint": len(top_Mint),
    }

    if verbose:
        print(f"  M(K):    dim={dim_M}, #facets={len(facets_M)}, #top facets={len(top_M)}")
        print(f"  Mint(K): dim={dim_Mint}, #facets={len(facets_Mint)}, #top facets={len(top_Mint)}")

    if compute_homology:
        if verbose:
            print("  Homology of M(K):")
        hM = reduced_homology(facets_M, verbose=verbose, **homology_kwargs)
        if verbose:
            print("  Homology of Mint(K):")
        hMint = reduced_homology(facets_Mint, verbose=verbose, **homology_kwargs)
        out["homology_M"] = hM
        out["homology_Mint"] = hMint

    return out
