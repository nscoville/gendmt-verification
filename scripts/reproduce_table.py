#!/usr/bin/env python3
"""
Reproduces every numerical claim in Table "tab: data" of the paper
"Generalized gradient vector fields and the complex of discrete Morse
intervals" -- and computes the previously-missing "triangle + two pendant
edges" row (Example "ex: leaves must share") and, as bonus data, the
homology of Mint(solid Delta^3), which the paper's Table left uncomputed.

Usage:
    python scripts/reproduce_table.py            # fast rows only (default)
    python scripts/reproduce_table.py --full      # include the two large
                                                    # tetrahedron examples
                                                    # (several minutes)

Torsion is verified exactly via integer Smith normal form except for the
largest boundary matrices in the tetrahedron examples, where it is skipped
by default for speed (free ranks / Betti numbers there ARE exact -- see
gendmt/homology.py docstring). Pass --torsion to also attempt exact torsion
on the large examples (slow, ~minutes, may still hit the size cap).
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gendmt import SimplicialComplex, analyze

EXAMPLES = [
    ("full triangle", [[0, 1, 2]], False),
    ("triangle + pendant edge", [[0, 1, 2], [2, 3]], False),
    ("triangle + two pendant edges", [[0, 1, 2], [2, 3], [1, 4]], False),
    ("two triangles, shared edge", [[0, 1, 2], [0, 1, 3]], False),
    ("boundary of tetrahedron", [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]], True),
    ("solid tetrahedron", [[0, 1, 2, 3]], True),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                     help="Include the two large (tetrahedron) examples.")
    ap.add_argument("--torsion", action="store_true",
                     help="Attempt exact torsion even on the largest boundary "
                          "matrices (slow; may still be skipped if too large).")
    args = ap.parse_args()

    homology_kwargs = {}
    if args.torsion:
        homology_kwargs["snf_size_limit"] = 2_000_000

    for name, maximal_faces, is_large in EXAMPLES:
        if is_large and not args.full:
            print(f"=== {name} === (skipped; rerun with --full)\n")
            continue
        K = SimplicialComplex(maximal_faces)
        t0 = time.time()
        analyze(K, name=name, compute_homology=True, verbose=True, **homology_kwargs)
        print(f"  [{time.time() - t0:.1f}s]\n")


if __name__ == "__main__":
    main()
