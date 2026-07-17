# gendmt: computational verification for "Generalized gradient vector fields and the complex of discrete Morse intervals"

This repository computes $\mathcal{M}(K)$ (the complex of discrete Morse
matchings) and $\mathcal{M}_{\mathrm{int}}(K)$ (the complex of discrete
Morse intervals) for a finite simplicial complex $K$, by brute-force
enumeration, and their integral reduced homology. It was used to compute
every entry of the paper's data table, including the previously-missing
row for the complex in Example "leaves must share" and, as bonus data,
the homology of $\mathcal{M}_{\mathrm{int}}(\Delta^3)$ (solid tetrahedron),
which the paper leaves uncomputed.

## Correspondence to the paper

| Paper concept | Code |
|---|---|
| Finite simplicial complex $K$ | `SimplicialComplex` (`gendmt/complex.py`) |
| Interval $[\alpha,\beta]$ of the face poset (Definition "interval") | `Interval` (`gendmt/intervals.py`) |
| Non-singular intervals, i.e. vertices of $\mathcal{M}_{\mathrm{int}}(K)$ | `all_nonsingular_intervals` |
| Rank-one intervals, i.e. vertices of $\mathcal{M}(K)$ | `all_nonsingular_intervals` filtered to `rank == 1` |
| Generalized Hasse diagram $\mathcal{H}_W$ (Definition "generalized Hasse") | built internally in `is_acyclic` (`gendmt/hasse.py`) |
| $W$ is a generalized gradient vector field $\iff$ $\mathcal{H}_W$ is acyclic (Theorem "main equivalence") | `is_acyclic` returning `True`/`False` — this single check is the computational core of the whole package |
| Faces of $\mathcal{M}_{\mathrm{int}}(K)$ / $\mathcal{M}(K)$ (Definition "gen morse complex") | `enumerate_facets` (`gendmt/enumerate_faces.py`) |
| Faces are downward closed (Lemma "downward closed") | used to justify computing only facets and taking their downward closure for homology |
| Integral reduced homology of $\mathcal{M}(K)$/$\mathcal{M}_{\mathrm{int}}(K)$ (Table "data") | `reduced_homology` (`gendmt/homology.py`) |

## Installation

```bash
pip install -r requirements.txt
```

## Usage

Reproduce the fast rows of the table:

```bash
python scripts/reproduce_table.py
```

Include the two large tetrahedron examples (several minutes):

```bash
python scripts/reproduce_table.py --full
```

Analyze your own complex:

```python
from gendmt import SimplicialComplex, analyze

K = SimplicialComplex([[0, 1, 2], [2, 3]])  # give maximal faces as vertex lists
analyze(K, name="my complex")
```

## Running the tests

```bash
pip install pytest
pytest tests/
```

The test suite pins down every value already verified against the paper's
table (including a regression test for a facet-enumeration bug described
below), so it should be run after any change to `enumerate_faces.py` or
`homology.py`.

## Known limitations

* **Torsion for the two largest examples is not fully verified.**
  `reduced_homology` computes free ranks (Betti numbers) exactly via modular
  arithmetic (cross-checked at two independent large primes) regardless of
  matrix size, but exact integer torsion requires a Smith normal form
  computation that we found does not scale past a certain matrix size in
  `sympy` (see benchmarks in the git history / PR description). Above that
  size, `reduced_homology` reports the string `"not computed (matrix too
  large for exact SNF in this run)"` in place of the torsion list, rather
  than silently assuming torsion-freeness. This currently affects eight
  nonzero homology groups, all in the $\partial\Delta^3$ and solid-$\Delta^3$
  rows. If you want to push on this further, `snf_size_limit` can be raised
  (see `--torsion` flag on the script), but be aware this can exhaust memory
  for the largest boundary matrices.

* **Facet enumeration is brute-force.** It is fast enough for every example
  in the paper (worst case here: ~209s for $\mathcal{M}_{\mathrm{int}}(\Delta^3)$
  homology) but will not scale indefinitely; there is no attempt at a
  smarter (e.g. shellability-based) enumeration.

* **A bug we hit and fixed, for anyone modifying `enumerate_facets`:** the
  search extends the current candidate set only with strictly
  larger-indexed intervals (to visit each subset exactly once). A naive
  "no larger-indexed extension exists" check is *not* a valid maximality
  test, since a set can be extendable only by a *smaller*-indexed interval
  reachable from a different branch of the search. `enumerate_facets`
  double-checks true maximality against *all* remaining candidates before
  recording a facet; `tests/test_known_values.py::test_facet_maximality_regression`
  guards against this regressing.

## License

Add your preferred license here before publishing.
