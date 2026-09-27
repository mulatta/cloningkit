from __future__ import annotations

import pytest

from cloningkit.io import circular_occurrences, normalize_dna, rotation_equivalent
from cloningkit.pareto import dominates, pareto_front_indices


def test_circular_occurrence_can_cross_origin() -> None:
    assert circular_occurrences("AACCGGTT", "TTAAC", True) == [6]
    assert circular_occurrences("AACCGGTT", "TTAAC", False) == []


def test_normalization_and_rotation() -> None:
    assert normalize_dna("ac gu") == "ACGT"
    assert rotation_equivalent("AACCGG", "GGAACC")
    assert not rotation_equivalent("AACCGG", "AACCGT")
    with pytest.raises(ValueError, match="only A/C/G/T"):
        normalize_dna("ACGN")


def test_pareto_minimization() -> None:
    vectors = [(1.0, 3.0), (2.0, 2.0), (3.0, 3.0), (1.0, 4.0)]
    assert dominates(vectors[0], vectors[3])
    assert pareto_front_indices(vectors) == [0, 1]
