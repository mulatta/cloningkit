from __future__ import annotations

import pytest

from cloningkit.conditions import ReactionConditions
from cloningkit.thermo import split_core, structural_risks, validate_pair

CONDITIONS = ReactionConditions(
    mv_conc=50.0,
    dv_conc=1.5,
    dntp_conc=0.6,
    dna_conc=50.0,
    annealing_temp_c=60.0,
)


def test_split_core_treats_only_five_prime_prefix_as_tail() -> None:
    assert split_core("CAGACATTCATCGGCCAGAATACCC", 22) == (
        "CAG",
        "ACATTCATCGGCCAGAATACCC",
    )
    with pytest.raises(ValueError, match="core length"):
        split_core("ACGT", 5)


def test_proposed_a142l_pair_reproduces_primer3_metrics() -> None:
    result = validate_pair(
        "GCACTGCTGTGCGATTTTTATCG",
        "CAGACATTCATCGGCCAGAATACCC",
        forward_core_length=23,
        reverse_core_length=22,
        conditions=CONDITIONS,
    )
    assert result["forward"]["core_tm_c"] == pytest.approx(61.5582117682)
    assert result["reverse"]["core_tm_c"] == pytest.approx(60.4872654520)
    assert result["pair"]["core_delta_tm_c"] == pytest.approx(1.0709463162)
    assert result["reverse"]["tail_5to3"] == "CAG"
    assert result["reverse"]["full_length_nt"] == 25
    assert (
        result["pair"]["end_stability_forward_3prime_against_reverse"][
            "ascii_structure"
        ]
        is None
    )
    risks = structural_risks(result)
    assert set(risks) == {
        "hairpin_risk_kcal_mol",
        "homodimer_risk_kcal_mol",
        "heterodimer_risk_kcal_mol",
        "end_dimer_risk_kcal_mol",
    }


def test_jlo_core_tm_is_distinct_from_full_tailed_oligo_tm() -> None:
    result = validate_pair(
        "CTGGCACTGCTGTGCGATTTTTA",
        "ACATTCATCGGCCAGAATAC",
        forward_core_length=20,
        reverse_core_length=20,
        conditions=CONDITIONS,
    )
    assert result["forward"]["tail_5to3"] == "CTG"
    assert result["forward"]["core_tm_c"] == pytest.approx(57.6697, abs=1e-3)
    assert result["reverse"]["core_tm_c"] == pytest.approx(55.6041, abs=1e-3)
    assert result["pair"]["core_delta_tm_c"] == pytest.approx(2.0656, abs=1e-3)


def test_conditions_and_full_oligo_limits_are_validated() -> None:
    with pytest.raises(ValueError, match="finite and nonnegative"):
        ReactionConditions(mv_conc=float("nan"))
    with pytest.raises(ValueError, match="between 0 and 100"):
        ReactionConditions(annealing_temp_c=-1.0)
    with pytest.raises(ValueError, match="at most 60 nt"):
        validate_pair("A" * 61, "ACGT" * 5, conditions=CONDITIONS)
