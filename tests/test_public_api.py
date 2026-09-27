from __future__ import annotations

import json
import random
from dataclasses import FrozenInstanceError, replace
from typing import Any, cast

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from cloningkit import (
    DesignConfig,
    DesignResult,
    MutagenesisCandidate,
    MutagenesisConfig,
    MutagenesisResult,
    MutationTarget,
    Oligo,
    PairValidationResult,
    PrimerPair,
    ProductVerificationConfig,
    ProductVerificationResult,
    ReactionConditions,
    SpecificityConfig,
    SpecificityResult,
    StructureResult,
    analyze_specificity,
    design_mutagenesis,
    design_primers,
    load_record,
    validate_pair,
    verify_product,
)
from cloningkit.design import design_standard_primers
from cloningkit.io import jsonable, reverse_complement
from cloningkit.mutagenesis import design_inverse_pcr_mutagenesis
from cloningkit.specificity import (
    analyze_template_specificity,
    simulate_product,
)
from cloningkit.thermo import validate_pair as validate_pair_payload

CONDITIONS = ReactionConditions(
    mv_conc=50.0,
    dv_conc=1.5,
    dntp_conc=0.6,
    dna_conc=50.0,
    annealing_temp_c=60.0,
)


def permissive_mutagenesis_config() -> MutagenesisConfig:
    return MutagenesisConfig(
        min_core_length=20,
        max_core_length=20,
        min_tm_c=0.0,
        max_tm_c=100.0,
        target_tm_c=60.0,
        max_delta_tm_c=100.0,
        min_gc_percent=0.0,
        max_gc_percent=100.0,
        max_homopolymer=100,
        max_three_prime_gc=5,
        max_hairpin_risk_kcal_mol=100.0,
        max_homodimer_risk_kcal_mol=100.0,
        max_heterodimer_risk_kcal_mol=100.0,
        max_end_dimer_risk_kcal_mol=100.0,
        codons=("CTG",),
        verify_top=10,
    )


def test_oligo_and_pair_are_normalized_frozen_value_objects() -> None:
    forward = Oligo("ac gu", core_length=3)
    reverse = Oligo("TTAA")
    pair = PrimerPair(forward=forward, reverse=reverse)

    assert forward.sequence_5to3 == "ACGT"
    assert forward.tail_5to3 == "A"
    assert forward.annealing_core_5to3 == "CGT"
    assert forward.resolved_core_length == 3
    assert pair.forward is forward
    assert not hasattr(forward, "__dict__")
    with pytest.raises(FrozenInstanceError):
        forward.sequence_5to3 = "AAAA"  # type: ignore[misc]
    with pytest.raises(ValueError, match="core length"):
        Oligo("ACGT", core_length=5)


def test_typed_validation_stage_supports_nested_access_and_serialization() -> None:
    pair = PrimerPair(
        forward=Oligo("GCACTGCTGTGCGATTTTTATCG", core_length=23),
        reverse=Oligo("CAGACATTCATCGGCCAGAATACCC", core_length=22),
    )
    result = validate_pair(pair, CONDITIONS)

    assert isinstance(result, PairValidationResult)
    assert isinstance(result.forward.hairpin, StructureResult)
    assert result.forward.core_tm_c == pytest.approx(61.5582117682)
    assert result.reverse.tail_5to3 == "CAG"
    assert result.pair.core_delta_tm_c == pytest.approx(1.0709463162)
    assert result.risks.heterodimer_risk_kcal_mol >= 0.0
    payload = cast(dict[str, Any], result.as_dict())
    assert payload == validate_pair_payload(
        pair.forward.sequence_5to3,
        pair.reverse.sequence_5to3,
        forward_core_length=pair.forward.resolved_core_length,
        reverse_core_length=pair.reverse.resolved_core_length,
        conditions=CONDITIONS,
    )
    assert payload["forward"]["full_sequence_5to3"] == pair.forward.sequence_5to3
    assert payload["pair"]["heterodimer"]["dg_kcal_mol"] == pytest.approx(
        result.pair.heterodimer.dg_kcal_mol
    )
    json.dumps(result, default=jsonable, allow_nan=False)


def test_typed_design_stage_returns_composable_primer_pairs() -> None:
    rng = random.Random(1)
    template = SeqRecord(
        Seq("".join(rng.choice("ACGT") for _ in range(1200))),
        id="linear-template",
    )
    config = DesignConfig(
        target_start_1based=500,
        target_length=40,
        product_min_bp=100,
        product_max_bp=300,
        min_gc_percent=20.0,
        max_gc_percent=80.0,
        num_return=2,
    )
    result = design_primers(template, config)

    assert isinstance(result, DesignResult)
    assert result.as_dict() == design_standard_primers(
        template,
        target_start_1based=config.target_start_1based,
        target_length=config.target_length,
        product_min=config.product_min_bp,
        product_max=config.product_max_bp,
        min_size=config.min_size_nt,
        opt_size=config.opt_size_nt,
        max_size=config.max_size_nt,
        min_tm=config.min_tm_c,
        opt_tm=config.opt_tm_c,
        max_tm=config.max_tm_c,
        min_gc=config.min_gc_percent,
        max_gc=config.max_gc_percent,
        max_delta_tm=config.max_delta_tm_c,
        num_return=config.num_return,
        conditions=config.conditions,
    )
    assert 1 <= len(result.pairs) <= 2
    pair = result.pairs[0].primer_pair
    assert isinstance(pair, PrimerPair)
    assert pair.forward.sequence_5to3 == result.pairs[0].forward.sequence_5to3
    payload = cast(dict[str, Any], result.as_dict())
    assert payload["pairs"][0]["rank"] == 0


def test_typed_specificity_and_verification_stages_share_primer_pair(
    circular_tadde_record,
) -> None:
    sequence = str(circular_tadde_record.seq)
    codon_start = 50 + (142 - 1) * 3
    codon_end = codon_start + 3
    expected = sequence[:codon_start] + "CTG" + sequence[codon_end:]
    pair = PrimerPair(
        forward=Oligo(sequence[codon_end : codon_end + 22], core_length=22),
        reverse=Oligo(
            "CAG" + reverse_complement(sequence[codon_start - 22 : codon_start]),
            core_length=22,
        ),
    )

    specificity = analyze_specificity(
        circular_tadde_record, pair, SpecificityConfig(pydna_limit=13)
    )
    verification = verify_product(
        circular_tadde_record,
        pair,
        ProductVerificationConfig(
            pydna_limit=13,
            circularize=True,
            expected_sequence=expected,
        ),
    )

    assert isinstance(specificity, SpecificityResult)
    assert specificity.as_dict() == analyze_template_specificity(
        circular_tadde_record,
        pair.forward.sequence_5to3,
        pair.reverse.sequence_5to3,
        forward_core_length=pair.forward.resolved_core_length,
        reverse_core_length=pair.reverse.resolved_core_length,
        limit=13,
    )
    raw_verification, _ = simulate_product(
        circular_tadde_record,
        pair.forward.sequence_5to3,
        pair.reverse.sequence_5to3,
        forward_core_length=pair.forward.resolved_core_length,
        reverse_core_length=pair.reverse.resolved_core_length,
        limit=13,
        circularize=True,
        expected_sequence=expected,
    )
    assert verification.as_dict() == raw_verification
    assert specificity.pydna.product_count == 1
    assert specificity.forward.exact_core_sites_1based
    assert isinstance(verification, ProductVerificationResult)
    assert verification.product.circular
    assert verification.expected_sequence_match is True
    payload = cast(dict[str, Any], verification.as_dict())
    assert "product" not in payload
    assert payload["primer_roles"] == {
        "forward": "forward",
        "reverse": "reverse",
    }


def test_typed_mutagenesis_stage_returns_candidates(circular_tadde_record) -> None:
    result = design_mutagenesis(
        circular_tadde_record,
        MutationTarget(
            feature_name="TadDE", aa_position_1based=142, target_amino_acid="l"
        ),
        permissive_mutagenesis_config(),
    )

    assert isinstance(result, MutagenesisResult)
    assert result.as_dict() == design_inverse_pcr_mutagenesis(
        circular_tadde_record,
        feature_name="TadDE",
        aa_position=142,
        target_aa="L",
        config=permissive_mutagenesis_config(),
    )
    assert result.target.target_amino_acid == "L"
    assert result.summary.enumerated == 2
    assert result.candidates
    assert isinstance(result.candidates[0], MutagenesisCandidate)
    assert isinstance(result.candidates[0].primer_pair, PrimerPair)
    payload = cast(dict[str, Any], result.as_dict())
    assert payload["target"]["source_coding_codon"] == "GCC"
    unverified = design_mutagenesis(
        circular_tadde_record,
        MutationTarget(
            feature_name="TadDE", aa_position_1based=142, target_amino_acid="L"
        ),
        replace(permissive_mutagenesis_config(), verify_top=0),
    )
    assert all(
        candidate.pydna_verification is None for candidate in unverified.candidates
    )


def test_configs_reject_invalid_values_before_stage_execution() -> None:
    with pytest.raises(ValueError, match="product bounds"):
        DesignConfig(product_min_bp=500, product_max_bp=100)
    with pytest.raises(ValueError, match="pydna limit"):
        SpecificityConfig(pydna_limit=0)
    with pytest.raises(ValueError, match="pydna limit"):
        ProductVerificationConfig(pydna_limit=0)
    with pytest.raises(ValueError, match="one-letter"):
        MutationTarget(
            feature_name="gene", aa_position_1based=2, target_amino_acid="Leu"
        )


def test_load_record_is_part_of_package_root() -> None:
    assert callable(load_record)
