from __future__ import annotations

import random

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from cloningkit.mutagenesis import (
    MutagenesisConfig,
    build_mutant_record,
    design_inverse_pcr_mutagenesis,
    find_feature,
)


def permissive_config(
    *, min_core_length: int = 20, max_core_length: int = 20
) -> MutagenesisConfig:
    return MutagenesisConfig(
        min_core_length=min_core_length,
        max_core_length=max_core_length,
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


def test_exhaustive_adjacent_search_and_pydna_verification(
    circular_tadde_record,
) -> None:
    result = design_inverse_pcr_mutagenesis(
        circular_tadde_record,
        feature_name="TadDE",
        aa_position=142,
        target_aa="L",
        config=permissive_config(),
    )
    assert result["target"]["source_coding_codon"] == "GCC"
    assert result["target"]["candidate_coding_codons"] == ["CTG"]
    assert result["summary"]["enumerated"] == 2
    assert result["summary"]["hard_filter_passed"] == 2
    assert result["summary"]["pareto_count"] >= 1
    verified = [
        candidate
        for candidate in result["candidates"]
        if candidate["pydna_verification"] is not None
    ]
    assert verified
    assert all(
        candidate["pydna_verification"]["expected_sequence_match"]
        for candidate in verified
    )
    assert result["interpretation"]["global_scope"].startswith("exhaustive only")


def test_origin_crossing_reverse_core_is_supported(origin_target_record) -> None:
    result = design_inverse_pcr_mutagenesis(
        origin_target_record,
        feature_name="origin-gene",
        aa_position=2,
        target_aa="L",
        config=permissive_config(),
    )
    assert result["target"]["genomic_codon_start_1based"] == 4
    assert result["summary"]["pydna_verified_count"] >= 1


def test_reverse_strand_target_uses_coding_orientation(reverse_tadde_record) -> None:
    result = design_inverse_pcr_mutagenesis(
        reverse_tadde_record,
        feature_name="reverse-TadDE",
        aa_position=142,
        target_aa="L",
        config=permissive_config(),
    )
    assert result["target"]["feature_strand"] == -1
    assert result["target"]["source_coding_codon"] == "GCC"
    assert result["target"]["genomic_codon_start_1based"] == 75
    assert result["summary"]["pydna_verified_count"] >= 1
    feature = find_feature(reverse_tadde_record, "reverse-TadDE")
    mutant = build_mutant_record(reverse_tadde_record, feature, 142, "CTG")
    assert str(feature.extract(mutant.seq).translate())[141] == "L"


def test_duplicate_search_axes_are_rejected(circular_tadde_record) -> None:
    from dataclasses import replace

    with pytest.raises(ValueError, match="tail sides must not contain duplicates"):
        design_inverse_pcr_mutagenesis(
            circular_tadde_record,
            feature_name="TadDE",
            aa_position=142,
            target_aa="L",
            config=replace(permissive_config(), tail_sides=("forward", "forward")),
        )
    with pytest.raises(ValueError, match="codons must not contain duplicates"):
        design_inverse_pcr_mutagenesis(
            circular_tadde_record,
            feature_name="TadDE",
            aa_position=142,
            target_aa="L",
            config=replace(permissive_config(), codons=("CTG", "ctg")),
        )


def test_codon_start_phase_is_applied() -> None:
    rng = random.Random(808)
    sequence = "".join(rng.choice("ACGT") for _ in range(400))
    sequence = sequence[:51] + "GCC" + sequence[54:]
    record = SeqRecord(Seq(sequence), id="phased-plasmid")
    record.annotations = {"molecule_type": "DNA", "topology": "circular"}
    record.features = [
        SeqFeature(
            SimpleLocation(50, 351, strand=1),
            type="CDS",
            qualifiers={
                "label": ["phased-gene"],
                "codon_start": ["2"],
                "transl_table": ["1"],
            },
        )
    ]
    result = design_inverse_pcr_mutagenesis(
        record,
        feature_name="phased-gene",
        aa_position=1,
        target_aa="L",
        config=permissive_config(),
    )
    assert result["target"]["feature_codon_start"] == 2
    assert result["target"]["feature_translation_table"] == 1
    assert result["target"]["source_coding_codon"] == "GCC"
    assert result["target"]["genomic_codon_start_1based"] == 52
    assert result["summary"]["pydna_verified_count"] >= 1


def test_reverse_strand_codon_start_phase_is_applied() -> None:
    rng = random.Random(810)
    sequence = "".join(rng.choice("ACGT") for _ in range(500))
    extracted = "AATGGCC" + "".join(rng.choice("ACGT") for _ in range(296))
    feature_start, feature_end = 50, 353
    genomic_feature = str(Seq(extracted).reverse_complement())
    sequence = sequence[:feature_start] + genomic_feature + sequence[feature_end:]
    record = SeqRecord(Seq(sequence), id="reverse-phased-plasmid")
    record.annotations = {"molecule_type": "DNA", "topology": "circular"}
    record.features = [
        SeqFeature(
            SimpleLocation(feature_start, feature_end, strand=-1),
            type="CDS",
            qualifiers={
                "label": ["reverse-phased-gene"],
                "codon_start": ["2"],
                "transl_table": ["1"],
            },
        )
    ]
    result = design_inverse_pcr_mutagenesis(
        record,
        feature_name="reverse-phased-gene",
        aa_position=2,
        target_aa="L",
        config=permissive_config(),
    )
    assert result["target"]["feature_strand"] == -1
    assert result["target"]["feature_codon_start"] == 2
    assert result["target"]["source_coding_codon"] == "GCC"
    assert result["target"]["genomic_codon_start_1based"] == 347
    assert result["summary"]["pydna_verified_count"] >= 1


def test_feature_translation_table_controls_allowed_target_codons() -> None:
    from dataclasses import replace

    rng = random.Random(809)
    sequence = "".join(rng.choice("ACGT") for _ in range(500))
    sequence = sequence[:90] + "ATGGCC" + sequence[96:]
    record = SeqRecord(Seq(sequence), id="table4-plasmid")
    record.annotations = {"molecule_type": "DNA", "topology": "circular"}
    record.features = [
        SeqFeature(
            SimpleLocation(90, 393, strand=1),
            type="CDS",
            qualifiers={
                "label": ["table4-gene"],
                "transl_table": ["4"],
            },
        )
    ]
    config = replace(permissive_config(), codons=("TGA",), verify_top=0)
    result = design_inverse_pcr_mutagenesis(
        record,
        feature_name="table4-gene",
        aa_position=2,
        target_aa="W",
        config=config,
    )
    assert result["target"]["feature_translation_table"] == 4
    assert result["target"]["candidate_coding_codons"] == ["TGA"]
    assert result["summary"]["enumerated"] == 2
    feature = find_feature(record, "table4-gene")
    mutant = build_mutant_record(record, feature, 2, "TGA")
    assert str(feature.extract(mutant.seq).translate(table=4))[1] == "W"


def test_build_mutant_record_changes_exactly_one_codon(circular_tadde_record) -> None:
    feature = find_feature(circular_tadde_record, "TadDE")
    mutant = build_mutant_record(circular_tadde_record, feature, 142, "CTG")
    differences = [
        index
        for index, (wild_type, changed) in enumerate(
            zip(circular_tadde_record.seq, mutant.seq, strict=True)
        )
        if wild_type != changed
    ]
    start = 50 + (142 - 1) * 3
    assert differences == [start, start + 1, start + 2]
    assert str(feature.extract(mutant.seq).translate())[141] == "L"


def test_mutagenesis_rejects_out_of_scope_or_invalid_annotation(
    circular_tadde_record,
) -> None:
    from copy import deepcopy
    from dataclasses import replace

    with pytest.raises(ValueError, match="residue-1/start-codon"):
        design_inverse_pcr_mutagenesis(
            circular_tadde_record,
            feature_name="TadDE",
            aa_position=1,
            target_aa="L",
            config=permissive_config(),
        )

    with pytest.raises(ValueError, match="must differ from the source"):
        design_inverse_pcr_mutagenesis(
            circular_tadde_record,
            feature_name="TadDE",
            aa_position=142,
            target_aa="A",
            config=permissive_config(),
        )

    invalid_table = deepcopy(circular_tadde_record)
    invalid_table.features[0].qualifiers["transl_table"] = ["999"]
    with pytest.raises(ValueError, match="supported DNA table"):
        design_inverse_pcr_mutagenesis(
            invalid_table,
            feature_name="TadDE",
            aa_position=142,
            target_aa="L",
            config=permissive_config(),
        )

    with pytest.raises(ValueError, match="annealing cores overlap"):
        design_inverse_pcr_mutagenesis(
            circular_tadde_record,
            feature_name="TadDE",
            aa_position=142,
            target_aa="L",
            config=replace(permissive_config(), max_core_length=400),
        )


def test_pbh001_a142l_primer_context_regression() -> None:
    rng = random.Random(1420)
    sequence = "".join(rng.choice("ACGT") for _ in range(700))
    feature_start = 50
    codon_start = feature_start + (142 - 1) * 3
    forward_core = "GCACTGCTGTGCGATTTTTATCG"
    reverse_core = "ACATTCATCGGCCAGAATACCC"
    upstream = str(Seq(reverse_core).reverse_complement())
    sequence = (
        sequence[: codon_start - len(upstream)]
        + upstream
        + "GCC"
        + forward_core
        + sequence[codon_start + 3 + len(forward_core) :]
    )
    record = SeqRecord(Seq(sequence), id="synthetic-pBH001-context")
    record.annotations = {"molecule_type": "DNA", "topology": "circular"}
    record.features = [
        SeqFeature(
            SimpleLocation(feature_start, 550, strand=1),
            type="misc_feature",
            qualifiers={"label": ["TadDE"]},
        )
    ]

    result = design_inverse_pcr_mutagenesis(
        record,
        feature_name="TadDE",
        aa_position=142,
        target_aa="L",
        config=permissive_config(min_core_length=22, max_core_length=23),
    )
    assert result["target"]["source_coding_codon"] == "GCC"
    assert result["target"]["feature_type"] == "misc_feature"
    candidate = next(
        candidate
        for candidate in result["candidates"]
        if candidate["id"] == "CTG-reverse-F23-R22"
    )
    assert candidate["forward"]["full_sequence_5to3"] == forward_core
    assert candidate["reverse"]["full_sequence_5to3"] == "CAG" + reverse_core
    assert candidate["pydna_verification"]["expected_sequence_match"] is True
