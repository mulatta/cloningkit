# Copyright (c) 2026 Seungwon Lee
# SPDX-License-Identifier: MIT
"""Exhaustive adjacent, non-overlapping inverse-PCR mutagenesis design."""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field
from typing import Any, Literal

import python_codon_tables
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation

from .conditions import ReactionConditions
from .io import circular_occurrences, is_circular, normalize_dna, reverse_complement
from .pareto import pareto_front_indices
from .provenance import tool_versions
from .serialization import JsonModel
from .specificity import simulate_product
from .thermo import structural_risks, validate_pair

type MutationTailSide = Literal["forward", "reverse"]


@dataclass(frozen=True, slots=True)
class MutagenesisConfig(JsonModel):
    """Search bounds and hard filters for adjacent inverse-PCR candidates."""

    min_core_length: int = 18
    max_core_length: int = 30
    min_tm_c: float = 55.0
    max_tm_c: float = 65.0
    target_tm_c: float = 60.0
    max_delta_tm_c: float = 2.5
    min_gc_percent: float = 35.0
    max_gc_percent: float = 65.0
    max_homopolymer: int = 5
    max_three_prime_gc: int = 3
    max_hairpin_risk_kcal_mol: float = 3.0
    max_homodimer_risk_kcal_mol: float = 9.0
    max_heterodimer_risk_kcal_mol: float = 9.0
    max_end_dimer_risk_kcal_mol: float = 5.0
    codon_table: str = "e_coli_316407"
    codons: tuple[str, ...] = ()
    tail_sides: tuple[MutationTailSide, ...] = ("forward", "reverse")
    pydna_limit: int = 13
    verify_top: int = 25
    include_failed: bool = False
    conditions: ReactionConditions = field(default_factory=ReactionConditions)


def _feature_names(feature: Any) -> set[str]:
    names: set[str] = set()
    for key in ("label", "gene", "product", "locus_tag"):
        names.update(str(value) for value in feature.qualifiers.get(key, []))
    return names


def find_feature(record: Any, name: str) -> Any:
    matches = [
        feature
        for feature in record.features
        if name.casefold() in {value.casefold() for value in _feature_names(feature)}
    ]
    if len(matches) != 1:
        raise ValueError(f"expected one feature named {name!r}, found {len(matches)}")
    feature = matches[0]
    if feature.location is None:
        raise ValueError("target feature must have a sequence location")
    if isinstance(feature.location, CompoundLocation):
        raise NotImplementedError(
            "compound/origin-spanning target features are not supported"
        )
    if feature.location.strand not in (-1, 1):
        raise ValueError("target feature must have strand +1 or -1")
    return feature


def _circular_slice(sequence: str, start: int, length: int) -> str:
    size = len(sequence)
    return "".join(sequence[(start + offset) % size] for offset in range(length))


def _target_codons(
    target_aa: str, config: MutagenesisConfig, translation_table: int
) -> list[tuple[str, float]]:
    target_aa = target_aa.upper()
    if len(target_aa) != 1:
        raise ValueError("target amino acid must be a one-letter code")
    try:
        genetic_code = CodonTable.unambiguous_dna_by_id[translation_table]
    except KeyError as error:
        raise ValueError(
            f"unsupported feature translation table {translation_table}"
        ) from error
    all_codons = sorted(
        codon
        for codon, amino_acid in genetic_code.forward_table.items()
        if amino_acid == target_aa
    )
    if not all_codons:
        raise ValueError(f"unsupported target amino acid {target_aa!r}")
    usage = python_codon_tables.get_codons_table(config.codon_table)
    frequencies = usage.get(target_aa, {})
    selected = (
        [codon.upper() for codon in config.codons] if config.codons else all_codons
    )
    invalid = sorted(set(selected) - set(all_codons))
    if invalid:
        raise ValueError(f"codons do not encode {target_aa}: {invalid}")
    return [(codon, float(frequencies.get(codon, 0.0))) for codon in selected]


def _target_context(record: Any, feature: Any, aa_position: int) -> dict[str, Any]:
    if aa_position < 1:
        raise ValueError("amino-acid position is 1-based and must be positive")
    extracted_sequence = normalize_dna(
        str(feature.extract(record.seq)), name="target CDS"
    )
    try:
        codon_start = int(feature.qualifiers.get("codon_start", [1])[0])
    except (IndexError, TypeError, ValueError) as error:
        raise ValueError("feature codon_start must be 1, 2, or 3") from error
    if codon_start not in (1, 2, 3):
        raise ValueError("feature codon_start must be 1, 2, or 3")
    if aa_position == 1 and codon_start == 1:
        raise ValueError(
            "residue-1/start-codon mutagenesis is outside this workflow because "
            "CDS initiator semantics differ from ordinary codon translation"
        )
    phase = codon_start - 1
    coding_sequence = extracted_sequence[phase:]
    coding_offset = (aa_position - 1) * 3
    feature_offset = phase + coding_offset
    if coding_offset + 3 > len(coding_sequence):
        raise ValueError("amino-acid position lies outside the target feature")
    codon = coding_sequence[coding_offset : coding_offset + 3]
    try:
        table_id = int(feature.qualifiers.get("transl_table", [1])[0])
        genetic_code = CodonTable.unambiguous_dna_by_id[table_id]
    except (IndexError, KeyError, TypeError, ValueError) as error:
        raise ValueError(
            "feature transl_table must name a supported DNA table"
        ) from error
    amino_acid = genetic_code.forward_table.get(codon)
    if amino_acid is None:
        amino_acid = "*" if codon in genetic_code.stop_codons else "X"
    start = int(feature.location.start)
    end = int(feature.location.end)
    if feature.location.strand == 1:
        genomic_start = start + feature_offset
        genomic_end = genomic_start + 3
    else:
        genomic_end = end - feature_offset
        genomic_start = genomic_end - 3
    return {
        "coding_sequence": coding_sequence,
        "coding_codon": codon,
        "amino_acid": amino_acid,
        "genomic_start_0based": genomic_start,
        "genomic_end_0based_exclusive": genomic_end,
        "strand": feature.location.strand,
        "translation_table": table_id,
        "codon_start": codon_start,
        "phase_bases": phase,
    }


def build_mutant_record(
    record: Any,
    feature: Any,
    aa_position: int,
    coding_codon: str,
) -> Any:
    """Return a deep-copied record with exactly one coding codon replaced."""
    context = _target_context(record, feature, aa_position)
    genomic_codon = (
        coding_codon if context["strand"] == 1 else reverse_complement(coding_codon)
    )
    start = context["genomic_start_0based"]
    end = context["genomic_end_0based_exclusive"]
    sequence = str(record.seq).upper()
    target_amino_acid = str(
        Seq(coding_codon).translate(table=context["translation_table"])
    )
    mutant = copy.deepcopy(record)
    mutant.seq = Seq(sequence[:start] + genomic_codon + sequence[end:])
    mutant.id = f"{record.id}_{context['amino_acid']}{aa_position}{target_amino_acid}"
    mutant.description = f"{record.description}; in-silico mutagenesis design"
    mutant.features.append(
        SeqFeature(
            SimpleLocation(start, end, strand=context["strand"]),
            type="misc_difference",
            qualifiers={
                "label": [f"{context['amino_acid']}{aa_position}{target_amino_acid}"],
                "note": [
                    f"designed coding codon {context['coding_codon']}->{coding_codon}; experimental validation required"
                ],
            },
        )
    )
    return mutant


def _compact_validation(validation: dict[str, Any]) -> dict[str, Any]:
    def compact_oligo(oligo: dict[str, Any]) -> dict[str, Any]:
        return {
            key: value
            for key, value in oligo.items()
            if key not in {"hairpin", "homodimer"}
        } | {
            "hairpin_tm_c": oligo["hairpin"]["tm_c"],
            "hairpin_dg_kcal_mol": oligo["hairpin"]["dg_kcal_mol"],
            "homodimer_tm_c": oligo["homodimer"]["tm_c"],
            "homodimer_dg_kcal_mol": oligo["homodimer"]["dg_kcal_mol"],
        }

    pair = validation["pair"]
    return {
        "forward": compact_oligo(validation["forward"]),
        "reverse": compact_oligo(validation["reverse"]),
        "pair": {
            "core_delta_tm_c": pair["core_delta_tm_c"],
            "heterodimer_tm_c": pair["heterodimer"]["tm_c"],
            "heterodimer_dg_kcal_mol": pair["heterodimer"]["dg_kcal_mol"],
            "end_forward_dg_kcal_mol": pair[
                "end_stability_forward_3prime_against_reverse"
            ]["dg_kcal_mol"],
            "end_reverse_dg_kcal_mol": pair[
                "end_stability_reverse_3prime_against_forward"
            ]["dg_kcal_mol"],
        },
    }


def _hard_filter_reasons(
    compact: dict[str, Any],
    risks: dict[str, float],
    forward_site_count: int,
    reverse_site_count: int,
    config: MutagenesisConfig,
) -> list[str]:
    fwd, rev, pair = compact["forward"], compact["reverse"], compact["pair"]
    reasons: list[str] = []
    for label, oligo in (("forward", fwd), ("reverse", rev)):
        if not config.min_tm_c <= oligo["core_tm_c"] <= config.max_tm_c:
            reasons.append(f"{label}_core_tm")
        if (
            not config.min_gc_percent
            <= oligo["core_gc_percent"]
            <= config.max_gc_percent
        ):
            reasons.append(f"{label}_core_gc")
        if oligo["longest_homopolymer"] > config.max_homopolymer:
            reasons.append(f"{label}_homopolymer")
        if oligo["three_prime_gc_count"] > config.max_three_prime_gc:
            reasons.append(f"{label}_three_prime_gc")
    if pair["core_delta_tm_c"] > config.max_delta_tm_c:
        reasons.append("core_delta_tm")
    thresholds = {
        "hairpin_risk_kcal_mol": config.max_hairpin_risk_kcal_mol,
        "homodimer_risk_kcal_mol": config.max_homodimer_risk_kcal_mol,
        "heterodimer_risk_kcal_mol": config.max_heterodimer_risk_kcal_mol,
        "end_dimer_risk_kcal_mol": config.max_end_dimer_risk_kcal_mol,
    }
    for name, threshold in thresholds.items():
        if risks[name] > threshold:
            reasons.append(name.removesuffix("_kcal_mol"))
    if forward_site_count != 1:
        reasons.append("forward_exact_site_count")
    if reverse_site_count != 1:
        reasons.append("reverse_exact_site_count")
    return reasons


def _objective_vector(candidate: dict[str, Any]) -> list[float]:
    objectives = candidate["objectives"]
    return [
        objectives["tm_target_deviation_c"],
        objectives["core_delta_tm_c"],
        objectives["total_oligo_length_nt"],
        objectives["hairpin_risk_kcal_mol"],
        objectives["dimer_risk_kcal_mol"],
        objectives["end_dimer_risk_kcal_mol"],
        objectives["codon_penalty"],
        objectives["nucleotide_substitutions"],
    ]


def design_inverse_pcr_mutagenesis(
    record: Any,
    *,
    feature_name: str,
    aa_position: int,
    target_aa: str,
    config: MutagenesisConfig | None = None,
) -> dict[str, Any]:
    """Exhaustively enumerate adjacent inverse-PCR pairs and return a Pareto set."""
    config = config or MutagenesisConfig()
    if not is_circular(record):
        raise ValueError(
            "adjacent whole-plasmid inverse PCR requires a circular template"
        )
    sequence = normalize_dna(str(record.seq), name="template")
    if config.min_core_length < 1 or config.max_core_length < config.min_core_length:
        raise ValueError("invalid core-length range")
    if 2 * config.max_core_length > len(sequence) - 3:
        raise ValueError(
            "max core length makes forward and reverse annealing cores overlap"
        )
    if config.max_core_length + 3 > 60:
        raise ValueError(
            "maximum core plus the 3-nt mutation tail exceeds Primer3's "
            "60-nt full-oligo thermodynamic limit"
        )
    if not all(
        math.isfinite(value)
        for value in (
            config.min_tm_c,
            config.target_tm_c,
            config.max_tm_c,
            config.min_gc_percent,
            config.max_gc_percent,
            config.max_delta_tm_c,
            config.max_hairpin_risk_kcal_mol,
            config.max_homodimer_risk_kcal_mol,
            config.max_heterodimer_risk_kcal_mol,
            config.max_end_dimer_risk_kcal_mol,
        )
    ):
        raise ValueError("thermodynamic bounds and risk limits must be finite")
    if not (
        config.min_tm_c <= config.target_tm_c <= config.max_tm_c
        and config.min_tm_c <= config.max_tm_c
    ):
        raise ValueError("Tm bounds must satisfy min <= target <= max")
    if not (0.0 <= config.min_gc_percent <= config.max_gc_percent <= 100.0):
        raise ValueError("GC bounds must satisfy 0 <= min <= max <= 100")
    if not config.tail_sides:
        raise ValueError("at least one mutation-tail side is required")
    if config.verify_top < 0:
        raise ValueError("verify_top must be nonnegative")
    if config.pydna_limit < 1:
        raise ValueError("pydna_limit must be positive")
    if config.pydna_limit > config.min_core_length:
        raise ValueError("pydna_limit must not exceed the minimum core length")
    if config.max_homopolymer < 1:
        raise ValueError("max_homopolymer must be positive")
    if not 0 <= config.max_three_prime_gc <= 5:
        raise ValueError("max_three_prime_gc must be between 0 and 5")
    if any(
        threshold < 0
        for threshold in (
            config.max_delta_tm_c,
            config.max_hairpin_risk_kcal_mol,
            config.max_homodimer_risk_kcal_mol,
            config.max_heterodimer_risk_kcal_mol,
            config.max_end_dimer_risk_kcal_mol,
        )
    ):
        raise ValueError("delta-Tm and structure-risk limits must be nonnegative")
    if any(side not in {"forward", "reverse"} for side in config.tail_sides):
        raise ValueError("tail sides must be 'forward' and/or 'reverse'")
    if len(set(config.tail_sides)) != len(config.tail_sides):
        raise ValueError("tail sides must not contain duplicates")
    if len({codon.upper() for codon in config.codons}) != len(config.codons):
        raise ValueError("codons must not contain duplicates")
    feature = find_feature(record, feature_name)
    context = _target_context(record, feature, aa_position)
    target_aa = target_aa.upper()
    if target_aa == context["amino_acid"]:
        raise ValueError(
            "target amino acid must differ from the source amino acid; "
            "synonymous edits are outside this workflow"
        )
    start = context["genomic_start_0based"]
    end = context["genomic_end_0based_exclusive"]
    codons = _target_codons(target_aa, config, context["translation_table"])
    candidates: list[dict[str, Any]] = []
    for coding_codon, frequency in codons:
        mutant_genomic_codon = (
            coding_codon if context["strand"] == 1 else reverse_complement(coding_codon)
        )
        expected_sequence = sequence[:start] + mutant_genomic_codon + sequence[end:]
        mutant_coding_sequence = str(feature.extract(Seq(expected_sequence))).upper()[
            context["phase_bases"] :
        ]
        coding_offset = (aa_position - 1) * 3
        recovered_codon = mutant_coding_sequence[coding_offset : coding_offset + 3]
        translated_target = str(
            Seq(recovered_codon).translate(table=context["translation_table"])
        )
        if recovered_codon != coding_codon or translated_target != target_aa:
            raise AssertionError("mutant coding-sequence translation invariant failed")
        substitution_count = sum(
            left != right
            for left, right in zip(context["coding_codon"], coding_codon, strict=True)
        )
        for side in config.tail_sides:
            for forward_length in range(
                config.min_core_length, config.max_core_length + 1
            ):
                forward_core = _circular_slice(sequence, end, forward_length)
                for reverse_length in range(
                    config.min_core_length, config.max_core_length + 1
                ):
                    upstream = _circular_slice(
                        sequence, start - reverse_length, reverse_length
                    )
                    reverse_core = reverse_complement(upstream)
                    if side == "forward":
                        forward = mutant_genomic_codon + forward_core
                        reverse = reverse_core
                    else:
                        forward = forward_core
                        reverse = (
                            reverse_complement(mutant_genomic_codon) + reverse_core
                        )
                    validation = validate_pair(
                        forward,
                        reverse,
                        forward_core_length=len(forward_core),
                        reverse_core_length=len(reverse_core),
                        conditions=config.conditions,
                    )
                    compact = _compact_validation(validation)
                    risks = structural_risks(validation)
                    forward_sites = circular_occurrences(sequence, forward_core, True)
                    reverse_sites = circular_occurrences(
                        sequence, reverse_complement(reverse_core), True
                    )
                    reasons = _hard_filter_reasons(
                        compact,
                        risks,
                        len(forward_sites),
                        len(reverse_sites),
                        config,
                    )
                    f_tm = compact["forward"]["core_tm_c"]
                    r_tm = compact["reverse"]["core_tm_c"]
                    dimer_risk = max(
                        risks["homodimer_risk_kcal_mol"],
                        risks["heterodimer_risk_kcal_mol"],
                    )
                    objectives = {
                        "tm_target_deviation_c": abs(f_tm - config.target_tm_c)
                        + abs(r_tm - config.target_tm_c),
                        "core_delta_tm_c": compact["pair"]["core_delta_tm_c"],
                        "total_oligo_length_nt": len(forward) + len(reverse),
                        "hairpin_risk_kcal_mol": risks["hairpin_risk_kcal_mol"],
                        "dimer_risk_kcal_mol": dimer_risk,
                        "end_dimer_risk_kcal_mol": risks["end_dimer_risk_kcal_mol"],
                        "codon_penalty": 1.0 - frequency,
                        "nucleotide_substitutions": substitution_count,
                    }
                    score = (
                        objectives["tm_target_deviation_c"]
                        + 2.0 * objectives["core_delta_tm_c"]
                        + 0.05 * objectives["total_oligo_length_nt"]
                        + objectives["hairpin_risk_kcal_mol"]
                        + 0.5 * objectives["dimer_risk_kcal_mol"]
                        + objectives["end_dimer_risk_kcal_mol"]
                        + 2.0 * objectives["codon_penalty"]
                        + 0.5 * objectives["nucleotide_substitutions"]
                    )
                    candidates.append(
                        {
                            "id": f"{coding_codon}-{side}-F{forward_length}-R{reverse_length}",
                            "coding_codon": coding_codon,
                            "codon_frequency_within_amino_acid": frequency,
                            "mutation_tail_side": side,
                            "nucleotide_substitutions": substitution_count,
                            "forward": compact["forward"],
                            "reverse": compact["reverse"],
                            "pair": compact["pair"],
                            "risks": risks,
                            "exact_site_counts": {
                                "forward": len(forward_sites),
                                "reverse": len(reverse_sites),
                            },
                            "hard_filter_pass": not reasons,
                            "hard_filter_reasons": reasons,
                            "objectives": objectives,
                            "weighted_score": score,
                            "pareto_optimal": False,
                            "pydna_verification": None,
                            "_expected_sequence": expected_sequence,
                        }
                    )
    passing = [candidate for candidate in candidates if candidate["hard_filter_pass"]]
    front_local_indices = pareto_front_indices(
        [_objective_vector(candidate) for candidate in passing]
    )
    for index in front_local_indices:
        passing[index]["pareto_optimal"] = True
    pareto_candidates = sorted(
        (candidate for candidate in passing if candidate["pareto_optimal"]),
        key=lambda candidate: (candidate["weighted_score"], candidate["id"]),
    )
    for candidate in pareto_candidates[: config.verify_top]:
        try:
            verification, _ = simulate_product(
                record,
                candidate["forward"]["full_sequence_5to3"],
                candidate["reverse"]["full_sequence_5to3"],
                forward_core_length=candidate["forward"]["core_length_nt"],
                reverse_core_length=candidate["reverse"]["core_length_nt"],
                limit=config.pydna_limit,
                circularize=True,
                expected_sequence=candidate["_expected_sequence"],
            )
            candidate["pydna_verification"] = verification
        except ValueError as error:  # pydna reports ambiguous/non-products this way
            candidate["pydna_verification"] = {
                "expected_sequence_match": False,
                "error": f"{type(error).__name__}: {error}",
            }
    output_candidates = candidates if config.include_failed else passing
    for candidate in candidates:
        candidate.pop("_expected_sequence", None)
    target_codons = [codon for codon, _ in codons]
    return {
        "template": {
            "id": record.id,
            "length_bp": len(record),
            "circular": True,
        },
        "target": {
            "feature": feature_name,
            "feature_type": feature.type,
            "feature_start_1based": int(feature.location.start) + 1,
            "feature_end_1based_inclusive": int(feature.location.end),
            "feature_strand": feature.location.strand,
            "feature_codon_start": context["codon_start"],
            "feature_translation_table": context["translation_table"],
            "aa_position_1based": aa_position,
            "source_amino_acid": context["amino_acid"],
            "target_amino_acid": target_aa,
            "source_coding_codon": context["coding_codon"],
            "candidate_coding_codons": target_codons,
            "genomic_codon_start_1based": start + 1,
            "genomic_codon_end_1based_inclusive": end,
        },
        "config": config.as_dict(),
        "summary": {
            "enumerated": len(candidates),
            "hard_filter_passed": len(passing),
            "pareto_count": len(pareto_candidates),
            "pydna_verified_count": sum(
                candidate["pydna_verification"] is not None
                and candidate["pydna_verification"].get("expected_sequence_match")
                is True
                for candidate in pareto_candidates
            ),
        },
        "pareto_candidate_ids_by_weighted_score": [
            candidate["id"] for candidate in pareto_candidates
        ],
        "candidates": output_candidates,
        "tool": tool_versions(),
        "interpretation": {
            "global_scope": (
                "exhaustive only within the declared adjacent inverse-PCR search bounds"
            ),
            "pareto": "nondominated under the explicitly reported minimization objectives",
            "weighted_score": (
                "transparent tie-breaking aid, not a biological or experimental global optimum"
            ),
        },
    }
