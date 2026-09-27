# Copyright (c) 2026 Seungwon Lee
# SPDX-License-Identifier: MIT
"""Typed outputs for stateless primer workflow stages."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from pydna.dseqrecord import Dseqrecord

from .conditions import ReactionConditions
from .models import Oligo, PrimerPair
from .provenance import ToolVersions
from .serialization import JsonModel
from .thermo import StructureResult

if TYPE_CHECKING:
    from .mutagenesis import MutagenesisConfig


@dataclass(frozen=True, slots=True)
class TemplateResult(JsonModel):
    id: str
    length_bp: int
    circular: bool
    sha256: str | None = field(default=None, metadata={"omit_none": True})


@dataclass(frozen=True, slots=True)
class DesignedOligo(JsonModel):
    sequence_5to3: str
    location_0based_primer3: tuple[int, int]
    tm_c: float
    gc_percent: float

    @property
    def oligo(self) -> Oligo:
        return Oligo(self.sequence_5to3)


@dataclass(frozen=True, slots=True)
class DesignedPrimerPair(JsonModel):
    rank: int
    penalty: float
    product_size_bp: int
    forward: DesignedOligo
    reverse: DesignedOligo

    @property
    def primer_pair(self) -> PrimerPair:
        return PrimerPair(self.forward.oligo, self.reverse.oligo)


@dataclass(frozen=True, slots=True)
class DesignExplain(JsonModel):
    left: str | None
    right: str | None
    pair: str | None


@dataclass(frozen=True, slots=True)
class DesignInterpretation(JsonModel):
    ranking: str
    circular_template: str | None


@dataclass(frozen=True, slots=True)
class DesignResult(JsonModel):
    template: TemplateResult
    conditions: ReactionConditions
    pairs: tuple[DesignedPrimerPair, ...]
    explain: DesignExplain
    interpretation: DesignInterpretation
    tool: ToolVersions


@dataclass(frozen=True, slots=True)
class OligoValidationResult(JsonModel):
    full_sequence_5to3: str
    tail_5to3: str
    annealing_core_5to3: str
    full_length_nt: int
    tail_length_nt: int
    core_length_nt: int
    core_tm_c: float
    full_oligo_tm_c: float
    core_gc_percent: float
    full_gc_percent: float
    longest_homopolymer: int
    three_prime_pentamer: str
    three_prime_gc_count: int
    hairpin: StructureResult
    homodimer: StructureResult

    @property
    def oligo(self) -> Oligo:
        return Oligo(self.full_sequence_5to3, self.core_length_nt)


@dataclass(frozen=True, slots=True)
class PairThermodynamicsResult(JsonModel):
    core_delta_tm_c: float
    heterodimer: StructureResult
    end_stability_forward_3prime_against_reverse: StructureResult
    end_stability_reverse_3prime_against_forward: StructureResult


@dataclass(frozen=True, slots=True)
class ValidationInterpretation(JsonModel):
    annealing_tm_basis: str
    structure_basis: str
    warning: str


@dataclass(frozen=True, slots=True)
class StructuralRisks(JsonModel):
    hairpin_risk_kcal_mol: float
    homodimer_risk_kcal_mol: float
    heterodimer_risk_kcal_mol: float
    end_dimer_risk_kcal_mol: float


@dataclass(frozen=True, slots=True)
class PairValidationResult(JsonModel):
    conditions: ReactionConditions
    forward: OligoValidationResult
    reverse: OligoValidationResult
    pair: PairThermodynamicsResult
    interpretation: ValidationInterpretation
    tool: ToolVersions

    @property
    def primer_pair(self) -> PrimerPair:
        return PrimerPair(self.forward.oligo, self.reverse.oligo)

    @property
    def risks(self) -> StructuralRisks:
        hairpin_dgs = (
            self.forward.hairpin.dg_kcal_mol,
            self.reverse.hairpin.dg_kcal_mol,
        )
        homodimer_dgs = (
            self.forward.homodimer.dg_kcal_mol,
            self.reverse.homodimer.dg_kcal_mol,
        )
        end_dgs = (
            self.pair.end_stability_forward_3prime_against_reverse.dg_kcal_mol,
            self.pair.end_stability_reverse_3prime_against_forward.dg_kcal_mol,
        )
        return StructuralRisks(
            hairpin_risk_kcal_mol=max(0.0, -min(hairpin_dgs)),
            homodimer_risk_kcal_mol=max(0.0, -min(homodimer_dgs)),
            heterodimer_risk_kcal_mol=max(0.0, -self.pair.heterodimer.dg_kcal_mol),
            end_dimer_risk_kcal_mol=max(0.0, -min(end_dgs)),
        )


@dataclass(frozen=True, slots=True)
class ExactSite(JsonModel):
    start: int
    end: int


@dataclass(frozen=True, slots=True)
class ForwardSpecificityResult(JsonModel):
    tail_5to3: str
    core_5to3: str
    exact_core_sites_1based: tuple[ExactSite, ...]

    @property
    def oligo(self) -> Oligo:
        return Oligo(self.tail_5to3 + self.core_5to3, len(self.core_5to3))


@dataclass(frozen=True, slots=True)
class ReverseSpecificityResult(JsonModel):
    tail_5to3: str
    core_5to3: str
    plus_strand_binding_sequence_5to3: str
    exact_core_sites_1based: tuple[ExactSite, ...]

    @property
    def oligo(self) -> Oligo:
        return Oligo(self.tail_5to3 + self.core_5to3, len(self.core_5to3))


@dataclass(frozen=True, slots=True)
class PydnaProduct(JsonModel):
    length_bp: int
    forward_primer_role: str
    reverse_primer_role: str


@dataclass(frozen=True, slots=True)
class PydnaSpecificityResult(JsonModel):
    limit: int
    report: str
    product_count: int
    product_lengths_bp: tuple[int, ...]
    products: tuple[PydnaProduct, ...]


@dataclass(frozen=True, slots=True)
class SpecificityResult(JsonModel):
    template: TemplateResult
    forward: ForwardSpecificityResult
    reverse: ReverseSpecificityResult
    pydna: PydnaSpecificityResult
    scope_warning: str
    tool: ToolVersions

    @property
    def primer_pair(self) -> PrimerPair:
        return PrimerPair(self.forward.oligo, self.reverse.oligo)


@dataclass(frozen=True, slots=True)
class PrimerRoles(JsonModel):
    forward: str
    reverse: str


@dataclass(frozen=True, slots=True)
class ProductVerificationReport(JsonModel):
    product_length_bp: int
    circular: bool
    sha256: str
    expected_sequence_match: bool | None
    figure: str
    annealing_report: str
    primer_roles: PrimerRoles
    scope_warning: str
    tool: ToolVersions


@dataclass(frozen=True, slots=True)
class ProductVerificationResult(ProductVerificationReport):
    product: Dseqrecord = field(
        repr=False,
        compare=False,
        metadata={"serialize": False},
    )


@dataclass(frozen=True, slots=True)
class ProductVerificationFailure(JsonModel):
    expected_sequence_match: bool
    error: str


@dataclass(frozen=True, slots=True)
class CompactOligoValidation(JsonModel):
    full_sequence_5to3: str
    tail_5to3: str
    annealing_core_5to3: str
    full_length_nt: int
    tail_length_nt: int
    core_length_nt: int
    core_tm_c: float
    full_oligo_tm_c: float
    core_gc_percent: float
    full_gc_percent: float
    longest_homopolymer: int
    three_prime_pentamer: str
    three_prime_gc_count: int
    hairpin_tm_c: float
    hairpin_dg_kcal_mol: float
    homodimer_tm_c: float
    homodimer_dg_kcal_mol: float

    @property
    def oligo(self) -> Oligo:
        return Oligo(self.full_sequence_5to3, self.core_length_nt)


@dataclass(frozen=True, slots=True)
class CompactPairValidation(JsonModel):
    core_delta_tm_c: float
    heterodimer_tm_c: float
    heterodimer_dg_kcal_mol: float
    end_forward_dg_kcal_mol: float
    end_reverse_dg_kcal_mol: float


@dataclass(frozen=True, slots=True)
class ExactSiteCounts(JsonModel):
    forward: int
    reverse: int


@dataclass(frozen=True, slots=True)
class MutagenesisObjectives(JsonModel):
    tm_target_deviation_c: float
    core_delta_tm_c: float
    total_oligo_length_nt: int
    hairpin_risk_kcal_mol: float
    dimer_risk_kcal_mol: float
    end_dimer_risk_kcal_mol: float
    codon_penalty: float
    nucleotide_substitutions: int


@dataclass(frozen=True, slots=True)
class MutagenesisCandidate(JsonModel):
    id: str
    coding_codon: str
    codon_frequency_within_amino_acid: float
    mutation_tail_side: str
    nucleotide_substitutions: int
    forward: CompactOligoValidation
    reverse: CompactOligoValidation
    pair: CompactPairValidation
    risks: StructuralRisks
    exact_site_counts: ExactSiteCounts
    hard_filter_pass: bool
    hard_filter_reasons: tuple[str, ...]
    objectives: MutagenesisObjectives
    weighted_score: float
    pareto_optimal: bool
    pydna_verification: ProductVerificationReport | ProductVerificationFailure | None

    @property
    def primer_pair(self) -> PrimerPair:
        return PrimerPair(self.forward.oligo, self.reverse.oligo)


@dataclass(frozen=True, slots=True)
class MutagenesisTargetResult(JsonModel):
    feature: str
    feature_type: str
    feature_start_1based: int
    feature_end_1based_inclusive: int
    feature_strand: int
    feature_codon_start: int
    feature_translation_table: int
    aa_position_1based: int
    source_amino_acid: str
    target_amino_acid: str
    source_coding_codon: str
    candidate_coding_codons: tuple[str, ...]
    genomic_codon_start_1based: int
    genomic_codon_end_1based_inclusive: int


@dataclass(frozen=True, slots=True)
class MutagenesisSummary(JsonModel):
    enumerated: int
    hard_filter_passed: int
    pareto_count: int
    pydna_verified_count: int


@dataclass(frozen=True, slots=True)
class MutagenesisInterpretation(JsonModel):
    global_scope: str
    pareto: str
    weighted_score: str


@dataclass(frozen=True, slots=True)
class MutagenesisResult(JsonModel):
    template: TemplateResult
    target: MutagenesisTargetResult
    config: MutagenesisConfig
    summary: MutagenesisSummary
    pareto_candidate_ids_by_weighted_score: tuple[str, ...]
    candidates: tuple[MutagenesisCandidate, ...]
    tool: ToolVersions
    interpretation: MutagenesisInterpretation


def _structure(payload: dict[str, Any]) -> StructureResult:
    return StructureResult(**payload)


def _tool(payload: dict[str, Any]) -> ToolVersions:
    return ToolVersions(**payload)


def _template(payload: dict[str, Any]) -> TemplateResult:
    return TemplateResult(**payload)


def _primer3_location(values: list[Any] | tuple[Any, ...]) -> tuple[int, int]:
    if len(values) != 2:
        raise ValueError("Primer3 location must contain start and length")
    return int(values[0]), int(values[1])


def design_result_from_payload(payload: dict[str, Any]) -> DesignResult:
    pairs = tuple(
        DesignedPrimerPair(
            rank=int(item["rank"]),
            penalty=float(item["penalty"]),
            product_size_bp=int(item["product_size_bp"]),
            forward=DesignedOligo(
                sequence_5to3=str(item["forward"]["sequence_5to3"]),
                location_0based_primer3=_primer3_location(
                    item["forward"]["location_0based_primer3"]
                ),
                tm_c=float(item["forward"]["tm_c"]),
                gc_percent=float(item["forward"]["gc_percent"]),
            ),
            reverse=DesignedOligo(
                sequence_5to3=str(item["reverse"]["sequence_5to3"]),
                location_0based_primer3=_primer3_location(
                    item["reverse"]["location_0based_primer3"]
                ),
                tm_c=float(item["reverse"]["tm_c"]),
                gc_percent=float(item["reverse"]["gc_percent"]),
            ),
        )
        for item in payload["pairs"]
    )
    return DesignResult(
        template=_template(payload["template"]),
        conditions=ReactionConditions(**payload["conditions"]),
        pairs=pairs,
        explain=DesignExplain(**payload["explain"]),
        interpretation=DesignInterpretation(**payload["interpretation"]),
        tool=_tool(payload["tool"]),
    )


def _oligo_validation(payload: dict[str, Any]) -> OligoValidationResult:
    values = dict(payload)
    values["hairpin"] = _structure(values["hairpin"])
    values["homodimer"] = _structure(values["homodimer"])
    return OligoValidationResult(**values)


def validation_result_from_payload(payload: dict[str, Any]) -> PairValidationResult:
    pair = dict(payload["pair"])
    pair["heterodimer"] = _structure(pair["heterodimer"])
    pair["end_stability_forward_3prime_against_reverse"] = _structure(
        pair["end_stability_forward_3prime_against_reverse"]
    )
    pair["end_stability_reverse_3prime_against_forward"] = _structure(
        pair["end_stability_reverse_3prime_against_forward"]
    )
    return PairValidationResult(
        conditions=ReactionConditions(**payload["conditions"]),
        forward=_oligo_validation(payload["forward"]),
        reverse=_oligo_validation(payload["reverse"]),
        pair=PairThermodynamicsResult(**pair),
        interpretation=ValidationInterpretation(**payload["interpretation"]),
        tool=_tool(payload["tool"]),
    )


def _sites(payload: list[dict[str, Any]]) -> tuple[ExactSite, ...]:
    return tuple(ExactSite(**site) for site in payload)


def specificity_result_from_payload(payload: dict[str, Any]) -> SpecificityResult:
    forward = dict(payload["forward"])
    forward["exact_core_sites_1based"] = _sites(forward["exact_core_sites_1based"])
    reverse = dict(payload["reverse"])
    reverse["exact_core_sites_1based"] = _sites(reverse["exact_core_sites_1based"])
    pydna = dict(payload["pydna"])
    pydna["product_lengths_bp"] = tuple(pydna["product_lengths_bp"])
    pydna["products"] = tuple(PydnaProduct(**item) for item in pydna["products"])
    return SpecificityResult(
        template=_template(payload["template"]),
        forward=ForwardSpecificityResult(**forward),
        reverse=ReverseSpecificityResult(**reverse),
        pydna=PydnaSpecificityResult(**pydna),
        scope_warning=str(payload["scope_warning"]),
        tool=_tool(payload["tool"]),
    )


def _verification_report(payload: dict[str, Any]) -> ProductVerificationReport:
    values = dict(payload)
    values["primer_roles"] = PrimerRoles(**values["primer_roles"])
    values["tool"] = _tool(values["tool"])
    return ProductVerificationReport(**values)


def verification_result_from_payload(
    payload: dict[str, Any], product: Dseqrecord
) -> ProductVerificationResult:
    report = _verification_report(payload)
    return ProductVerificationResult(
        product_length_bp=report.product_length_bp,
        circular=report.circular,
        sha256=report.sha256,
        expected_sequence_match=report.expected_sequence_match,
        figure=report.figure,
        annealing_report=report.annealing_report,
        primer_roles=report.primer_roles,
        scope_warning=report.scope_warning,
        tool=report.tool,
        product=product,
    )


def _candidate_verification(
    payload: dict[str, Any] | None,
) -> ProductVerificationReport | ProductVerificationFailure | None:
    if payload is None:
        return None
    if "error" in payload:
        return ProductVerificationFailure(**payload)
    return _verification_report(payload)


def _mutagenesis_candidate(payload: dict[str, Any]) -> MutagenesisCandidate:
    return MutagenesisCandidate(
        id=str(payload["id"]),
        coding_codon=str(payload["coding_codon"]),
        codon_frequency_within_amino_acid=float(
            payload["codon_frequency_within_amino_acid"]
        ),
        mutation_tail_side=str(payload["mutation_tail_side"]),
        nucleotide_substitutions=int(payload["nucleotide_substitutions"]),
        forward=CompactOligoValidation(**payload["forward"]),
        reverse=CompactOligoValidation(**payload["reverse"]),
        pair=CompactPairValidation(**payload["pair"]),
        risks=StructuralRisks(**payload["risks"]),
        exact_site_counts=ExactSiteCounts(**payload["exact_site_counts"]),
        hard_filter_pass=bool(payload["hard_filter_pass"]),
        hard_filter_reasons=tuple(payload["hard_filter_reasons"]),
        objectives=MutagenesisObjectives(**payload["objectives"]),
        weighted_score=float(payload["weighted_score"]),
        pareto_optimal=bool(payload["pareto_optimal"]),
        pydna_verification=_candidate_verification(payload["pydna_verification"]),
    )


def mutagenesis_result_from_payload(
    payload: dict[str, Any], config: MutagenesisConfig
) -> MutagenesisResult:
    target = dict(payload["target"])
    target["candidate_coding_codons"] = tuple(target["candidate_coding_codons"])
    return MutagenesisResult(
        template=_template(payload["template"]),
        target=MutagenesisTargetResult(**target),
        config=config,
        summary=MutagenesisSummary(**payload["summary"]),
        pareto_candidate_ids_by_weighted_score=tuple(
            payload["pareto_candidate_ids_by_weighted_score"]
        ),
        candidates=tuple(
            _mutagenesis_candidate(candidate) for candidate in payload["candidates"]
        ),
        tool=_tool(payload["tool"]),
        interpretation=MutagenesisInterpretation(**payload["interpretation"]),
    )
