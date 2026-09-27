# Copyright (c) 2026 Seungwon Lee
# SPDX-License-Identifier: MIT
"""Stateless typed composition API for Primer Workbench workflows."""

from __future__ import annotations

from Bio.SeqRecord import SeqRecord

from .conditions import ReactionConditions
from .design import design_standard_primers as _design_standard_primers
from .models import (
    DesignConfig,
    MutationTarget,
    PrimerPair,
    ProductVerificationConfig,
    SpecificityConfig,
)
from .mutagenesis import (
    MutagenesisConfig,
)
from .mutagenesis import (
    design_inverse_pcr_mutagenesis as _design_inverse_pcr_mutagenesis,
)
from .results import (
    DesignResult,
    MutagenesisResult,
    PairValidationResult,
    ProductVerificationResult,
    SpecificityResult,
    design_result_from_payload,
    mutagenesis_result_from_payload,
    specificity_result_from_payload,
    validation_result_from_payload,
    verification_result_from_payload,
)
from .specificity import (
    analyze_template_specificity as _analyze_template_specificity,
)
from .specificity import simulate_product as _simulate_product
from .thermo import validate_pair as _validate_pair


def design_primers(
    record: SeqRecord, config: DesignConfig | None = None
) -> DesignResult:
    """Design conventional PCR pairs under one immutable configuration."""
    config = config or DesignConfig()
    payload = _design_standard_primers(
        record,
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
    return design_result_from_payload(payload)


def validate_pair(
    pair: PrimerPair, conditions: ReactionConditions | None = None
) -> PairValidationResult:
    """Validate core Tm and full-oligo structures for one named pair."""
    conditions = conditions or ReactionConditions()
    payload = _validate_pair(
        pair.forward.sequence_5to3,
        pair.reverse.sequence_5to3,
        forward_core_length=pair.forward.resolved_core_length,
        reverse_core_length=pair.reverse.resolved_core_length,
        conditions=conditions,
    )
    return validation_result_from_payload(payload)


def analyze_specificity(
    record: SeqRecord,
    pair: PrimerPair,
    config: SpecificityConfig | None = None,
) -> SpecificityResult:
    """Analyze exact binding sites and all pydna products on one template."""
    config = config or SpecificityConfig()
    payload = _analyze_template_specificity(
        record,
        pair.forward.sequence_5to3,
        pair.reverse.sequence_5to3,
        forward_core_length=pair.forward.resolved_core_length,
        reverse_core_length=pair.reverse.resolved_core_length,
        limit=config.pydna_limit,
    )
    return specificity_result_from_payload(payload)


def verify_product(
    record: SeqRecord,
    pair: PrimerPair,
    config: ProductVerificationConfig | None = None,
) -> ProductVerificationResult:
    """Simulate one role-correct PCR product and retain its sequence artifact."""
    config = config or ProductVerificationConfig()
    payload, product = _simulate_product(
        record,
        pair.forward.sequence_5to3,
        pair.reverse.sequence_5to3,
        forward_core_length=pair.forward.resolved_core_length,
        reverse_core_length=pair.reverse.resolved_core_length,
        limit=config.pydna_limit,
        circularize=config.circularize,
        expected_sequence=config.expected_sequence,
    )
    return verification_result_from_payload(payload, product)


def design_mutagenesis(
    record: SeqRecord,
    target: MutationTarget,
    config: MutagenesisConfig | None = None,
) -> MutagenesisResult:
    """Compose exhaustive inverse-PCR design into an immutable typed result."""
    config = config or MutagenesisConfig()
    payload = _design_inverse_pcr_mutagenesis(
        record,
        feature_name=target.feature_name,
        aa_position=target.aa_position_1based,
        target_aa=target.target_amino_acid,
        config=config,
    )
    return mutagenesis_result_from_payload(payload, config)
