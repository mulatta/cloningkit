# Copyright (c) 2026 Seungwon Lee
# SPDX-License-Identifier: MIT
"""Typed inputs and configuration for stateless primer workflow stages."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .conditions import ReactionConditions
from .io import normalize_dna
from .serialization import JsonModel


@dataclass(frozen=True, slots=True)
class Oligo(JsonModel):
    """Full 5'-to-3' oligo with a 3'-terminal template-binding core."""

    sequence_5to3: str
    core_length: int | None = None

    def __post_init__(self) -> None:
        sequence = normalize_dna(self.sequence_5to3, name="oligo")
        object.__setattr__(self, "sequence_5to3", sequence)
        if self.core_length is not None and not 1 <= self.core_length <= len(sequence):
            raise ValueError("core length must be between 1 and full oligo length")

    @property
    def resolved_core_length(self) -> int:
        return len(self.sequence_5to3) if self.core_length is None else self.core_length

    @property
    def tail_5to3(self) -> str:
        return self.sequence_5to3[: -self.resolved_core_length]

    @property
    def annealing_core_5to3(self) -> str:
        return self.sequence_5to3[-self.resolved_core_length :]


@dataclass(frozen=True, slots=True)
class PrimerPair(JsonModel):
    """Named forward and reverse oligos in their declared PCR roles."""

    forward: Oligo
    reverse: Oligo

    def __post_init__(self) -> None:
        if not isinstance(self.forward, Oligo) or not isinstance(self.reverse, Oligo):
            raise TypeError("forward and reverse must be Oligo values")


@dataclass(frozen=True, slots=True)
class DesignConfig(JsonModel):
    """Primer3 bounds for conventional inward-facing primer design."""

    target_start_1based: int | None = None
    target_length: int | None = None
    product_min_bp: int = 100
    product_max_bp: int = 1000
    min_size_nt: int = 18
    opt_size_nt: int = 20
    max_size_nt: int = 25
    min_tm_c: float = 58.0
    opt_tm_c: float = 60.0
    max_tm_c: float = 62.0
    min_gc_percent: float = 40.0
    max_gc_percent: float = 60.0
    max_delta_tm_c: float = 2.0
    num_return: int = 5
    conditions: ReactionConditions = field(default_factory=ReactionConditions)

    def __post_init__(self) -> None:
        if not 1 <= self.product_min_bp <= self.product_max_bp:
            raise ValueError("product bounds must satisfy 1 <= min <= max")
        if not 1 <= self.min_size_nt <= self.opt_size_nt <= self.max_size_nt:
            raise ValueError("primer-size bounds must satisfy 1 <= min <= opt <= max")
        if not all(
            math.isfinite(value)
            for value in (self.min_tm_c, self.opt_tm_c, self.max_tm_c)
        ):
            raise ValueError("Tm bounds must be finite")
        if not self.min_tm_c <= self.opt_tm_c <= self.max_tm_c:
            raise ValueError("Tm bounds must satisfy min <= opt <= max")
        if not 0 <= self.min_gc_percent <= self.max_gc_percent <= 100:
            raise ValueError("GC bounds must satisfy 0 <= min <= max <= 100")
        if not math.isfinite(self.max_delta_tm_c) or self.max_delta_tm_c < 0:
            raise ValueError("maximum pair delta-Tm must be finite and nonnegative")
        if self.num_return < 1:
            raise ValueError("num_return must be positive")
        if (self.target_start_1based is None) != (self.target_length is None):
            raise ValueError("target start and target length must be supplied together")
        if self.target_start_1based is not None and self.target_start_1based < 1:
            raise ValueError("target start and length must be positive")
        if self.target_length is not None and self.target_length < 1:
            raise ValueError("target start and length must be positive")
        if not isinstance(self.conditions, ReactionConditions):
            raise TypeError("conditions must be ReactionConditions")


@dataclass(frozen=True, slots=True)
class SpecificityConfig(JsonModel):
    """Exact-template and pydna product-enumeration settings."""

    pydna_limit: int = 13

    def __post_init__(self) -> None:
        if self.pydna_limit < 1:
            raise ValueError("pydna limit must be positive")


@dataclass(frozen=True, slots=True)
class ProductVerificationConfig(JsonModel):
    """Settings for unique-product simulation and optional sequence comparison."""

    pydna_limit: int = 13
    circularize: bool = False
    expected_sequence: str | None = None

    def __post_init__(self) -> None:
        if self.pydna_limit < 1:
            raise ValueError("pydna limit must be positive")


@dataclass(frozen=True, slots=True)
class MutationTarget(JsonModel):
    """One annotated residue substitution requested from mutagenesis design."""

    feature_name: str
    aa_position_1based: int
    target_amino_acid: str

    def __post_init__(self) -> None:
        if not self.feature_name:
            raise ValueError("feature name must not be empty")
        if self.aa_position_1based < 1:
            raise ValueError("amino-acid position is 1-based and must be positive")
        target_amino_acid = self.target_amino_acid.upper()
        if len(target_amino_acid) != 1:
            raise ValueError("target amino acid must be a one-letter code")
        object.__setattr__(self, "target_amino_acid", target_amino_acid)
