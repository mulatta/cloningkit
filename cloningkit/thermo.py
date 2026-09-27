# Copyright (c) 2026 Seungwon Lee
# SPDX-License-Identifier: MIT
"""Primer3-py thermodynamic validation with explicit core/full separation."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise
from typing import Any

import primer3

from .conditions import ReactionConditions
from .io import normalize_dna
from .provenance import tool_versions
from .serialization import JsonModel


@dataclass(frozen=True, slots=True)
class StructureResult(JsonModel):
    structure_found: bool
    tm_c: float
    dg_kcal_mol: float
    dh_kcal_mol: float
    ds_cal_mol_k: float
    ascii_structure: str | None = None

    @classmethod
    def from_primer3(cls, result: Any) -> StructureResult:
        return cls(
            structure_found=bool(result.structure_found),
            tm_c=float(result.tm),
            dg_kcal_mol=float(result.dg) / 1000.0,
            dh_kcal_mol=float(result.dh) / 1000.0,
            ds_cal_mol_k=float(result.ds),
            ascii_structure=getattr(result, "ascii_structure", None),
        )


def gc_percent(sequence: str) -> float:
    sequence = normalize_dna(sequence)
    return 100.0 * sum(base in "GC" for base in sequence) / len(sequence)


def longest_homopolymer(sequence: str) -> int:
    sequence = normalize_dna(sequence)
    longest = current = 1
    for previous, base in pairwise(sequence):
        current = current + 1 if previous == base else 1
        longest = max(longest, current)
    return longest


def split_core(full_sequence: str, core_length: int | None) -> tuple[str, str]:
    full_sequence = normalize_dna(full_sequence, name="primer")
    if core_length is None:
        core_length = len(full_sequence)
    if not 1 <= core_length <= len(full_sequence):
        raise ValueError("core length must be between 1 and full oligo length")
    return full_sequence[:-core_length], full_sequence[-core_length:]


def _oligo_metrics(
    full_sequence: str,
    core_length: int | None,
    conditions: ReactionConditions,
) -> dict[str, Any]:
    tail, core = split_core(full_sequence, core_length)
    structure_kwargs = conditions.structure_kwargs()
    hairpin = StructureResult.from_primer3(
        primer3.calc_hairpin(full_sequence, output_structure=True, **structure_kwargs)
    )
    homodimer = StructureResult.from_primer3(
        primer3.calc_homodimer(full_sequence, output_structure=True, **structure_kwargs)
    )
    return {
        "full_sequence_5to3": full_sequence,
        "tail_5to3": tail,
        "annealing_core_5to3": core,
        "full_length_nt": len(full_sequence),
        "tail_length_nt": len(tail),
        "core_length_nt": len(core),
        "core_tm_c": float(primer3.calc_tm(core, **conditions.tm_kwargs())),
        "full_oligo_tm_c": float(
            primer3.calc_tm(full_sequence, **conditions.tm_kwargs())
        ),
        "core_gc_percent": gc_percent(core),
        "full_gc_percent": gc_percent(full_sequence),
        "longest_homopolymer": longest_homopolymer(full_sequence),
        "three_prime_pentamer": core[-5:],
        "three_prime_gc_count": sum(base in "GC" for base in core[-5:]),
        "hairpin": hairpin.as_dict(),
        "homodimer": homodimer.as_dict(),
    }


def validate_pair(
    forward: str,
    reverse: str,
    *,
    forward_core_length: int | None = None,
    reverse_core_length: int | None = None,
    conditions: ReactionConditions | None = None,
) -> dict[str, Any]:
    """Validate a primer pair using core Tm and full-oligo structures.

    The annealing core is always the 3' suffix of the full oligo. Non-templated
    5' tails are excluded from pair Tm matching and included in all structure
    calculations.
    """
    conditions = conditions or ReactionConditions()
    forward = normalize_dna(forward, name="forward primer")
    reverse = normalize_dna(reverse, name="reverse primer")
    if len(forward) > 60 or len(reverse) > 60:
        raise ValueError(
            "Primer3 thermodynamic structure calculations require each full "
            "oligo to be at most 60 nt"
        )
    fwd = _oligo_metrics(forward, forward_core_length, conditions)
    rev = _oligo_metrics(reverse, reverse_core_length, conditions)
    kwargs = conditions.structure_kwargs()
    heterodimer = StructureResult.from_primer3(
        primer3.calc_heterodimer(forward, reverse, output_structure=True, **kwargs)
    )
    end_fwd = StructureResult.from_primer3(
        primer3.calc_end_stability(forward, reverse, **kwargs)
    )
    end_rev = StructureResult.from_primer3(
        primer3.calc_end_stability(reverse, forward, **kwargs)
    )
    delta_tm = abs(float(fwd["core_tm_c"]) - float(rev["core_tm_c"]))
    return {
        "conditions": conditions.as_dict(),
        "forward": fwd,
        "reverse": rev,
        "pair": {
            "core_delta_tm_c": delta_tm,
            "heterodimer": heterodimer.as_dict(),
            "end_stability_forward_3prime_against_reverse": end_fwd.as_dict(),
            "end_stability_reverse_3prime_against_forward": end_rev.as_dict(),
        },
        "interpretation": {
            "annealing_tm_basis": "template-binding 3-prime core only",
            "structure_basis": "full oligo including non-annealing 5-prime tails",
            "warning": "model-dependent prediction; inspect structures at the actual annealing temperature",
        },
        "tool": tool_versions(),
    }


def structural_risks(validation: dict[str, Any]) -> dict[str, float]:
    """Return one-sided risk magnitudes; zero means no favorable folding at Ta."""
    fwd = validation["forward"]
    rev = validation["reverse"]
    pair = validation["pair"]
    hairpin_dgs = [fwd["hairpin"]["dg_kcal_mol"], rev["hairpin"]["dg_kcal_mol"]]
    homo_dgs = [
        fwd["homodimer"]["dg_kcal_mol"],
        rev["homodimer"]["dg_kcal_mol"],
    ]
    end_dgs = [
        pair["end_stability_forward_3prime_against_reverse"]["dg_kcal_mol"],
        pair["end_stability_reverse_3prime_against_forward"]["dg_kcal_mol"],
    ]
    return {
        "hairpin_risk_kcal_mol": max(0.0, -min(hairpin_dgs)),
        "homodimer_risk_kcal_mol": max(0.0, -min(homo_dgs)),
        "heterodimer_risk_kcal_mol": max(0.0, -pair["heterodimer"]["dg_kcal_mol"]),
        "end_dimer_risk_kcal_mol": max(0.0, -min(end_dgs)),
    }
