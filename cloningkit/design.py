# Copyright (c) 2026 Seungwon Lee
# SPDX-License-Identifier: MIT
"""Conventional inward-facing Primer3 design."""

from __future__ import annotations

import math
from typing import Any

import primer3

from .conditions import ReactionConditions
from .io import is_circular, normalize_dna
from .provenance import tool_versions


def _plain(value: Any) -> Any:
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_plain(item) for item in value]
    return value


def design_standard_primers(
    record: Any,
    *,
    target_start_1based: int | None = None,
    target_length: int | None = None,
    product_min: int = 100,
    product_max: int = 1000,
    min_size: int = 18,
    opt_size: int = 20,
    max_size: int = 25,
    min_tm: float = 58.0,
    opt_tm: float = 60.0,
    max_tm: float = 62.0,
    min_gc: float = 40.0,
    max_gc: float = 60.0,
    max_delta_tm: float = 2.0,
    num_return: int = 5,
    conditions: ReactionConditions | None = None,
) -> dict[str, Any]:
    """Design ranked conventional PCR pairs on a linear template."""
    circular = is_circular(record)
    sequence = normalize_dna(str(record.seq), name="template")
    if not 1 <= product_min <= product_max:
        raise ValueError("product bounds must satisfy 1 <= min <= max")
    if not 1 <= min_size <= opt_size <= max_size:
        raise ValueError("primer-size bounds must satisfy 1 <= min <= opt <= max")
    if not all(math.isfinite(value) for value in (min_tm, opt_tm, max_tm)):
        raise ValueError("Tm bounds must be finite")
    if not min_tm <= opt_tm <= max_tm:
        raise ValueError("Tm bounds must satisfy min <= opt <= max")
    if not 0 <= min_gc <= max_gc <= 100:
        raise ValueError("GC bounds must satisfy 0 <= min <= max <= 100")
    if not math.isfinite(max_delta_tm) or max_delta_tm < 0:
        raise ValueError("maximum pair delta-Tm must be finite and nonnegative")
    if num_return < 1:
        raise ValueError("num_return must be positive")
    if (target_start_1based is None) != (target_length is None):
        raise ValueError("target start and target length must be supplied together")
    seq_args: dict[str, Any] = {
        "SEQUENCE_ID": record.id,
        "SEQUENCE_TEMPLATE": sequence,
    }
    if target_start_1based is not None and target_length is not None:
        if target_start_1based < 1 or target_length < 1:
            raise ValueError("target start and length must be positive")
        if target_start_1based - 1 + target_length > len(sequence):
            raise ValueError("target interval lies outside the template")
        seq_args["SEQUENCE_TARGET"] = [target_start_1based - 1, target_length]
    conditions = conditions or ReactionConditions()
    if conditions.salt_corrections_method == "owczarzy":
        raise ValueError(
            "Primer3 design cannot combine PRIMER_ANNEALING_TEMP with "
            "Owczarzy salt correction; use SantaLucia/Schildkraut or run "
            "separate validation under the intended model"
        )
    global_args = {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_NUM_RETURN": num_return,
        "PRIMER_OPT_SIZE": opt_size,
        "PRIMER_MIN_SIZE": min_size,
        "PRIMER_MAX_SIZE": max_size,
        "PRIMER_OPT_TM": opt_tm,
        "PRIMER_MIN_TM": min_tm,
        "PRIMER_MAX_TM": max_tm,
        "PRIMER_PAIR_MAX_DIFF_TM": max_delta_tm,
        "PRIMER_MIN_GC": min_gc,
        "PRIMER_MAX_GC": max_gc,
        "PRIMER_PRODUCT_SIZE_RANGE": [[product_min, product_max]],
        "PRIMER_SALT_MONOVALENT": conditions.mv_conc,
        "PRIMER_SALT_DIVALENT": conditions.dv_conc,
        "PRIMER_DNTP_CONC": conditions.dntp_conc,
        "PRIMER_DNA_CONC": conditions.dna_conc,
        "PRIMER_TM_FORMULA": {"breslauer": 0, "santalucia": 1}[conditions.tm_method],
        "PRIMER_SALT_CORRECTIONS": {
            "schildkraut": 0,
            "santalucia": 1,
            "owczarzy": 2,
        }[conditions.salt_corrections_method],
        "PRIMER_ANNEALING_TEMP": conditions.annealing_temp_c,
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
        "PRIMER_EXPLAIN_FLAG": 1,
    }
    raw = primer3.design_primers(seq_args, global_args)
    count = int(raw.get("PRIMER_PAIR_NUM_RETURNED", 0))
    pairs = []
    for index in range(count):
        pairs.append(
            {
                "rank": index,
                "penalty": raw[f"PRIMER_PAIR_{index}_PENALTY"],
                "product_size_bp": raw[f"PRIMER_PAIR_{index}_PRODUCT_SIZE"],
                "forward": {
                    "sequence_5to3": raw[f"PRIMER_LEFT_{index}_SEQUENCE"],
                    "location_0based_primer3": _plain(raw[f"PRIMER_LEFT_{index}"]),
                    "tm_c": raw[f"PRIMER_LEFT_{index}_TM"],
                    "gc_percent": raw[f"PRIMER_LEFT_{index}_GC_PERCENT"],
                },
                "reverse": {
                    "sequence_5to3": raw[f"PRIMER_RIGHT_{index}_SEQUENCE"],
                    "location_0based_primer3": _plain(raw[f"PRIMER_RIGHT_{index}"]),
                    "tm_c": raw[f"PRIMER_RIGHT_{index}_TM"],
                    "gc_percent": raw[f"PRIMER_RIGHT_{index}_GC_PERCENT"],
                },
            }
        )
    return {
        "template": {"id": record.id, "length_bp": len(record), "circular": circular},
        "conditions": conditions.as_dict(),
        "pairs": pairs,
        "explain": {
            "left": raw.get("PRIMER_LEFT_EXPLAIN"),
            "right": raw.get("PRIMER_RIGHT_EXPLAIN"),
            "pair": raw.get("PRIMER_PAIR_EXPLAIN"),
        },
        "interpretation": {
            "ranking": (
                "rank 0 minimizes the configured Primer3 penalty on this template; "
                "it is not a genome-specificity or experimental guarantee"
            ),
            "circular_template": (
                "Primer3 searches the supplied linear representation and does not "
                "design origin-spanning products; rotate the record when needed"
                if circular
                else None
            ),
        },
        "tool": tool_versions(),
    }
