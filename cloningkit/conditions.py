# Copyright (c) 2026 Seungwon Lee
# SPDX-License-Identifier: MIT
"""Reaction-condition models shared by all thermodynamic calculations."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

from .serialization import JsonModel

type TmMethod = Literal["breslauer", "santalucia"]
type SaltCorrectionMethod = Literal["schildkraut", "santalucia", "owczarzy"]


@dataclass(frozen=True, slots=True)
class ReactionConditions(JsonModel):
    """Primer3 reaction conditions in the units used by primer3-py.

    Monovalent, divalent, and dNTP concentrations are millimolar; oligo
    concentration is nanomolar; temperatures are degrees Celsius.
    """

    mv_conc: float = 50.0
    dv_conc: float = 1.5
    dntp_conc: float = 0.6
    dna_conc: float = 50.0
    annealing_temp_c: float = 60.0
    tm_method: TmMethod = "santalucia"
    salt_corrections_method: SaltCorrectionMethod = "santalucia"

    def __post_init__(self) -> None:
        concentrations = {
            "mv_conc": self.mv_conc,
            "dv_conc": self.dv_conc,
            "dntp_conc": self.dntp_conc,
        }
        if any(
            not math.isfinite(value) or value < 0 for value in concentrations.values()
        ):
            raise ValueError(
                "salt and dNTP concentrations must be finite and nonnegative"
            )
        if not math.isfinite(self.dna_conc) or self.dna_conc <= 0:
            raise ValueError("oligo concentration must be finite and positive")
        if (
            not math.isfinite(self.annealing_temp_c)
            or not 0 <= self.annealing_temp_c <= 100
        ):
            raise ValueError(
                "annealing temperature must be finite and between 0 and 100 C"
            )
        if self.tm_method not in {"breslauer", "santalucia"}:
            raise ValueError("unsupported Primer3 Tm method")
        if self.salt_corrections_method not in {
            "schildkraut",
            "santalucia",
            "owczarzy",
        }:
            raise ValueError("unsupported Primer3 salt-correction method")

    def tm_kwargs(self) -> dict[str, float | str]:
        return {
            "mv_conc": self.mv_conc,
            "dv_conc": self.dv_conc,
            "dntp_conc": self.dntp_conc,
            "dna_conc": self.dna_conc,
            "tm_method": self.tm_method,
            "salt_corrections_method": self.salt_corrections_method,
        }

    def structure_kwargs(self) -> dict[str, float]:
        return {
            "mv_conc": self.mv_conc,
            "dv_conc": self.dv_conc,
            "dntp_conc": self.dntp_conc,
            "dna_conc": self.dna_conc,
            "temp_c": self.annealing_temp_c,
        }
