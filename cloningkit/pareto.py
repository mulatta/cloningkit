# Copyright (c) 2026 Seungwon Lee
# SPDX-License-Identifier: MIT
"""Transparent Pareto filtering utilities."""

from __future__ import annotations

from collections.abc import Sequence


def dominates(left: Sequence[float], right: Sequence[float]) -> bool:
    """Return true when left is no worse in all objectives and better in one."""
    if len(left) != len(right):
        raise ValueError("objective vectors must have equal lengths")
    return all(a <= b for a, b in zip(left, right, strict=True)) and any(
        a < b for a, b in zip(left, right, strict=True)
    )


def pareto_front_indices(vectors: Sequence[Sequence[float]]) -> list[int]:
    """Return nondominated indices for minimization objectives."""
    front: list[int] = []
    for index, vector in enumerate(vectors):
        if not any(
            other_index != index and dominates(other, vector)
            for other_index, other in enumerate(vectors)
        ):
            front.append(index)
    return front
