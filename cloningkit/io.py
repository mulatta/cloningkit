# Copyright (c) 2026 Seungwon Lee
# SPDX-License-Identifier: MIT
"""Input/output and sequence-coordinate helpers."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import cast

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .serialization import JsonModel


def load_record(path: str | Path) -> SeqRecord:
    """Load exactly one FASTA or GenBank record and preserve its annotations."""
    source = Path(path)
    suffix = source.suffix.lower()
    fmt = "genbank" if suffix in {".gb", ".gbk", ".genbank"} else "fasta"
    records = list(SeqIO.parse(source, fmt))
    if len(records) != 1:
        raise ValueError(f"expected exactly one {fmt} record, found {len(records)}")
    record = cast(SeqRecord, records[0])
    record.annotations.setdefault("molecule_type", "DNA")
    return record


def is_circular(record: SeqRecord) -> bool:
    return str(record.annotations.get("topology", "")).lower() == "circular"


def reverse_complement(sequence: str) -> str:
    return str(Seq(sequence).reverse_complement()).upper()


def normalize_dna(sequence: str, *, name: str = "sequence") -> str:
    normalized = "".join(sequence.split()).upper().replace("U", "T")
    invalid = sorted(set(normalized) - set("ACGT"))
    if not normalized or invalid:
        raise ValueError(f"{name} must contain only A/C/G/T; invalid={invalid}")
    return normalized


def circular_occurrences(sequence: str, query: str, circular: bool) -> list[int]:
    """Return zero-based starts of exact query occurrences, including origin spans."""
    sequence = normalize_dna(sequence)
    query = normalize_dna(query, name="query")
    if len(query) > len(sequence):
        return []
    haystack = sequence + (sequence[: len(query) - 1] if circular else "")
    return [i for i in range(len(sequence)) if haystack.startswith(query, i)]


def rotation_equivalent(left: str, right: str) -> bool:
    left = normalize_dna(left, name="left")
    right = normalize_dna(right, name="right")
    return len(left) == len(right) and right in left + left


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(normalize_dna(sequence).encode()).hexdigest()


def jsonable(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, JsonModel):
        return value.as_dict()
    raise TypeError(f"cannot JSON-encode {type(value).__name__}")
