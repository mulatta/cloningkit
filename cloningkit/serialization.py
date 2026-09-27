# Copyright (c) 2026 Seungwon Lee
# SPDX-License-Identifier: MIT
"""Lossless JSON serialization shared by typed public models."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, fields

type JsonScalar = str | int | float | bool | None
type JsonValue = JsonScalar | list[JsonValue] | dict[str, JsonValue]
type JsonObject = dict[str, JsonValue]


@dataclass(frozen=True, slots=True)
class JsonModel:
    """Serialize frozen dataclasses without exposing mutable result dictionaries."""

    def as_dict(self) -> JsonObject:
        output: JsonObject = {}
        for item in fields(self):
            if item.metadata.get("serialize", True) is False:
                continue
            value = getattr(self, item.name)
            if value is None and item.metadata.get("omit_none", False) is True:
                continue
            key = str(item.metadata.get("json_name", item.name))
            output[key] = to_jsonable(value)
        return output


def to_jsonable(value: object) -> JsonValue:
    """Return JSON-compatible values while preserving tuple/list field semantics."""
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, JsonModel):
        return value.as_dict()
    if isinstance(value, Mapping):
        return {str(key): to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(item) for item in value]
    raise TypeError(f"cannot JSON-encode {type(value).__name__}")
