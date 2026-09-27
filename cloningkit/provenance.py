# Copyright (c) 2026 Seungwon Lee
# SPDX-License-Identifier: MIT
"""Runtime provenance reported with machine-readable results."""

from __future__ import annotations

from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version

import Bio
import primer3
import pydna
import python_codon_tables

from .serialization import JsonModel, JsonObject


@dataclass(frozen=True, slots=True)
class ToolVersions(JsonModel):
    cloningkit_version: str
    python_codon_tables_version: str
    biopython_version: str
    primer3_py_version: str
    pydna_version: str


def runtime_versions() -> ToolVersions:
    """Return typed versions for every engine represented in result payloads."""
    try:
        cloningkit_version = version("cloningkit")
    except PackageNotFoundError:
        cloningkit_version = "unknown"
    return ToolVersions(
        cloningkit_version=cloningkit_version,
        python_codon_tables_version=str(python_codon_tables.__version__),
        biopython_version=str(Bio.__version__),
        primer3_py_version=str(primer3.__version__),
        pydna_version=str(pydna.__version__),
    )


def tool_versions() -> JsonObject:
    """Keep legacy dictionary-producing engines stable during API migration."""
    return runtime_versions().as_dict()
