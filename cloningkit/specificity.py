# Copyright (c) 2026 Seungwon Lee
# SPDX-License-Identifier: MIT
"""Exact-template specificity and pydna PCR product simulation."""

from __future__ import annotations

import copy
from typing import Any

from pydna.amplify import Anneal
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer

from .io import (
    circular_occurrences,
    is_circular,
    normalize_dna,
    reverse_complement,
    rotation_equivalent,
    sequence_sha256,
)
from .provenance import tool_versions
from .thermo import split_core


def _pydna_template(record: Any) -> Dseqrecord:
    return Dseqrecord(record, circular=is_circular(record))


def _enumerate_pydna_products(
    record: Any,
    forward: str,
    reverse: str,
    *,
    forward_core_length: int,
    reverse_core_length: int,
    limit: int,
) -> tuple[str, list[Any]]:
    """Enumerate primer-mixture products without pydna's circular position mutation.

    pydna 5.5.16 mutates binding positions while iterating multiple circular
    site pairs. Isolating every already-discovered forward/reverse binding pair
    preserves its product model while keeping lengths origin-independent.
    """
    if limit < 1:
        raise ValueError("pydna limit must be positive")
    normalize_dna(str(record.seq), name="template")
    template = _pydna_template(record)
    primers = [
        Primer(
            forward,
            id="forward",
            name="forward",
            description="forward",
            footprint=forward_core_length,
        ),
        Primer(
            reverse,
            id="reverse",
            name="reverse",
            description="reverse",
            footprint=reverse_core_length,
        ),
    ]
    anneal = Anneal(primers, template, limit=limit)
    report = str(anneal.report())
    products: list[Any] = []
    for forward_binding in anneal.forward_primers:
        for reverse_binding in anneal.reverse_primers:
            if (
                not template.circular
                and forward_binding.position > reverse_binding.position
            ):
                continue
            isolated = Anneal([], anneal.template, limit=limit)
            isolated.forward_primers = [copy.deepcopy(forward_binding)]
            isolated.reverse_primers = [copy.deepcopy(reverse_binding)]
            pair_products = list(isolated.products)
            if len(pair_products) != 1:
                raise AssertionError(
                    "isolated pydna binding pair did not produce exactly one product"
                )
            products.extend(pair_products)
    return report, products


def _product_roles(product: Any) -> tuple[str, str]:
    return str(product.forward_primer.id), str(product.reverse_primer.id)


def analyze_template_specificity(
    record: Any,
    forward: str,
    reverse: str,
    *,
    forward_core_length: int | None = None,
    reverse_core_length: int | None = None,
    limit: int = 13,
) -> dict[str, Any]:
    """Inspect exact circular hits and all pydna products on one template."""
    forward = normalize_dna(forward, name="forward primer")
    reverse = normalize_dna(reverse, name="reverse primer")
    f_tail, f_core = split_core(forward, forward_core_length)
    r_tail, r_core = split_core(reverse, reverse_core_length)
    circular = is_circular(record)
    sequence = str(record.seq).upper()
    forward_starts = circular_occurrences(sequence, f_core, circular)
    reverse_binding_sequence = reverse_complement(r_core)
    reverse_starts = circular_occurrences(sequence, reverse_binding_sequence, circular)
    report, products = _enumerate_pydna_products(
        record,
        forward,
        reverse,
        forward_core_length=len(f_core),
        reverse_core_length=len(r_core),
        limit=limit,
    )
    product_details = [
        {
            "length_bp": len(product),
            "forward_primer_role": _product_roles(product)[0],
            "reverse_primer_role": _product_roles(product)[1],
        }
        for product in products
    ]
    return {
        "template": {
            "id": record.id,
            "length_bp": len(record),
            "circular": circular,
            "sha256": sequence_sha256(sequence),
        },
        "forward": {
            "tail_5to3": f_tail,
            "core_5to3": f_core,
            "exact_core_sites_1based": [
                {
                    "start": start + 1,
                    "end": ((start + len(f_core) - 1) % len(record)) + 1,
                }
                for start in forward_starts
            ],
        },
        "reverse": {
            "tail_5to3": r_tail,
            "core_5to3": r_core,
            "plus_strand_binding_sequence_5to3": reverse_binding_sequence,
            "exact_core_sites_1based": [
                {
                    "start": start + 1,
                    "end": ((start + len(reverse_binding_sequence) - 1) % len(record))
                    + 1,
                }
                for start in reverse_starts
            ],
        },
        "pydna": {
            "limit": limit,
            "report": report,
            "product_count": len(products),
            "product_lengths_bp": [len(product) for product in products],
            "products": product_details,
        },
        "scope_warning": (
            "exact supplied-template analysis only; not mismatch-tolerant or genome-wide"
        ),
        "tool": tool_versions(),
    }


def simulate_product(
    record: Any,
    forward: str,
    reverse: str,
    *,
    forward_core_length: int | None = None,
    reverse_core_length: int | None = None,
    limit: int = 13,
    circularize: bool = False,
    expected_sequence: str | None = None,
) -> tuple[dict[str, Any], Any]:
    """Simulate a unique pydna PCR product and optional blunt circularization."""
    forward = normalize_dna(forward, name="forward primer")
    reverse = normalize_dna(reverse, name="reverse primer")
    _, f_core = split_core(forward, forward_core_length)
    _, r_core = split_core(reverse, reverse_core_length)
    report, products = _enumerate_pydna_products(
        record,
        forward,
        reverse,
        forward_core_length=len(f_core),
        reverse_core_length=len(r_core),
        limit=limit,
    )
    if len(products) != 1:
        raise ValueError(
            "expected exactly one pydna PCR product in the primer mixture, "
            f"found {len(products)}"
        )
    product = products[0]
    roles = _product_roles(product)
    if roles != ("forward", "reverse"):
        raise ValueError(
            "unique pydna product does not use the named forward and reverse "
            f"primers in their declared roles; observed={roles}"
        )
    if circularize:
        try:
            final_product = product.looped()
        except TypeError as error:
            raise ValueError(
                f"pydna could not circularize the PCR product: {error}"
            ) from error
    else:
        final_product = product
    final_sequence = str(final_product.seq).upper()
    equivalent = None
    if expected_sequence is not None:
        expected_sequence = normalize_dna(expected_sequence, name="expected sequence")
        equivalent = (
            rotation_equivalent(final_sequence, expected_sequence)
            if circularize
            else final_sequence == expected_sequence
        )
    result = {
        "product_length_bp": len(final_product),
        "circular": bool(final_product.circular),
        "sha256": sequence_sha256(final_sequence),
        "expected_sequence_match": equivalent,
        "figure": product.figure(),
        "annealing_report": report,
        "primer_roles": {
            "forward": roles[0],
            "reverse": roles[1],
        },
        "scope_warning": "sequence/topology simulation, not an experimental yield prediction",
        "tool": tool_versions(),
    }
    return result, final_product
