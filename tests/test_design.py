from __future__ import annotations

import random

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from cloningkit.conditions import ReactionConditions
from cloningkit.design import design_standard_primers


def test_standard_design_returns_ranked_pairs() -> None:
    rng = random.Random(1)
    sequence = "".join(rng.choice("ACGT") for _ in range(1200))
    template = SeqRecord(Seq(sequence), id="linear-template")
    template.annotations["molecule_type"] = "DNA"
    result = design_standard_primers(
        template,
        target_start_1based=500,
        target_length=40,
        product_min=100,
        product_max=300,
        min_gc=20,
        max_gc=80,
        num_return=3,
    )
    assert 1 <= len(result["pairs"]) <= 3
    assert [pair["rank"] for pair in result["pairs"]] == list(
        range(len(result["pairs"]))
    )
    assert result["template"]["circular"] is False


def test_circular_design_warns_about_origin(circular_tadde_record) -> None:
    result = design_standard_primers(
        circular_tadde_record,
        target_start_1based=300,
        target_length=10,
        product_min=80,
        product_max=250,
        min_gc=20,
        max_gc=80,
        num_return=1,
    )
    assert result["template"]["circular"] is True
    assert "origin-spanning" in result["interpretation"]["circular_template"]


def test_design_rejects_owczarzy_with_primer_binding_temperature() -> None:
    template = SeqRecord(Seq("ACGT" * 300), id="linear-template")
    with pytest.raises(ValueError, match="cannot combine PRIMER_ANNEALING_TEMP"):
        design_standard_primers(
            template,
            conditions=ReactionConditions(salt_corrections_method="owczarzy"),
        )
