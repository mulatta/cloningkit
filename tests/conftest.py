from __future__ import annotations

import random

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord


@pytest.fixture
def circular_tadde_record() -> SeqRecord:
    rng = random.Random(142)
    sequence = "".join(rng.choice("ACGT") for _ in range(700))
    feature_start = 50
    feature_end = 500
    codon_start = feature_start + (142 - 1) * 3
    sequence = sequence[:codon_start] + "GCC" + sequence[codon_start + 3 :]
    record = SeqRecord(
        Seq(sequence), id="synthetic-plasmid", description="test fixture"
    )
    record.annotations = {"molecule_type": "DNA", "topology": "circular"}
    record.features = [
        SeqFeature(
            SimpleLocation(feature_start, feature_end, strand=1),
            type="CDS",
            qualifiers={"label": ["TadDE"], "transl_table": ["1"]},
        )
    ]
    return record


@pytest.fixture
def origin_target_record() -> SeqRecord:
    rng = random.Random(7)
    sequence = "ATGGCC" + "".join(rng.choice("ACGT") for _ in range(594))
    record = SeqRecord(Seq(sequence), id="origin-plasmid")
    record.annotations = {"molecule_type": "DNA", "topology": "circular"}
    record.features = [
        SeqFeature(
            SimpleLocation(0, 450, strand=1),
            type="CDS",
            qualifiers={"label": ["origin-gene"], "transl_table": ["1"]},
        )
    ]
    return record


@pytest.fixture
def reverse_tadde_record() -> SeqRecord:
    rng = random.Random(314)
    sequence = "".join(rng.choice("ACGT") for _ in range(700))
    feature_start = 50
    feature_end = 500
    genomic_end = feature_end - (142 - 1) * 3
    genomic_start = genomic_end - 3
    sequence = sequence[:genomic_start] + "GGC" + sequence[genomic_end:]
    record = SeqRecord(Seq(sequence), id="reverse-plasmid")
    record.annotations = {"molecule_type": "DNA", "topology": "circular"}
    record.features = [
        SeqFeature(
            SimpleLocation(feature_start, feature_end, strand=-1),
            type="CDS",
            qualifiers={"label": ["reverse-TadDE"], "transl_table": ["1"]},
        )
    ]
    return record
