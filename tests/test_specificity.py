from __future__ import annotations

import random

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from cloningkit.io import reverse_complement
from cloningkit.specificity import analyze_template_specificity, simulate_product


def test_specificity_reports_multiple_exact_sites() -> None:
    motif = "ACGTACGTACGTACGTACGT"
    sequence = motif + "G" * 30 + motif + "C" * 30
    record = SeqRecord(Seq(sequence), id="repeated")
    record.annotations = {"molecule_type": "DNA", "topology": "circular"}
    reverse_core = reverse_complement("C" * 20)
    result = analyze_template_specificity(record, motif, reverse_core, limit=13)
    assert len(result["forward"]["exact_core_sites_1based"]) == 2
    assert result["scope_warning"].startswith("exact supplied-template")


def test_adjacent_inverse_pcr_product_matches_expected(circular_tadde_record) -> None:
    sequence = str(circular_tadde_record.seq)
    codon_start = 50 + (142 - 1) * 3
    codon_end = codon_start + 3
    expected = sequence[:codon_start] + "CTG" + sequence[codon_end:]
    forward_core = sequence[codon_end : codon_end + 22]
    upstream = sequence[codon_start - 22 : codon_start]
    reverse_core = reverse_complement(upstream)
    result, product = simulate_product(
        circular_tadde_record,
        forward_core,
        "CAG" + reverse_core,
        forward_core_length=22,
        reverse_core_length=22,
        circularize=True,
        expected_sequence=expected,
    )
    assert product.circular
    assert len(product) == len(circular_tadde_record)
    assert result["expected_sequence_match"] is True
    assert result["primer_roles"] == {"forward": "forward", "reverse": "reverse"}


def test_invalid_junction_is_detected(circular_tadde_record) -> None:
    sequence = str(circular_tadde_record.seq)
    codon_start = 50 + (142 - 1) * 3
    codon_end = codon_start + 3
    expected = sequence[:codon_start] + "CTG" + sequence[codon_end:]
    forward_core = sequence[codon_end : codon_end + 22]
    reverse_core = reverse_complement(sequence[codon_start - 22 : codon_start])
    result, _ = simulate_product(
        circular_tadde_record,
        forward_core,
        "CAA" + reverse_core,
        forward_core_length=22,
        reverse_core_length=22,
        circularize=True,
        expected_sequence=expected,
    )
    assert result["expected_sequence_match"] is False


def test_circular_product_lengths_are_rotation_invariant() -> None:
    rng = random.Random(99)
    sequence = "".join(rng.choice("ACGT") for _ in range(220))
    forward = "ACGTTGCAAGTCCGATCGTA"
    reverse_binding = "TTGACCGTACGATGCTAGCA"
    sequence = (
        sequence[:30]
        + forward
        + sequence[50:90]
        + reverse_binding
        + sequence[110:160]
        + reverse_binding
        + sequence[180:]
    )
    reverse = reverse_complement(reverse_binding)

    lengths_by_rotation = []
    for rotation in (0, 100):
        rotated = sequence[rotation:] + sequence[:rotation]
        record = SeqRecord(Seq(rotated), id=f"rotated-{rotation}")
        record.annotations = {"molecule_type": "DNA", "topology": "circular"}
        result = analyze_template_specificity(record, forward, reverse, limit=13)
        lengths_by_rotation.append(sorted(result["pydna"]["product_lengths_bp"]))
        assert {
            (product["forward_primer_role"], product["reverse_primer_role"])
            for product in result["pydna"]["products"]
        } == {("forward", "reverse")}

    assert lengths_by_rotation == [[80, 150], [80, 150]]


def test_verify_product_rejects_unique_homoprimer_product() -> None:
    import pytest

    forward = "ACGTTGCAAGTCCGATCGTA"
    reverse = "TTTTGGGGAAAACCCCTTTG"
    sequence = "G" * 10 + forward + "G" * 70 + reverse_complement(forward) + "G" * 30
    record = SeqRecord(Seq(sequence), id="homoprimer-template")
    record.annotations = {"molecule_type": "DNA"}

    specificity = analyze_template_specificity(record, forward, reverse, limit=13)
    assert specificity["pydna"]["products"] == [
        {
            "length_bp": 110,
            "forward_primer_role": "forward",
            "reverse_primer_role": "forward",
        }
    ]
    with pytest.raises(ValueError, match="does not use the named forward and reverse"):
        simulate_product(record, forward, reverse, limit=13)
