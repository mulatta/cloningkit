from __future__ import annotations

import json
from pathlib import Path

import primer3
from Bio import SeqIO

from cloningkit.cli import main


def invoke(capsys, arguments: list[str]) -> tuple[int, dict]:
    code = 0
    try:
        main(arguments)
    except SystemExit as error:
        code = error.code if isinstance(error.code, int) else 1
    output = capsys.readouterr()
    payload = output.out if code == 0 else output.err
    return code, json.loads(payload)


def test_validate_cli_emits_versioned_envelope(capsys) -> None:
    code, payload = invoke(
        capsys,
        [
            "validate",
            "--forward",
            "GCACTGCTGTGCGATTTTTATCG",
            "--reverse",
            "CAGACATTCATCGGCCAGAATACCC",
            "--forward-core-length",
            "23",
            "--reverse-core-length",
            "22",
        ],
    )
    assert code == 0
    assert payload["schema_version"] == "cloningkit/v1"
    assert payload["command"] == "validate"
    assert payload["ok"] is True
    assert payload["result"]["tool"]["primer3_py_version"] == primer3.__version__
    assert payload["result"]["tool"]["pydna_version"]


def test_cli_domain_error_is_machine_readable(capsys) -> None:
    code, payload = invoke(
        capsys,
        [
            "validate",
            "--forward",
            "ACGN",
            "--reverse",
            "ACGT",
        ],
    )
    assert code == 2
    assert payload["ok"] is False
    assert payload["error"]["type"] == "ValueError"
    assert payload["error"]["message"] == (
        "forward primer must contain only A/C/G/T; invalid=['N']"
    )


def test_verify_product_cli_writes_genbank(
    capsys, tmp_path: Path, circular_tadde_record
) -> None:
    template = tmp_path / "template.gb"
    SeqIO.write(circular_tadde_record, template, "genbank")
    sequence = str(circular_tadde_record.seq)
    start = 50 + (142 - 1) * 3
    end = start + 3
    forward = sequence[end : end + 20]
    from cloningkit.io import reverse_complement

    reverse = "CAG" + reverse_complement(sequence[start - 20 : start])
    product = tmp_path / "product.gb"
    code, payload = invoke(
        capsys,
        [
            "verify-product",
            "--input",
            str(template),
            "--forward",
            forward,
            "--reverse",
            reverse,
            "--forward-core-length",
            "20",
            "--reverse-core-length",
            "20",
            "--circularize",
            "--output-product",
            str(product),
        ],
    )
    assert code == 0
    assert payload["result"]["circular"] is True
    written = SeqIO.read(product, "genbank")
    assert len(written) == len(circular_tadde_record)
    assert written.annotations["topology"] == "circular"


def test_output_mutant_requires_one_explicit_codon(capsys, tmp_path: Path) -> None:
    code, payload = invoke(
        capsys,
        [
            "mutagenesis",
            "--input",
            str(tmp_path / "not-read.gb"),
            "--feature",
            "gene",
            "--aa-position",
            "1",
            "--target-aa",
            "M",
            "--output-mutant",
            str(tmp_path / "mutant.gb"),
        ],
    )
    assert code == 2
    assert payload["error"]["type"] == "ValueError"
    assert "exactly one explicit --codon" in payload["error"]["message"]
