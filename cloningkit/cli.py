# Copyright (c) 2026 Seungwon Lee
# SPDX-License-Identifier: MIT
"""Machine-readable command-line interface for Primer Workbench."""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any, NoReturn

from Bio import SeqIO

from . import __version__
from .api import (
    analyze_specificity,
    design_mutagenesis,
    design_primers,
    validate_pair,
    verify_product,
)
from .conditions import ReactionConditions
from .io import load_record, normalize_dna
from .models import (
    DesignConfig,
    MutationTarget,
    Oligo,
    PrimerPair,
    ProductVerificationConfig,
    SpecificityConfig,
)
from .mutagenesis import MutagenesisConfig, build_mutant_record, find_feature

SCHEMA_VERSION = "cloningkit/v1"


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def _finite_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise argparse.ArgumentTypeError("must be finite")
    return parsed


def _nonnegative_float(value: str) -> float:
    parsed = _finite_float(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be nonnegative")
    return parsed


def _positive_float(value: str) -> float:
    parsed = _finite_float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def _temperature_float(value: str) -> float:
    parsed = _finite_float(value)
    if not 0 <= parsed <= 100:
        raise argparse.ArgumentTypeError("must be between 0 and 100 C")
    return parsed


def _conditions(namespace: argparse.Namespace) -> ReactionConditions:
    return ReactionConditions(
        mv_conc=namespace.mv_conc,
        dv_conc=namespace.dv_conc,
        dntp_conc=namespace.dntp_conc,
        dna_conc=namespace.dna_conc,
        annealing_temp_c=namespace.annealing_temp,
        tm_method=namespace.tm_method,
        salt_corrections_method=namespace.salt_correction,
    )


def _add_conditions(
    parser: argparse.ArgumentParser,
    *,
    annealing_help: str = "structure-evaluation temperature, C (default: 60)",
) -> None:
    group = parser.add_argument_group("reaction conditions")
    group.add_argument(
        "--mv-conc",
        type=_nonnegative_float,
        default=50.0,
        help="monovalent salt, mM (default: 50)",
    )
    group.add_argument(
        "--dv-conc",
        type=_nonnegative_float,
        default=1.5,
        help="divalent salt, mM (default: 1.5)",
    )
    group.add_argument(
        "--dntp-conc",
        type=_nonnegative_float,
        default=0.6,
        help="total dNTP, mM (default: 0.6)",
    )
    group.add_argument(
        "--dna-conc",
        type=_positive_float,
        default=50.0,
        help="oligo concentration, nM (default: 50)",
    )
    group.add_argument(
        "--annealing-temp",
        type=_temperature_float,
        default=60.0,
        help=annealing_help,
    )
    group.add_argument(
        "--tm-method", choices=("breslauer", "santalucia"), default="santalucia"
    )
    group.add_argument(
        "--salt-correction",
        choices=("schildkraut", "santalucia", "owczarzy"),
        default="santalucia",
    )


def _add_pair(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--forward", required=True, help="forward full oligo, 5' to 3'")
    parser.add_argument("--reverse", required=True, help="reverse full oligo, 5' to 3'")
    parser.add_argument(
        "--forward-core-length",
        type=_positive_int,
        help="3'-terminal template-binding length; default is full oligo",
    )
    parser.add_argument(
        "--reverse-core-length",
        type=_positive_int,
        help="3'-terminal template-binding length; default is full oligo",
    )


def _primer_pair(
    namespace: argparse.Namespace, *, thermodynamic: bool = False
) -> PrimerPair:
    forward = normalize_dna(namespace.forward, name="forward primer")
    reverse = normalize_dna(namespace.reverse, name="reverse primer")
    # Sequence and length checks stay ahead of core parsing in the v1 CLI contract.
    if thermodynamic and (len(forward) > 60 or len(reverse) > 60):
        raise ValueError(
            "Primer3 thermodynamic structure calculations require each full "
            "oligo to be at most 60 nt"
        )
    return PrimerPair(
        forward=Oligo(forward, namespace.forward_core_length),
        reverse=Oligo(reverse, namespace.reverse_core_length),
    )


def _cmd_design(namespace: argparse.Namespace) -> dict[str, Any]:
    record = load_record(namespace.input)
    conditions = _conditions(namespace)
    # Keep template errors ahead of design-bound errors in the CLI domain contract.
    normalize_dna(str(record.seq), name="template")
    config = DesignConfig(
        target_start_1based=namespace.target_start,
        target_length=namespace.target_length,
        product_min_bp=namespace.product_min,
        product_max_bp=namespace.product_max,
        min_size_nt=namespace.min_size,
        opt_size_nt=namespace.opt_size,
        max_size_nt=namespace.max_size,
        min_tm_c=namespace.min_tm,
        opt_tm_c=namespace.opt_tm,
        max_tm_c=namespace.max_tm,
        min_gc_percent=namespace.min_gc,
        max_gc_percent=namespace.max_gc,
        max_delta_tm_c=namespace.max_delta_tm,
        num_return=namespace.num_return,
        conditions=conditions,
    )
    return design_primers(record, config).as_dict()


def _cmd_validate(namespace: argparse.Namespace) -> dict[str, Any]:
    conditions = _conditions(namespace)
    return validate_pair(
        _primer_pair(namespace, thermodynamic=True), conditions
    ).as_dict()


def _cmd_specificity(namespace: argparse.Namespace) -> dict[str, Any]:
    return analyze_specificity(
        load_record(namespace.input),
        _primer_pair(namespace),
        SpecificityConfig(pydna_limit=namespace.limit),
    ).as_dict()


def _cmd_verify(namespace: argparse.Namespace) -> dict[str, Any]:
    record = load_record(namespace.input)
    expected_sequence = None
    if namespace.expected:
        expected_sequence = str(load_record(namespace.expected).seq)
    verification = verify_product(
        record,
        _primer_pair(namespace),
        ProductVerificationConfig(
            pydna_limit=namespace.limit,
            circularize=namespace.circularize,
            expected_sequence=expected_sequence,
        ),
    )
    result = verification.as_dict()
    product = verification.product
    if namespace.output_product:
        output = Path(namespace.output_product)
        output.parent.mkdir(parents=True, exist_ok=True)
        fmt = (
            "genbank"
            if output.suffix.lower() in {".gb", ".gbk", ".genbank"}
            else "fasta"
        )
        product.annotations.setdefault("molecule_type", "DNA")
        if product.circular:
            product.annotations["topology"] = "circular"
        SeqIO.write(product, output, fmt)
        result["output_product"] = str(output)
    return result


def _cmd_mutagenesis(namespace: argparse.Namespace) -> dict[str, Any]:
    if namespace.output_mutant and len(namespace.codon or ()) != 1:
        raise ValueError("--output-mutant requires exactly one explicit --codon")
    conditions = _conditions(namespace)
    config = MutagenesisConfig(
        min_core_length=namespace.min_core_length,
        max_core_length=namespace.max_core_length,
        min_tm_c=namespace.min_tm,
        max_tm_c=namespace.max_tm,
        target_tm_c=namespace.target_tm,
        max_delta_tm_c=namespace.max_delta_tm,
        min_gc_percent=namespace.min_gc,
        max_gc_percent=namespace.max_gc,
        max_homopolymer=namespace.max_homopolymer,
        max_three_prime_gc=namespace.max_three_prime_gc,
        max_hairpin_risk_kcal_mol=namespace.max_hairpin_risk,
        max_homodimer_risk_kcal_mol=namespace.max_homodimer_risk,
        max_heterodimer_risk_kcal_mol=namespace.max_heterodimer_risk,
        max_end_dimer_risk_kcal_mol=namespace.max_end_dimer_risk,
        codon_table=namespace.codon_table,
        codons=tuple(namespace.codon or ()),
        tail_sides=tuple(namespace.tail_side),
        pydna_limit=namespace.limit,
        verify_top=namespace.verify_top,
        include_failed=namespace.include_failed,
        conditions=conditions,
    )
    record = load_record(namespace.input)
    target = MutationTarget(
        feature_name=namespace.feature,
        aa_position_1based=namespace.aa_position,
        target_amino_acid=namespace.target_aa,
    )
    typed_result = design_mutagenesis(record, target, config)
    result = typed_result.as_dict()
    if namespace.output_mutant:
        codons = typed_result.target.candidate_coding_codons
        if len(codons) != 1:
            raise AssertionError("explicit codon did not define one candidate codon")
        feature = find_feature(record, namespace.feature)
        mutant = build_mutant_record(record, feature, namespace.aa_position, codons[0])
        output = Path(namespace.output_mutant)
        output.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(mutant, output, "genbank")
        result["output_mutant"] = str(output)
    return result


def _add_output_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--output-json",
        type=Path,
        help="write the JSON envelope to this file instead of stdout",
    )
    parser.add_argument("--pretty", action="store_true", help="indent JSON output")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cloningkit",
        description="Reproducible primer design and sequence-level verification",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    design = subparsers.add_parser(
        "design", help="design conventional inward-facing PCR pairs with Primer3"
    )
    design.add_argument(
        "--input",
        required=True,
        type=Path,
        help="single FASTA/GenBank template (linear)",
    )
    design.add_argument(
        "--target-start", type=_positive_int, help="target start, 1-based inclusive"
    )
    design.add_argument("--target-length", type=_positive_int)
    design.add_argument("--product-min", type=_positive_int, default=100)
    design.add_argument("--product-max", type=_positive_int, default=1000)
    design.add_argument("--min-size", type=_positive_int, default=18)
    design.add_argument("--opt-size", type=_positive_int, default=20)
    design.add_argument("--max-size", type=_positive_int, default=25)
    design.add_argument("--min-tm", type=_finite_float, default=58.0)
    design.add_argument("--opt-tm", type=_finite_float, default=60.0)
    design.add_argument("--max-tm", type=_finite_float, default=62.0)
    design.add_argument("--min-gc", type=_finite_float, default=40.0)
    design.add_argument("--max-gc", type=_finite_float, default=60.0)
    design.add_argument("--max-delta-tm", type=_nonnegative_float, default=2.0)
    design.add_argument("--num-return", type=_positive_int, default=5)
    _add_conditions(
        design,
        annealing_help=(
            "Primer3 primer-binding annealing temperature, C (default: 60); "
            "incompatible with Owczarzy salt correction"
        ),
    )
    _add_output_args(design)
    design.set_defaults(func=_cmd_design)

    validate = subparsers.add_parser(
        "validate", help="calculate core Tm and full-oligo structures with primer3-py"
    )
    _add_pair(validate)
    _add_conditions(validate)
    _add_output_args(validate)
    validate.set_defaults(func=_cmd_validate)

    specificity = subparsers.add_parser(
        "specificity",
        help="find exact sites and enumerate pydna PCR products on one template",
    )
    specificity.add_argument("--input", required=True, type=Path)
    _add_pair(specificity)
    specificity.add_argument(
        "--limit",
        type=_positive_int,
        default=13,
        help="pydna minimum 3' annealing length",
    )
    _add_output_args(specificity)
    specificity.set_defaults(func=_cmd_specificity)

    mutagenesis = subparsers.add_parser(
        "mutagenesis", help="exhaustively search adjacent inverse-PCR mutagenesis pairs"
    )
    mutagenesis.add_argument("--input", required=True, type=Path)
    mutagenesis.add_argument(
        "--feature", required=True, help="unique feature label/gene/product/locus_tag"
    )
    mutagenesis.add_argument(
        "--aa-position",
        required=True,
        type=_positive_int,
        help="1-based residue within the feature",
    )
    mutagenesis.add_argument(
        "--target-aa", required=True, help="one-letter target amino acid"
    )
    mutagenesis.add_argument(
        "--codon-table", default="e_coli_316407", help="python-codon-tables table name"
    )
    mutagenesis.add_argument(
        "--codon", action="append", help="restrict target coding codon; repeatable"
    )
    mutagenesis.add_argument(
        "--tail-side",
        action="append",
        choices=("forward", "reverse"),
        default=None,
        help="mutation-tail side; repeat for both (default: both)",
    )
    mutagenesis.add_argument("--min-core-length", type=_positive_int, default=18)
    mutagenesis.add_argument("--max-core-length", type=_positive_int, default=30)
    mutagenesis.add_argument("--min-tm", type=_finite_float, default=55.0)
    mutagenesis.add_argument("--target-tm", type=_finite_float, default=60.0)
    mutagenesis.add_argument("--max-tm", type=_finite_float, default=65.0)
    mutagenesis.add_argument("--max-delta-tm", type=_nonnegative_float, default=2.5)
    mutagenesis.add_argument("--min-gc", type=_finite_float, default=35.0)
    mutagenesis.add_argument("--max-gc", type=_finite_float, default=65.0)
    mutagenesis.add_argument("--max-homopolymer", type=_positive_int, default=5)
    mutagenesis.add_argument("--max-three-prime-gc", type=_positive_int, default=3)
    mutagenesis.add_argument(
        "--max-hairpin-risk",
        type=_nonnegative_float,
        default=3.0,
        help="maximum -dG magnitude at Ta, kcal/mol",
    )
    mutagenesis.add_argument(
        "--max-homodimer-risk",
        type=_nonnegative_float,
        default=9.0,
        help="maximum -dG magnitude at Ta, kcal/mol",
    )
    mutagenesis.add_argument(
        "--max-heterodimer-risk",
        type=_nonnegative_float,
        default=9.0,
        help="maximum -dG magnitude at Ta, kcal/mol",
    )
    mutagenesis.add_argument(
        "--max-end-dimer-risk",
        type=_nonnegative_float,
        default=5.0,
        help="maximum end-stability -dG magnitude, kcal/mol",
    )
    mutagenesis.add_argument("--limit", type=_positive_int, default=13)
    mutagenesis.add_argument(
        "--verify-top",
        type=int,
        default=25,
        help="number of score-sorted Pareto pairs to verify with pydna; may be zero",
    )
    mutagenesis.add_argument(
        "--include-failed",
        action="store_true",
        help="include hard-filter failures in candidates",
    )
    mutagenesis.add_argument(
        "--output-mutant",
        type=Path,
        help="write expected mutant GenBank; requires one --codon",
    )
    _add_conditions(mutagenesis)
    _add_output_args(mutagenesis)
    mutagenesis.set_defaults(func=_cmd_mutagenesis)

    verify = subparsers.add_parser(
        "verify-product",
        help="simulate a unique PCR product and optional blunt circularization",
    )
    verify.add_argument("--input", required=True, type=Path)
    _add_pair(verify)
    verify.add_argument("--limit", type=_positive_int, default=13)
    verify.add_argument("--circularize", action="store_true")
    verify.add_argument("--expected", type=Path, help="expected FASTA/GenBank sequence")
    verify.add_argument(
        "--output-product",
        type=Path,
        help="write product as FASTA or GenBank by extension",
    )
    _add_output_args(verify)
    verify.set_defaults(func=_cmd_verify)
    return parser


def _emit(envelope: dict[str, Any], namespace: argparse.Namespace) -> None:
    serialized = json.dumps(
        envelope,
        indent=2 if namespace.pretty else None,
        sort_keys=True,
        allow_nan=False,
    )
    output = getattr(namespace, "output_json", None)
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized + "\n")
    else:
        print(serialized)


def _fail(namespace: argparse.Namespace, error: Exception) -> NoReturn:
    envelope = {
        "schema_version": SCHEMA_VERSION,
        "command": getattr(namespace, "command", None),
        "ok": False,
        "error": {"type": type(error).__name__, "message": str(error)},
    }
    print(json.dumps(envelope, sort_keys=True), file=sys.stderr)
    raise SystemExit(2)


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_parser()
    namespace = parser.parse_args(argv)
    if namespace.command == "mutagenesis" and namespace.tail_side is None:
        namespace.tail_side = ["forward", "reverse"]
    if namespace.command == "mutagenesis" and namespace.verify_top < 0:
        parser.error("--verify-top must be nonnegative")
    try:
        result = namespace.func(namespace)
        envelope = {
            "schema_version": SCHEMA_VERSION,
            "command": namespace.command,
            "ok": True,
            "result": result,
        }
        _emit(envelope, namespace)
    except (ValueError, OSError, RuntimeError, AssertionError) as error:
        _fail(namespace, error)


if __name__ == "__main__":
    main()
