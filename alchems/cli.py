from __future__ import annotations

import argparse
from pathlib import Path

from alchems.alchemical_rules import alchemical, classify_alchemical
from alchems.alchemical_rules import unwrap_alchemical
from alchems.composite_rules import extract
from alchems.composite_rules import unwrap as unwrap_composite
from alchems.scoring import overlap


def add_rule_extraction_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--environment-atom-count", type=int, default=1)
    parser.add_argument("--include-rings", action="store_true")
    parser.add_argument("--keep-leaving-groups", action="store_true", default=True)
    parser.add_argument(
        "--drop-leaving-groups",
        dest="keep_leaving_groups",
        action="store_false",
    )
    parser.add_argument("--keep-incoming-groups", action="store_true")
    parser.add_argument("--reactor-validation", action="store_true")


def add_composite_extraction_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--routes-json", type=Path, required=True)
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help=(
            "Output prefix/path. Separate files are written as "
            "<prefix>_t<size>_composite_rules.tsv."
        ),
    )
    add_rule_extraction_arguments(parser)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--min-length", type=int, default=2)
    parser.add_argument("--max-length", type=int, default=5)
    parser.add_argument("--ignore-errors", action="store_true")
    parser.add_argument("--progress-interval", type=int, default=250)


def add_alchemical_extraction_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--composite-rule-tsv",
        type=Path,
        nargs="+",
        required=True,
        help="One or more composite rule TSV files or directories containing them.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help=(
            "Output TSV file or output directory. When a directory is given, "
            "<prefix>_alchemical_rules.tsv and sidecar files are written there."
        ),
    )
    parser.add_argument("--output-smi", "--output_smi", type=Path, default=None)
    parser.add_argument("--summary", type=Path, default=None)
    parser.add_argument("--errors", type=Path, default=None)
    add_rule_extraction_arguments(parser)
    parser.add_argument("--limit-rows", type=int, default=None)
    parser.add_argument("--limit-applications", type=int, default=None)
    parser.add_argument("--ignore-errors", action="store_true")
    parser.add_argument("--progress-interval", type=int, default=250)


def add_alchemical_classification_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--alchemical-rules-tsv", type=Path, required=True)
    parser.add_argument("--default-rules-tsv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, default=None)


def add_composite_unwrap_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--smiles", required=True, help="Target molecule SMILES.")
    rule_source = parser.add_mutually_exclusive_group(required=True)
    rule_source.add_argument("--composite-rule", help="Composite rule string.")
    rule_source.add_argument(
        "--composite-rule-tsv",
        type=Path,
        help="TSV containing a Composite_rule column.",
    )
    parser.add_argument("--row", type=int, default=0, help="0-based TSV row index.")
    parser.add_argument("--route-id", type=int, default=0)
    parser.add_argument("--output-json", type=Path, default=None)
    parser.add_argument("--output-svg", type=Path, default=None)
    parser.add_argument("--labeled", action="store_true")
    parser.add_argument("--do-not-mark-leaves-in-stock", action="store_true")


def add_alchemical_unwrap_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--smiles", required=True, help="Target molecule SMILES.")
    rule_source = parser.add_mutually_exclusive_group(required=True)
    rule_source.add_argument("--alchemical-rule", help="Alchemical rule SMARTS.")
    rule_source.add_argument(
        "--alchemical-rule-tsv",
        type=Path,
        help="TSV containing an Alchemical_rule column.",
    )
    parser.add_argument("--row", type=int, default=0, help="0-based TSV row index.")
    parser.add_argument("--route-id", type=int, default=0)
    parser.add_argument("--output-json", type=Path, default=None)
    parser.add_argument("--output-svg", type=Path, default=None)
    parser.add_argument("--labeled", action="store_true")
    parser.add_argument("--do-not-mark-leaves-in-stock", action="store_true")


def add_scoring_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--extracted-tsv",
        "--composite-rule-tsv",
        type=Path,
        nargs="+",
        required=True,
        dest="extracted_tsv",
        help="One or more extracted composite rule TSV files or directories.",
    )
    parser.add_argument(
        "--reference-routes-json",
        "--routes-json",
        type=Path,
        required=True,
        dest="reference_routes_json",
        help="Reference route JSON. Composite rules are extracted from this file.",
    )
    parser.add_argument("--output", type=Path, required=True)
    add_rule_extraction_arguments(parser)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--min-length", type=int, default=2)
    parser.add_argument("--max-length", type=int, default=5)
    parser.add_argument("--ignore-errors", action="store_true")
    parser.add_argument("--progress-interval", type=int, default=250)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="alchems")
    subparsers = parser.add_subparsers(dest="command", required=True)

    extract_parser = subparsers.add_parser(
        "extract-composite-rules",
        aliases=["extract-composites"],
        help="Extract composite rules from route JSON.",
    )
    add_composite_extraction_arguments(extract_parser)
    extract_parser.set_defaults(func=extract.run)

    alchemical_parser = subparsers.add_parser(
        "extract-alchemical-rules",
        aliases=["extract-alchemicals"],
        help="Extract alchemical rules from composite rule TSVs.",
    )
    add_alchemical_extraction_arguments(alchemical_parser)
    alchemical_parser.set_defaults(func=alchemical.run)

    classify_parser = subparsers.add_parser(
        "classify-alchemical-rules",
        aliases=["classify-alchemicals"],
        help="Classify alchemical rules as positive or negative.",
    )
    add_alchemical_classification_arguments(classify_parser)
    classify_parser.set_defaults(func=classify_alchemical.run)

    unwrap_composite_parser = subparsers.add_parser(
        "unwrap-composite-rule",
        aliases=["unwrap-composite"],
        help="Apply a composite rule sequence to a target molecule.",
    )
    add_composite_unwrap_arguments(unwrap_composite_parser)
    unwrap_composite_parser.set_defaults(func=unwrap_composite.run)

    unwrap_alchemical_parser = subparsers.add_parser(
        "unwrap-alchemical-rule",
        aliases=["unwrap-alchemical"],
        help="Apply one alchemical rule to a target molecule.",
    )
    add_alchemical_unwrap_arguments(unwrap_alchemical_parser)
    unwrap_alchemical_parser.set_defaults(func=unwrap_alchemical.run)

    score_parser = subparsers.add_parser(
        "score-composite-overlap",
        aliases=["score-overlap"],
        help="Score extracted composite rules against reference route JSON.",
    )
    add_scoring_arguments(score_parser)
    score_parser.set_defaults(func=overlap.run)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


def extract_composite_rules(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Extract SynPlanner/chython composite rules from route trees."
    )
    add_composite_extraction_arguments(parser)
    return extract.run(parser.parse_args(argv))


def extract_alchemical_rules(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Collapse composite-rule unwrappings into pseudo-reactions and "
            "extract alchemical rules."
        )
    )
    add_alchemical_extraction_arguments(parser)
    return alchemical.run(parser.parse_args(argv))


def classify_alchemical_rules(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Classify alchemical rules as negative if their QueryCGR matches a "
            "default non-alchemical rule, otherwise positive."
        )
    )
    add_alchemical_classification_arguments(parser)
    return classify_alchemical.run(parser.parse_args(argv))


def score_composite_overlap(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Score extracted composite rules against reference route JSON."
    )
    add_scoring_arguments(parser)
    return overlap.run(parser.parse_args(argv))


def unwrap_composite_rule(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Sequentially apply a composite rule to unwrap a target molecule."
    )
    add_composite_unwrap_arguments(parser)
    return unwrap_composite.run(parser.parse_args(argv))


def unwrap_alchemical_rule(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Apply one alchemical rule to a target molecule."
    )
    add_alchemical_unwrap_arguments(parser)
    return unwrap_alchemical.run(parser.parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
