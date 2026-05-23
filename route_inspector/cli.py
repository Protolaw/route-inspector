from __future__ import annotations

import argparse
import sys
from pathlib import Path

from route_inspector.alchemical_rules import alchemical, classify_alchemical
from route_inspector.alchemical_rules import unwrap_alchemical
from route_inspector.composite_rules import extract
from route_inspector.composite_rules import unwrap as unwrap_composite
from route_inspector.io import (
    read_json,
    resolve_existing_path,
    setup_runtime_cache_dirs,
    write_standard_sidecars,
)
from route_inspector.protection.analysis import (
    ProtectionAnalysisConfig,
    analyze_protection_in_routes,
    load_composite_rule_index,
)
from route_inspector.protection.chython_rules import load_chython_protection_rules
from route_inspector.protection.outputs import (
    dataset_prefix_from_routes_path,
    write_protection_outputs,
)
from route_inspector.scoring import overlap


def add_rule_extraction_arguments(parser: argparse.ArgumentParser) -> None:
    """Add rule extraction options to an argparse parser.

    The parser mutation is kept in one place so terminal commands, tests, and notebook
    examples expose the same options and defaults.
    """
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


def add_parallel_arguments(parser: argparse.ArgumentParser) -> None:
    """Add parallel options to an argparse parser.

    The parser mutation is kept in one place so terminal commands, tests, and notebook
    examples expose the same options and defaults.
    """
    parser.add_argument(
        "--n_cpu",
        "--n-cpu",
        type=int,
        default=1,
        dest="n_cpu",
        help="Number of worker processes to use. Use 1 for sequential execution.",
    )


def add_composite_extraction_arguments(parser: argparse.ArgumentParser) -> None:
    """Add composite extraction options to an argparse parser.

    The parser mutation is kept in one place so terminal commands, tests, and notebook
    examples expose the same options and defaults.
    """
    parser.add_argument("--routes-json", type=Path, required=True)
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Output prefix/path. Separate files are written as "
            "<prefix>_t<size>_composite_rules.tsv."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Directory where standard <dataset>_t<size> rule files and sidecars "
            "are written. The dataset is inferred from --routes-json."
        ),
    )
    add_rule_extraction_arguments(parser)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--min-length", type=int, default=2)
    parser.add_argument("--max-length", type=int, default=5)
    parser.add_argument(
        "--routes-without-composites-output",
        type=Path,
        default=None,
        help=(
            "Optional path for routes that were processed successfully but "
            "did not produce any composite rules, as a get_route_svg_from_json-"
            "readable JSON lookup. Defaults to "
            "<output-prefix>_routes_without_composite_rules.json."
        ),
    )
    parser.add_argument(
        "--skip-routes-without-composites-output",
        action="store_true",
        help=(
            "Do not write the routes_without_composite_rules JSON sidecar. "
            "This is much faster and lighter for very large route collections; "
            "summary counts and reasons are still reported."
        ),
    )
    parser.add_argument(
        "--unique-reactions-first",
        action="store_true",
        help=(
            "Use a two-pass extraction plan: scan routes for unique normalized "
            "reaction SMILES, extract each unique reaction rule once, then "
            "compose route composite rules from the precomputed rule cache."
        ),
    )
    parser.add_argument("--ignore-errors", action="store_true")
    parser.add_argument("--progress-interval", type=int, default=250)
    add_parallel_arguments(parser)
    parser.add_argument(
        "--worker-chunksize",
        type=int,
        default=16,
        help=(
            "Number of routes sent to each worker task during parallel composite "
            "extraction. Larger values reduce multiprocessing overhead."
        ),
    )
    parser.add_argument(
        "--max-pending-chunks",
        type=int,
        default=None,
        help=(
            "Maximum queued worker chunks for parallel composite extraction. "
            "Defaults to 2 * n_cpu."
        ),
    )


def add_alchemical_extraction_arguments(parser: argparse.ArgumentParser) -> None:
    """Add alchemical extraction options to an argparse parser.

    The parser mutation is kept in one place so terminal commands, tests, and notebook
    examples expose the same options and defaults.
    """
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
        default=None,
        help=(
            "Output TSV file or output directory. When a directory is given, "
            "<prefix>_alchemical_rules.tsv and sidecar files are written there."
        ),
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--output-smi", "--output_smi", type=Path, default=None)
    parser.add_argument("--summary", type=Path, default=None)
    parser.add_argument("--errors", type=Path, default=None)
    add_rule_extraction_arguments(parser)
    parser.add_argument("--limit-rows", type=int, default=None)
    parser.add_argument("--limit-applications", type=int, default=None)
    parser.add_argument("--ignore-errors", action="store_true")
    parser.add_argument("--progress-interval", type=int, default=250)
    add_parallel_arguments(parser)


def add_alchemical_classification_arguments(parser: argparse.ArgumentParser) -> None:
    """Add alchemical classification options to an argparse parser.

    The parser mutation is kept in one place so terminal commands, tests, and notebook
    examples expose the same options and defaults.
    """
    parser.add_argument("--alchemical-rules-tsv", type=Path, required=True)
    parser.add_argument("--default-rules-tsv", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--summary", type=Path, default=None)
    add_parallel_arguments(parser)


def add_composite_unwrap_arguments(parser: argparse.ArgumentParser) -> None:
    """Add composite unwrap options to an argparse parser.

    The parser mutation is kept in one place so terminal commands, tests, and notebook
    examples expose the same options and defaults.
    """
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
    """Add alchemical unwrap options to an argparse parser.

    The parser mutation is kept in one place so terminal commands, tests, and notebook
    examples expose the same options and defaults.
    """
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
    """Add scoring options to an argparse parser.

    The parser mutation is kept in one place so terminal commands, tests, and notebook
    examples expose the same options and defaults.
    """
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
    parser.add_argument(
        "--classification-tsv",
        type=Path,
        nargs="+",
        default=None,
        help=(
            "Optional classified alchemical rule TSV files or directories. "
            "Their Composite_rules links are used to compute pos_overlap and "
            "neg_overlap."
        ),
    )
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    add_rule_extraction_arguments(parser)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--min-length", type=int, default=2)
    parser.add_argument("--max-length", type=int, default=5)
    parser.add_argument("--ignore-errors", action="store_true")
    parser.add_argument("--progress-interval", type=int, default=250)
    add_parallel_arguments(parser)


def add_protection_arguments(parser: argparse.ArgumentParser) -> None:
    """Add protection options to an argparse parser.

    The parser mutation is kept in one place so terminal commands, tests, and notebook
    examples expose the same options and defaults.
    """
    parser.add_argument("--routes-json", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--composite-rule-tsv",
        type=Path,
        nargs="+",
        default=None,
        help="Optional composite rule TSV files or directories.",
    )
    parser.add_argument("--alchemical-rules-tsv", type=Path, nargs="+", default=None)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--route-ids-file", type=Path, default=None)
    parser.add_argument("--limit", "--max-routes", type=int, default=None)
    parser.add_argument("--min-composite-size", type=int, default=None)
    parser.add_argument("--max-composite-size", type=int, default=None)
    parser.add_argument("--similarity-threshold", type=float, default=None)
    parser.add_argument("--include-multicenter", action="store_true")
    parser.add_argument("--deprotection-first", action="store_true")
    parser.add_argument("--querycgr-compare", action="store_true")
    parser.add_argument("--write-debug-json", action="store_true")
    parser.add_argument("--write-debug-svg", action="store_true")
    parser.add_argument("--ignore-errors", action="store_true")
    parser.add_argument(
        "--progress-interval",
        type=int,
        default=100,
        help="Print route-level progress every N processed routes. Use 0 to disable.",
    )
    add_parallel_arguments(parser)


def add_preprocess_routes_arguments(parser: argparse.ArgumentParser) -> None:
    """Add PaRoutes preprocessing options to an argparse parser."""
    parser.add_argument("--input-dir", type=Path, default=Path("data/raw"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/clean"))
    parser.add_argument(
        "--summary-dir",
        type=Path,
        default=None,
        help=(
            "Optional output root for generated preprocessing reports. When set, "
            "files are written under <summary-dir>/<dataset>/00_preprocess/."
        ),
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=["n1_routes.json", "n5_routes.json"],
        help=(
            "Dataset filenames to preprocess. Missing underscore names are resolved "
            "against the PaRoutes hyphenated filenames when needed."
        ),
    )
    add_rule_extraction_arguments(parser)
    parser.add_argument("--protection-config", type=Path, default=None)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--ignore-errors", action="store_true")
    parser.add_argument(
        "--progress-interval",
        type=int,
        default=100,
        help="Print route-level progress every N processed routes. Use 0 to disable.",
    )
    add_parallel_arguments(parser)


def read_protection_route_ids(path: Path | None) -> set[str] | None:
    """Read route IDs selected for protection analysis from a text file.

    The function translates CLI arguments into the typed inputs expected by the
    underlying analysis module.
    """
    if path is None:
        return None
    route_ids = set()
    with path.open(encoding="utf-8") as file:
        for line in file:
            value = line.strip()
            if value:
                route_ids.add(value)
    return route_ids


def run_protection_analysis(args: argparse.Namespace) -> int:
    """Run protection analysis using configured inputs.

    The function translates CLI arguments into the typed inputs expected by the
    underlying analysis module.
    """
    setup_runtime_cache_dirs()
    routes_path = resolve_existing_path(args.routes_json)
    config_path = resolve_existing_path(args.config) if args.config else None
    config = ProtectionAnalysisConfig.from_yaml(config_path).with_cli_overrides(args)

    print(f"[analyze-protection] loading routes: {routes_path}", file=sys.stderr, flush=True)
    routes_json = read_json(routes_path)

    print("[analyze-protection] loading composite rule index", file=sys.stderr, flush=True)
    composite_rule_index = load_composite_rule_index(args.composite_rule_tsv)
    print(
        "[analyze-protection] composite rule families: "
        f"{len(composite_rule_index)}",
        file=sys.stderr,
        flush=True,
    )

    print("[analyze-protection] loading chython protection rules", file=sys.stderr, flush=True)
    protection_rules = load_chython_protection_rules()
    source_counts: dict[str, int] = {}
    for rule in protection_rules.values():
        source_counts[rule.source] = source_counts.get(rule.source, 0) + 1
    print(
        "[analyze-protection] protection rules: "
        f"{len(protection_rules)} from {source_counts}",
        file=sys.stderr,
        flush=True,
    )

    route_ids_path = (
        resolve_existing_path(args.route_ids_file) if args.route_ids_file else None
    )
    result = analyze_protection_in_routes(
        routes_json,
        composite_rule_index=composite_rule_index,
        config=config,
        protection_rules=protection_rules,
        limit=args.limit,
        route_ids=read_protection_route_ids(route_ids_path),
        progress_interval=args.progress_interval,
        n_cpu=args.n_cpu,
    )
    output_info = write_protection_outputs(
        result,
        args.output_dir,
        dataset_prefix=dataset_prefix_from_routes_path(routes_path),
    )
    protection_summary = dict(result.summary)
    protection_summary["output_files"] = output_info["output_files"]
    write_standard_sidecars(
        args.output_dir,
        command_name="analyze-protection",
        summary=protection_summary,
        errors=result.summary.get("errors", []),
        input_files=[
            routes_path,
            *(args.composite_rule_tsv or []),
            *(args.alchemical_rules_tsv or []),
        ],
        output_files=output_info["output_files"],
        config_path=config_path,
        cli_args=args,
    )
    print(
        "[analyze-protection] done: "
        f"{result.summary['n_routes']} routes, "
        f"{result.summary['n_deprotection_events']} deprotection events",
        file=sys.stderr,
        flush=True,
    )
    print(
        f"[analyze-protection] summary: {output_info['output_files']['summary']}",
        file=sys.stderr,
        flush=True,
    )
    return 0


def run_preprocess_routes(args: argparse.Namespace) -> int:
    """Run PaRoutes preprocessing using configured inputs."""
    from route_inspector import preprocess_routes

    return preprocess_routes.run(args)


def build_parser() -> argparse.ArgumentParser:
    """Build parser from normalized inputs.

    The returned parser wires every route-inspector subcommand to its implementation
    without importing heavy chemistry modules during help output.
    """
    parser = argparse.ArgumentParser(prog="route-inspector")
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

    protection_parser = subparsers.add_parser(
        "analyze-protection",
        aliases=["analyze-protecting-groups"],
        help="Analyze route-level protection/deprotection strategies.",
    )
    add_protection_arguments(protection_parser)
    protection_parser.set_defaults(func=run_protection_analysis)

    preprocess_parser = subparsers.add_parser(
        "preprocess-routes",
        help="Normalize PaRoutes datasets and split protection-related multicenter reactions.",
    )
    add_preprocess_routes_arguments(preprocess_parser)
    preprocess_parser.set_defaults(func=run_preprocess_routes)

    return parser


def main(argv: list[str] | None = None) -> int:
    """Dispatch the route-inspector CLI to the selected subcommand.

    This is the package entry point used by `python -m route_inspector.cli` and by the
    installed console script.
    """
    args = build_parser().parse_args(argv)
    return args.func(args)


def extract_composite_rules(argv: list[str] | None = None) -> int:
    """Extract composite rules from mapped route data.

    This wrapper lets notebooks and external Python code invoke the same subcommand
    implementation without re-creating argparse setup manually.
    """
    parser = argparse.ArgumentParser(
        description="Extract SynPlanner/chython composite rules from route trees."
    )
    add_composite_extraction_arguments(parser)
    return extract.run(parser.parse_args(argv))


def extract_alchemical_rules(argv: list[str] | None = None) -> int:
    """Extract alchemical rules from mapped route data.

    This wrapper lets notebooks and external Python code invoke the same subcommand
    implementation without re-creating argparse setup manually.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Collapse composite-rule unwrappings into pseudo-reactions and "
            "extract alchemical rules."
        )
    )
    add_alchemical_extraction_arguments(parser)
    return alchemical.run(parser.parse_args(argv))


def classify_alchemical_rules(argv: list[str] | None = None) -> int:
    """Classify alchemical rules against reference rules.

    This wrapper lets notebooks and external Python code invoke the same subcommand
    implementation without re-creating argparse setup manually.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Classify alchemical rules as negative if their QueryCGR matches a "
            "default non-alchemical rule, otherwise positive."
        )
    )
    add_alchemical_classification_arguments(parser)
    return classify_alchemical.run(parser.parse_args(argv))


def score_composite_overlap(argv: list[str] | None = None) -> int:
    """Score composite overlap against reference routes.

    This wrapper lets notebooks and external Python code invoke the same subcommand
    implementation without re-creating argparse setup manually.
    """
    parser = argparse.ArgumentParser(
        description="Score extracted composite rules against reference route JSON."
    )
    add_scoring_arguments(parser)
    return overlap.run(parser.parse_args(argv))


def unwrap_composite_rule(argv: list[str] | None = None) -> int:
    """Unwrap composite rule into a retrosynthetic route.

    This wrapper lets notebooks and external Python code invoke the same subcommand
    implementation without re-creating argparse setup manually.
    """
    parser = argparse.ArgumentParser(
        description="Sequentially apply a composite rule to unwrap a target molecule."
    )
    add_composite_unwrap_arguments(parser)
    return unwrap_composite.run(parser.parse_args(argv))


def unwrap_alchemical_rule(argv: list[str] | None = None) -> int:
    """Unwrap alchemical rule into a retrosynthetic route.

    This wrapper lets notebooks and external Python code invoke the same subcommand
    implementation without re-creating argparse setup manually.
    """
    parser = argparse.ArgumentParser(
        description="Apply one alchemical rule to a target molecule."
    )
    add_alchemical_unwrap_arguments(parser)
    return unwrap_alchemical.run(parser.parse_args(argv))


def analyze_protection(argv: list[str] | None = None) -> int:
    """Run protection analysis from the CLI wrapper.

    This wrapper lets notebooks and external Python code invoke the same subcommand
    implementation without re-creating argparse setup manually.
    """
    parser = argparse.ArgumentParser(
        description="Analyze route-level protection/deprotection strategies."
    )
    add_protection_arguments(parser)
    return run_protection_analysis(parser.parse_args(argv))


def preprocess_routes(argv: list[str] | None = None) -> int:
    """Preprocess PaRoutes datasets from the CLI wrapper."""
    parser = argparse.ArgumentParser(
        description=(
            "Normalize PaRoutes trees and split protection-related multicenter "
            "reactions."
        )
    )
    add_preprocess_routes_arguments(parser)
    return run_preprocess_routes(parser.parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
