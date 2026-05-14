from __future__ import annotations

import argparse
import copy
import json
import sys
import traceback
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

PROJECT_ROOT = Path(__file__).resolve().parents[2]

if __package__ in (None, "") and str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from alchems.io import (
    read_json,
    resolve_existing_path,
    setup_runtime_cache_dirs,
    write_composite_errors as write_errors,
    write_composite_rules,
    write_composite_summary as write_summary,
)


@dataclass(frozen=True)
class ReactionRuleStep:
    """A route reaction annotated with its extracted rule and reaction center."""

    rule_smarts: str
    center_atoms: frozenset[int]
    reaction_smiles: str
    target_smiles: str = ""


@dataclass
class RouteProcessingStats:
    routes_seen: int = 0
    routes_with_composites: int = 0
    reactions_seen: int = 0
    reaction_rule_cache_hits: int = 0
    reaction_rule_cache_misses: int = 0
    skipped_reactions: int = 0
    errors: int = 0


class RuleExtractionError(Exception):
    """Raised when a reaction cannot produce exactly one usable rule."""


def reaction_smiles_from_node(node: dict[str, Any]) -> str:
    metadata = node.get("metadata") or {}
    smiles = (
        node.get("smiles")
        or metadata.get("smiles")
        or metadata.get("mapped_reaction_smiles")
        or metadata.get("rsmi")
    )
    if not smiles:
        raise ValueError("reaction node has no mapped reaction SMILES")
    return smiles


def normalize_route_tree(route: dict[str, Any]) -> dict[str, Any]:
    """Return a route copy whose reaction nodes contain mapped reaction SMILES."""

    route = copy.deepcopy(route)

    def visit(node: dict[str, Any]) -> None:
        if node.get("type") == "reaction":
            node["smiles"] = reaction_smiles_from_node(node)
        for child in node.get("children", []) or []:
            if isinstance(child, dict):
                visit(child)

    visit(route)
    return route


def route_items(routes_json: Any) -> Iterable[tuple[Any, dict[str, Any]]]:
    if isinstance(routes_json, list):
        for route_id, route in enumerate(routes_json):
            yield route_id, route
        return

    if isinstance(routes_json, dict):
        for route_id, route in routes_json.items():
            yield route_id, route
        return

    raise TypeError(f"unsupported routes JSON root: {type(routes_json)!r}")


def child_reaction_nodes(reaction_node: dict[str, Any]) -> list[dict[str, Any]]:
    children = []
    for mol_node in reaction_node.get("children", []) or []:
        if mol_node.get("type") != "mol":
            continue
        for child in mol_node.get("children", []) or []:
            if child.get("type") == "reaction":
                children.append(child)
    return children


def root_reaction_nodes(route: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        child
        for child in route.get("children", []) or []
        if isinstance(child, dict) and child.get("type") == "reaction"
    ]


def reaction_paths_from_node(
    reaction_node: dict[str, Any],
    step_by_reaction_smiles: dict[str, ReactionRuleStep],
) -> list[list[ReactionRuleStep]]:
    step = step_by_reaction_smiles[reaction_smiles_from_node(reaction_node)]
    children = child_reaction_nodes(reaction_node)
    if not children:
        return [[step]]

    paths: list[list[ReactionRuleStep]] = []
    for child in children:
        for suffix in reaction_paths_from_node(child, step_by_reaction_smiles):
            paths.append([step] + suffix)
    return paths


def adjacent_centers_overlap(left: ReactionRuleStep, right: ReactionRuleStep) -> bool:
    return bool(left.center_atoms & right.center_atoms)


def valid_composite_sequences(
    path: list[ReactionRuleStep],
    *,
    min_length: int,
    max_length: int | None,
) -> Iterable[tuple[str, ...]]:
    for sequence, _target_smiles in valid_composite_sequence_occurrences(
        path,
        min_length=min_length,
        max_length=max_length,
    ):
        yield sequence


def valid_composite_sequence_occurrences(
    path: list[ReactionRuleStep],
    *,
    min_length: int,
    max_length: int | None,
) -> Iterable[tuple[tuple[str, ...], str]]:
    if len(path) < min_length:
        return

    segment: list[ReactionRuleStep] = [path[0]]

    def emit_segment(
        steps: list[ReactionRuleStep],
    ) -> Iterable[tuple[tuple[str, ...], str]]:
        if len(steps) < min_length:
            return
        upper = len(steps) if max_length is None else min(len(steps), max_length)
        for start in range(len(steps)):
            for end in range(start + min_length, min(len(steps), start + upper) + 1):
                yield (
                    tuple(step.rule_smarts for step in steps[start:end]),
                    steps[start].target_smiles,
                )

    for step in path[1:]:
        if adjacent_centers_overlap(segment[-1], step):
            segment.append(step)
            continue
        yield from emit_segment(segment)
        segment = [step]

    yield from emit_segment(segment)


class SynPlannerRuleExtractor:
    def __init__(self, config: Any):
        from synplan.chem.data.standardizing import RemoveReagentsStandardizer

        self.config = config
        self.standardizer = RemoveReagentsStandardizer()
        self.cache: dict[str, ReactionRuleStep | None] = {}

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "SynPlannerRuleExtractor":
        from synplan.utils.config import RuleExtractionConfig

        if args.config:
            config = RuleExtractionConfig.from_yaml(str(resolve_existing_path(args.config)))
        else:
            config = RuleExtractionConfig(
                min_popularity=1,
                single_product_only=True,
                environment_atom_count=args.environment_atom_count,
                multicenter_rules=True,
                include_rings=args.include_rings,
                include_func_groups=False,
                keep_leaving_groups=args.keep_leaving_groups,
                keep_incoming_groups=args.keep_incoming_groups,
                keep_reagents=False,
                reactor_validation=args.reactor_validation,
            )
        return cls(config)

    def extract(self, reaction_smiles: str) -> tuple[ReactionRuleStep | None, bool]:
        """Return `(step, cache_hit)` for one mapped reaction SMILES."""

        if reaction_smiles in self.cache:
            return self.cache[reaction_smiles], True

        from chython import smiles as parse_smiles
        from synplan.chem.reaction_rules.extraction import (
            _rule_to_reactor_smarts,
            extract_rules,
        )

        reaction = parse_smiles(reaction_smiles)
        standardized = self.standardizer(reaction)
        center_atoms = frozenset((~standardized).center_atoms)
        rules, skipped = extract_rules(self.config, reaction)
        if skipped or not rules:
            self.cache[reaction_smiles] = None
            return None, False
        if len(rules) != 1:
            raise RuleExtractionError(
                "composite extraction expects one multicenter rule per reaction; "
                f"got {len(rules)} rules"
            )

        rule_smarts = _rule_to_reactor_smarts(rules[0])
        step = ReactionRuleStep(
            rule_smarts=rule_smarts,
            center_atoms=center_atoms,
            reaction_smiles=reaction_smiles,
        )
        self.cache[reaction_smiles] = step
        return step, False


def collect_reaction_contexts(route: dict[str, Any]) -> list[tuple[str, str]]:
    contexts: list[tuple[str, str]] = []

    def visit(node: dict[str, Any]) -> None:
        if node.get("type") == "mol":
            target_smiles = node.get("smiles", "")
            for child in node.get("children", []) or []:
                if isinstance(child, dict) and child.get("type") == "reaction":
                    contexts.append((reaction_smiles_from_node(child), target_smiles))
        for child in node.get("children", []) or []:
            if isinstance(child, dict):
                visit(child)

    visit(route)
    return contexts


def collect_reaction_smiles(route: dict[str, Any]) -> list[str]:
    return [reaction_smiles for reaction_smiles, _ in collect_reaction_contexts(route)]


def extract_route_composites(
    route: dict[str, Any],
    rule_extractor: SynPlannerRuleExtractor,
    *,
    min_length: int,
    max_length: int | None,
    stats: RouteProcessingStats,
) -> dict[tuple[str, ...], set[str]]:
    route = normalize_route_tree(route)
    step_by_reaction_smiles: dict[str, ReactionRuleStep] = {}

    for reaction_smiles, target_smiles in collect_reaction_contexts(route):
        stats.reactions_seen += 1
        step, cache_hit = rule_extractor.extract(reaction_smiles)
        if cache_hit:
            stats.reaction_rule_cache_hits += 1
        else:
            stats.reaction_rule_cache_misses += 1
        if step is None:
            stats.skipped_reactions += 1
            continue
        step_by_reaction_smiles[reaction_smiles] = ReactionRuleStep(
            rule_smarts=step.rule_smarts,
            center_atoms=step.center_atoms,
            reaction_smiles=step.reaction_smiles,
            target_smiles=target_smiles,
        )

    sequences: dict[tuple[str, ...], set[str]] = defaultdict(set)
    for root in root_reaction_nodes(route):
        try:
            paths = reaction_paths_from_node(root, step_by_reaction_smiles)
        except KeyError:
            continue
        for path in paths:
            for sequence, target_smiles in valid_composite_sequence_occurrences(
                path,
                min_length=min_length,
                max_length=max_length,
            ):
                sequences[sequence].add(target_smiles)
    return sequences


def run(args: argparse.Namespace) -> int:
    setup_runtime_cache_dirs()
    if args.min_length < 2:
        raise ValueError("--min-length must be at least 2")
    if args.max_length is not None and args.max_length <= 0:
        args.max_length = None
    if args.max_length is not None and args.max_length < args.min_length:
        raise ValueError("--max-length must be greater than or equal to --min-length")

    rule_extractor = SynPlannerRuleExtractor.from_args(args)

    routes_json = read_json(args.routes_json)

    references_by_sequence: dict[tuple[str, ...], set[Any]] = defaultdict(set)
    target_molecules_by_sequence: dict[tuple[str, ...], dict[Any, set[str]]] = (
        defaultdict(lambda: defaultdict(set))
    )
    errors: list[dict[str, Any]] = []
    stats = RouteProcessingStats()

    for index, (route_id, route) in enumerate(route_items(routes_json), start=1):
        if args.limit is not None and index > args.limit:
            break
        stats.routes_seen += 1
        try:
            route_sequences = extract_route_composites(
                route,
                rule_extractor,
                min_length=args.min_length,
                max_length=args.max_length,
                stats=stats,
            )
            if route_sequences:
                stats.routes_with_composites += 1
            for sequence, target_molecules in route_sequences.items():
                references_by_sequence[sequence].add(route_id)
                target_molecules_by_sequence[sequence][route_id].update(
                    target_molecules
                )
        except Exception as exc:
            stats.errors += 1
            if not args.ignore_errors:
                raise
            errors.append(
                {
                    "route_id": route_id,
                    "stage": "extract_route_composites",
                    "error_type": type(exc).__qualname__,
                    "message": str(exc) or traceback.format_exc(limit=1).strip(),
                }
            )

        if args.progress_interval and index % args.progress_interval == 0:
            print(
                f"processed routes={index} composite_rules={len(references_by_sequence)} "
                f"errors={stats.errors}",
                flush=True,
            )

    output_summary = write_composite_rules(
        args.output,
        references_by_sequence,
        target_molecules_by_sequence=target_molecules_by_sequence,
    )
    write_errors(args.output, errors)

    summary = {
        "routes_json": str(args.routes_json),
        "routes_seen": stats.routes_seen,
        "routes_with_composite_rules": stats.routes_with_composites,
        "reactions_seen": stats.reactions_seen,
        "reaction_rule_cache_hits": stats.reaction_rule_cache_hits,
        "reaction_rule_cache_misses": stats.reaction_rule_cache_misses,
        "skipped_reactions": stats.skipped_reactions,
        "errors": stats.errors,
        "unique_composite_rules": len(references_by_sequence),
        "target_molecule_occurrences": sum(
            len(targets)
            for route_targets in target_molecules_by_sequence.values()
            for targets in route_targets.values()
        ),
        "min_length": args.min_length,
        "max_length": args.max_length,
        "output_prefix": str(args.output.with_suffix("")),
        **output_summary,
    }
    summary_path = write_summary(args.output, summary)
    summary["summary_file"] = str(summary_path)
    write_summary(args.output, summary)

    print(json.dumps(summary, indent=2), flush=True)
    return 0
