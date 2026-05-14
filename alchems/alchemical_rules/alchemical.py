from __future__ import annotations

import argparse
import json
import sys
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

PROJECT_ROOT = Path(__file__).resolve().parents[2]

if __package__ in (None, "") and str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from alchems.composite_rules.unwrap import (
    RuleApplicationError,
    split_composite_rule,
    unwrap_rule_sequence,
)
from alchems.io import (
    expand_composite_rule_tsv_paths,
    iter_composite_rule_applications,
    resolve_alchemical_output_paths,
    resolve_existing_path,
    setup_runtime_cache_dirs,
    write_alchemical_errors,
    write_alchemical_rules_tsv,
    write_json,
    write_pseudo_reactions_smi,
)


@dataclass(frozen=True)
class ExtractedAlchemicalRule:
    rule_smarts: str
    cgr_key: str


@dataclass
class PseudoReactionRecord:
    pseudo_reaction_id: str
    alchemical_cgr: str
    reaction_smiles: str
    source_tsv: str
    source_row: int
    route_ids: tuple[str, ...]
    target_smiles: str
    composite_size: int
    composite_rule: str


@dataclass
class AlchemicalRuleAggregate:
    rule_smarts: str
    cgr_key: str
    route_ids: set[str] = field(default_factory=set)
    target_molecules: set[str] = field(default_factory=set)
    composite_rules: set[str] = field(default_factory=set)
    composite_sizes: set[int] = field(default_factory=set)
    source_rows: set[str] = field(default_factory=set)
    pseudo_reaction_ids: list[str] = field(default_factory=list)


@dataclass
class AlchemicalCollectionStats:
    composite_rows_seen: int = 0
    applications_seen: int = 0
    pseudo_reactions_built: int = 0
    alchemical_rules_extracted: int = 0
    skipped_unwrap_applications: int = 0
    skipped_rule_extractions: int = 0
    skipped_rule_extraction_errors: int = 0
    errors: int = 0


def is_standardization_error(exc: Exception) -> bool:
    return type(exc).__name__ == "StandardizationError"


def compose_pseudo_reaction_smiles(
    target_smiles: str,
    composite_rule: str,
) -> str:
    from chython.containers.reaction import ReactionContainer

    unwrapped = unwrap_rule_sequence(
        target_smiles,
        split_composite_rule(composite_rule),
        route_id=0,
        rule_key_prefix="composite",
        mark_leaves_in_stock=True,
    )
    reactants, product = normalize_pseudo_reaction_mapping(
        unwrapped.leaf_molecules,
        unwrapped.target_molecule,
    )
    pseudo_reaction = ReactionContainer(
        reactants=tuple(reactants),
        products=(product,),
    )
    return format(pseudo_reaction, "m")


def normalize_pseudo_reaction_mapping(
    reactants: Iterable[Any],
    product: Any,
) -> tuple[list[Any], Any]:
    product = product.copy()
    product_atoms = {
        atom_number: atom.atomic_number for atom_number, atom in product.atoms()
    }
    used_reactant_numbers: set[int] = set()
    next_atom_number = max(product_atoms, default=0) + 1
    normalized_reactants = []

    for reactant in reactants:
        reactant = reactant.copy()
        remapping: dict[int, int] = {}
        for atom_number, atom in reactant.atoms():
            product_atomic_number = product_atoms.get(atom_number)
            can_keep_number = (
                atom_number not in used_reactant_numbers
                and (
                    product_atomic_number is None
                    or product_atomic_number == atom.atomic_number
                )
            )
            if can_keep_number:
                used_reactant_numbers.add(atom_number)
                continue

            while (
                next_atom_number in product_atoms
                or next_atom_number in used_reactant_numbers
            ):
                next_atom_number += 1
            remapping[atom_number] = next_atom_number
            used_reactant_numbers.add(next_atom_number)
            next_atom_number += 1

        if remapping:
            reactant.remap(remapping)
        normalized_reactants.append(reactant)

    return normalized_reactants, product


def rule_cgr_key(rule_smarts: str) -> str:
    from chython import smarts
    from chython.containers.reaction import ReactionContainer
    from chython.reactor import Reactor

    try:
        return str(~smarts(rule_smarts))
    except Exception:
        reactor = Reactor.from_smarts(rule_smarts, delete_atoms=False)
        reaction = ReactionContainer(
            reactor.__dict__["_patterns"],
            reactor.__dict__["_products"],
        )
        return str(~reaction)


class AlchemicalRuleExtractor:
    def __init__(self, config: Any):
        self.config = config
        self.cache: dict[str, ExtractedAlchemicalRule | None] = {}

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "AlchemicalRuleExtractor":
        from synplan.utils.config import RuleExtractionConfig

        if args.config:
            config = RuleExtractionConfig.from_yaml(
                str(resolve_existing_path(args.config))
            )
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

    def extract(self, reaction_smiles: str) -> ExtractedAlchemicalRule | None:
        if reaction_smiles in self.cache:
            return self.cache[reaction_smiles]

        from chython import smiles as parse_smiles
        from synplan.chem.reaction_rules.extraction import (
            _rule_to_reactor_smarts,
            extract_rules,
        )

        reaction = parse_smiles(reaction_smiles)
        rules, skipped = extract_rules(self.config, reaction)
        if skipped or not rules:
            self.cache[reaction_smiles] = None
            return None
        if len(rules) != 1:
            raise ValueError(f"expected one alchemical rule, got {len(rules)}")

        rule = rules[0]
        extracted = ExtractedAlchemicalRule(
            rule_smarts=_rule_to_reactor_smarts(rule),
            cgr_key=str(~rule),
        )
        self.cache[reaction_smiles] = extracted
        return extracted


def collect_alchemical_rules(
    composite_rule_tsvs: list[Path],
    extractor: AlchemicalRuleExtractor,
    *,
    limit_rows: int | None = None,
    limit_applications: int | None = None,
    ignore_errors: bool = False,
    progress_interval: int = 250,
) -> tuple[
    dict[str, AlchemicalRuleAggregate],
    list[PseudoReactionRecord],
    AlchemicalCollectionStats,
    list[dict[str, Any]],
]:
    aggregates: dict[str, AlchemicalRuleAggregate] = {}
    pseudo_reactions: list[PseudoReactionRecord] = []
    errors: list[dict[str, Any]] = []
    stats = AlchemicalCollectionStats()
    rows_seen: set[tuple[Path, int]] = set()

    for application in iter_composite_rule_applications(composite_rule_tsvs):
        if (
            limit_applications is not None
            and stats.applications_seen >= limit_applications
        ):
            break

        row_key = (application.source_tsv, application.row_index)
        if row_key not in rows_seen:
            if limit_rows is not None and stats.composite_rows_seen >= limit_rows:
                break
            rows_seen.add(row_key)
            stats.composite_rows_seen += 1

        stats.applications_seen += 1

        try:
            pseudo_reaction_smiles = compose_pseudo_reaction_smiles(
                application.target_smiles,
                application.composite_rule,
            )
            stats.pseudo_reactions_built += 1
            extracted = extractor.extract(pseudo_reaction_smiles)
            if extracted is None:
                stats.skipped_rule_extractions += 1
                continue
            stats.alchemical_rules_extracted += 1

            pseudo_reaction_id = f"p{len(pseudo_reactions)}"
            pseudo_reactions.append(
                PseudoReactionRecord(
                    pseudo_reaction_id=pseudo_reaction_id,
                    alchemical_cgr=extracted.cgr_key,
                    reaction_smiles=pseudo_reaction_smiles,
                    source_tsv=str(application.source_tsv),
                    source_row=application.row_index,
                    route_ids=application.route_ids,
                    target_smiles=application.target_smiles,
                    composite_size=application.composite_size,
                    composite_rule=application.composite_rule,
                )
            )

            aggregate = aggregates.get(extracted.cgr_key)
            if aggregate is None:
                aggregate = AlchemicalRuleAggregate(
                    rule_smarts=extracted.rule_smarts,
                    cgr_key=extracted.cgr_key,
                )
                aggregates[extracted.cgr_key] = aggregate

            aggregate.route_ids.update(application.route_ids)
            aggregate.target_molecules.add(application.target_smiles)
            aggregate.composite_rules.add(application.composite_rule)
            aggregate.composite_sizes.add(application.composite_size)
            aggregate.source_rows.add(
                f"{application.source_tsv.name}:{application.row_index}"
            )
            aggregate.pseudo_reaction_ids.append(pseudo_reaction_id)
        except Exception as exc:
            if isinstance(exc, RuleApplicationError):
                stats.skipped_unwrap_applications += 1
                continue
            if is_standardization_error(exc):
                stats.skipped_rule_extraction_errors += 1
                continue

            stats.errors += 1
            errors.append(
                {
                    "source_tsv": str(application.source_tsv),
                    "row_index": application.row_index,
                    "target_smiles": application.target_smiles,
                    "error_type": type(exc).__qualname__,
                    "message": str(exc) or traceback.format_exc(limit=1).strip(),
                }
            )
            if not ignore_errors:
                raise

        if progress_interval and stats.applications_seen % progress_interval == 0:
            print(
                "processed applications="
                f"{stats.applications_seen} alchemical_rules={len(aggregates)} "
                f"skipped_unwrap={stats.skipped_unwrap_applications} "
                f"errors={stats.errors}",
                flush=True,
            )

    return aggregates, pseudo_reactions, stats, errors


def run(args: argparse.Namespace) -> int:
    setup_runtime_cache_dirs()
    composite_rule_tsvs = expand_composite_rule_tsv_paths(args.composite_rule_tsv)
    extractor = AlchemicalRuleExtractor.from_args(args)
    aggregates, pseudo_reactions, stats, errors = collect_alchemical_rules(
        composite_rule_tsvs,
        extractor,
        limit_rows=args.limit_rows,
        limit_applications=args.limit_applications,
        ignore_errors=args.ignore_errors,
        progress_interval=args.progress_interval,
    )

    rules_path, smi_path, summary_path, error_path = resolve_alchemical_output_paths(
        args.output,
        composite_rule_tsvs,
        output_smi=args.output_smi,
        summary=args.summary,
        errors=args.errors,
    )

    output_stats = write_alchemical_rules_tsv(rules_path, aggregates)
    write_pseudo_reactions_smi(smi_path, pseudo_reactions, aggregates)
    write_alchemical_errors(error_path, errors)

    summary = {
        "composite_rule_tsv": [str(path) for path in composite_rule_tsvs],
        "output": str(rules_path),
        "pseudo_reactions_smi": str(smi_path),
        "errors_file": str(error_path) if errors else None,
        "composite_rows_seen": stats.composite_rows_seen,
        "applications_seen": stats.applications_seen,
        "pseudo_reactions_built": stats.pseudo_reactions_built,
        "alchemical_rules_extracted": stats.alchemical_rules_extracted,
        "skipped_unwrap_applications": stats.skipped_unwrap_applications,
        "skipped_rule_extractions": stats.skipped_rule_extractions,
        "skipped_rule_extraction_errors": stats.skipped_rule_extraction_errors,
        "errors": stats.errors,
        "unique_alchemical_rules": len(aggregates),
        **output_stats,
    }
    write_json(summary_path, summary)
    summary["summary_file"] = str(summary_path)
    write_json(summary_path, summary)
    print(json.dumps(summary, indent=2), flush=True)
    return 0

