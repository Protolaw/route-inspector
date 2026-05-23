from __future__ import annotations

import argparse
import json
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

PROJECT_ROOT = Path(__file__).resolve().parents[2]

if __package__ in (None, "") and str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from route_inspector.composite_rules.unwrap import (
    RuleApplicationError,
    split_composite_rule,
    unwrap_rule_sequence,
)
from route_inspector.io import (
    expand_composite_rule_tsv_paths,
    iter_composite_rule_applications,
    normalize_n_cpu,
    resolve_alchemical_output_paths,
    resolve_existing_path,
    setup_runtime_cache_dirs,
    write_alchemical_errors,
    write_alchemical_rules_tsv,
    write_json,
    write_pseudo_reactions_smi,
    write_standard_sidecars,
)


_ALCHEMICAL_WORKER_EXTRACTOR: AlchemicalRuleExtractor | None = None


@dataclass(frozen=True)
class ExtractedAlchemicalRule:
    rule_smarts: str
    cgr_key: str
    query_cgr: Any


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
    query_cgr: Any | None = None
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
    """Return whether standardization error matches the expected condition.

    This step supports alchemical rule extraction from composite-rule applications and
    keeps duplicate rules merged by structural identity.
    """
    return type(exc).__name__ == "StandardizationError"


def compose_pseudo_reaction_smiles(
    target_smiles: str,
    composite_rule: str,
) -> str:
    """Compose pseudo reaction SMILES from normalized inputs.

    The pseudo-reaction collapses an unwrapped route into target-to-stock reactants
    before SynPlanner rule extraction is run again.
    """
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
    """Normalize pseudo reaction mapping for route-inspector processing.

    The pseudo-reaction collapses an unwrapped route into target-to-stock reactants
    before SynPlanner rule extraction is run again.
    """
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


def rule_query_cgr(rule_smarts: str) -> Any:
    """Return rule query CGR used for rule extraction or comparison.

    The comparison uses composed QueryCGR structure rather than raw SMARTS strings,
    which avoids treating equivalent atom-map numbering as different.
    """
    from chython import smarts
    from chython.containers.reaction import ReactionContainer
    from chython.reactor import Reactor

    try:
        return smarts(rule_smarts).compose()
    except Exception:
        reactor = Reactor.from_smarts(rule_smarts, delete_atoms=False)
        reaction = ReactionContainer(
            reactor.__dict__["_patterns"],
            reactor.__dict__["_products"],
        )
        return reaction.compose()


def rule_cgr_key(rule_smarts: str) -> str:
    """Return rule CGR key used for rule extraction or comparison.

    The comparison uses composed QueryCGR structure rather than raw SMARTS strings,
    which avoids treating equivalent atom-map numbering as different.
    """
    return json.dumps(
        query_cgr_coarse_signature(rule_query_cgr(rule_smarts)),
        sort_keys=True,
        default=repr,
    )


def freeze_query_value(value: Any) -> Any:
    """Convert query value into a hashable representation.

    The comparison uses composed QueryCGR structure rather than raw SMARTS strings,
    which avoids treating equivalent atom-map numbering as different.
    """
    if isinstance(value, (tuple, list)):
        return tuple(freeze_query_value(item) for item in value)
    if isinstance(value, set):
        return tuple(sorted(freeze_query_value(item) for item in value))
    if isinstance(value, (int, float, str, bool, type(None))):
        return value
    return repr(value)


def query_atom_signature(atom: Any) -> tuple[Any, ...]:
    """Return query atom signature used for structural rule comparison.

    The comparison uses composed QueryCGR structure rather than raw SMARTS strings,
    which avoids treating equivalent atom-map numbering as different.
    """
    return (
        atom.atomic_number,
        repr(freeze_query_value(atom.charge)),
        repr(freeze_query_value(atom.p_charge)),
        repr(freeze_query_value(atom.is_radical)),
        repr(freeze_query_value(atom.p_is_radical)),
        repr(freeze_query_value(atom.neighbors)),
        repr(freeze_query_value(atom.p_neighbors)),
        repr(freeze_query_value(atom.hybridization)),
        repr(freeze_query_value(atom.p_hybridization)),
        repr(freeze_query_value(atom.isotope)),
    )


def query_bond_signature(bond: Any) -> tuple[Any, ...]:
    """Return query bond signature used for structural rule comparison.

    The comparison uses composed QueryCGR structure rather than raw SMARTS strings,
    which avoids treating equivalent atom-map numbering as different.
    """
    return (
        repr(freeze_query_value(bond.order)),
        repr(freeze_query_value(bond.p_order)),
    )


def query_cgr_adjacency(query_cgr: Any) -> dict[int, dict[int, tuple[Any, ...]]]:
    """Return query CGR adjacency used for structural rule comparison.

    The comparison uses composed QueryCGR structure rather than raw SMARTS strings,
    which avoids treating equivalent atom-map numbering as different.
    """
    adjacency = {atom_number: {} for atom_number, _atom in query_cgr.atoms()}
    for atom_1, atom_2, bond in query_cgr.bonds():
        signature = query_bond_signature(bond)
        adjacency[atom_1][atom_2] = signature
        adjacency[atom_2][atom_1] = signature
    return adjacency


def query_cgr_coarse_signature(query_cgr: Any) -> tuple[Any, ...]:
    """Return query CGR coarse signature used for structural rule comparison.

    The comparison uses composed QueryCGR structure rather than raw SMARTS strings,
    which avoids treating equivalent atom-map numbering as different.
    """
    atom_signatures = {
        atom_number: query_atom_signature(atom)
        for atom_number, atom in query_cgr.atoms()
    }
    bond_signatures = []
    for atom_1, atom_2, bond in query_cgr.bonds():
        endpoints = sorted((atom_signatures[atom_1], atom_signatures[atom_2]))
        bond_signatures.append((tuple(endpoints), query_bond_signature(bond)))
    return (
        tuple(sorted(atom_signatures.values())),
        tuple(sorted(bond_signatures)),
    )


def query_cgr_isomorphic(left: Any, right: Any) -> bool:
    """Return query CGR isomorphic used for structural rule comparison.

    The comparison uses composed QueryCGR structure rather than raw SMARTS strings,
    which avoids treating equivalent atom-map numbering as different.
    """
    if left.atoms_count != right.atoms_count:
        return False
    if left.bonds_count != right.bonds_count:
        return False

    left_atoms = {
        atom_number: query_atom_signature(atom) for atom_number, atom in left.atoms()
    }
    right_atoms = {
        atom_number: query_atom_signature(atom) for atom_number, atom in right.atoms()
    }
    left_adjacency = query_cgr_adjacency(left)
    right_adjacency = query_cgr_adjacency(right)

    candidates = {
        left_atom: [
            right_atom
            for right_atom, right_signature in right_atoms.items()
            if right_signature == left_signature
            and len(right_adjacency[right_atom]) == len(left_adjacency[left_atom])
        ]
        for left_atom, left_signature in left_atoms.items()
    }
    if any(not atom_candidates for atom_candidates in candidates.values()):
        return False

    order = sorted(
        left_atoms,
        key=lambda atom_number: (
            len(candidates[atom_number]),
            -len(left_adjacency[atom_number]),
        ),
    )
    mapping: dict[int, int] = {}
    used_right_atoms: set[int] = set()

    def mapping_is_consistent(left_atom: int, right_atom: int) -> bool:
        """Return whether a candidate atom mapping is consistent with prior matches.

        This step supports alchemical rule extraction from composite-rule applications
        and keeps duplicate rules merged by structural identity.
        """
        for mapped_left_atom, mapped_right_atom in mapping.items():
            left_bond = left_adjacency[left_atom].get(mapped_left_atom)
            right_bond = right_adjacency[right_atom].get(mapped_right_atom)
            if (left_bond is None) != (right_bond is None):
                return False
            if left_bond is not None and left_bond != right_bond:
                return False
        return True

    def search(index: int) -> bool:
        """Search candidate atom mappings recursively.

        This step supports alchemical rule extraction from composite-rule applications
        and keeps duplicate rules merged by structural identity.
        """
        if index == len(order):
            return True
        left_atom = order[index]
        for right_atom in candidates[left_atom]:
            if right_atom in used_right_atoms:
                continue
            if not mapping_is_consistent(left_atom, right_atom):
                continue
            mapping[left_atom] = right_atom
            used_right_atoms.add(right_atom)
            if search(index + 1):
                return True
            used_right_atoms.remove(right_atom)
            del mapping[left_atom]
        return False

    return search(0)


def matching_aggregate(
    aggregate_buckets: dict[tuple[Any, ...], list[AlchemicalRuleAggregate]],
    query_cgr: Any,
) -> AlchemicalRuleAggregate | None:
    """Find an existing alchemical aggregate with an isomorphic QueryCGR.

    The aggregation links each alchemical rule back to the composite rules, targets, and
    route IDs that produced it.
    """
    for aggregate in aggregate_buckets.get(query_cgr_coarse_signature(query_cgr), []):
        if aggregate.query_cgr is not None and query_cgr_isomorphic(
            aggregate.query_cgr,
            query_cgr,
        ):
            return aggregate
    return None


def collection_error_row(
    application: Any,
    *,
    stage: str,
    exc: Exception,
    alchemical_rule: str = "",
) -> dict[str, Any]:
    """Build a TSV row describing a failed alchemical-rule application.

    This step supports alchemical rule extraction from composite-rule applications and
    keeps duplicate rules merged by structural identity.
    """
    return {
        "source_tsv": str(application.source_tsv),
        "row_index": application.row_index,
        "Target_smiles": application.target_smiles,
        "Composite_rule": application.composite_rule,
        "Composite_size": application.composite_size,
        "Route_ids": ",".join(application.route_ids),
    }


class AlchemicalRuleExtractor:
    def __init__(self, config: Any):
        """Initialize this object with its resolved configuration.

        This step supports alchemical rule extraction from composite-rule applications
        and keeps duplicate rules merged by structural identity.
        """
        self.config = config
        self.cache: dict[str, ExtractedAlchemicalRule | None] = {}

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "AlchemicalRuleExtractor":
        """Build an extractor instance from parsed CLI arguments.

        This step supports alchemical rule extraction from composite-rule applications
        and keeps duplicate rules merged by structural identity.
        """
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
        """Extract one rule from a mapped reaction SMILES string.

        This step supports alchemical rule extraction from composite-rule applications
        and keeps duplicate rules merged by structural identity.
        """
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
        query_cgr = rule.compose()
        rule_smarts = _rule_to_reactor_smarts(rule)
        extracted = ExtractedAlchemicalRule(
            rule_smarts=rule_smarts,
            cgr_key=rule_cgr_key(rule_smarts),
            query_cgr=query_cgr,
        )
        self.cache[reaction_smiles] = extracted
        return extracted


def alchemical_extractor_args_dict(args: argparse.Namespace) -> dict[str, Any]:
    """Serialize alchemical extractor settings for worker initialization.

    This step supports alchemical rule extraction from composite-rule applications and
    keeps duplicate rules merged by structural identity.
    """
    return {
        "config": str(args.config) if getattr(args, "config", None) else None,
        "environment_atom_count": getattr(args, "environment_atom_count", 1),
        "include_rings": getattr(args, "include_rings", False),
        "keep_leaving_groups": getattr(args, "keep_leaving_groups", True),
        "keep_incoming_groups": getattr(args, "keep_incoming_groups", False),
        "reactor_validation": getattr(args, "reactor_validation", False),
    }


def _init_alchemical_worker(extractor_args: dict[str, Any]) -> None:
    """Run the worker entry point for init alchemical.

    This step supports alchemical rule extraction from composite-rule applications and
    keeps duplicate rules merged by structural identity.
    """
    global _ALCHEMICAL_WORKER_EXTRACTOR
    setup_runtime_cache_dirs()
    _ALCHEMICAL_WORKER_EXTRACTOR = AlchemicalRuleExtractor.from_args(
        argparse.Namespace(**extractor_args)
    )


def _alchemical_application_worker(item: tuple[int, Any]) -> dict[str, Any]:
    """Run the worker entry point for alchemical application.

    This step supports alchemical rule extraction from composite-rule applications and
    keeps duplicate rules merged by structural identity.
    """
    _application_index, application = item
    try:
        if _ALCHEMICAL_WORKER_EXTRACTOR is None:
            raise RuntimeError("alchemical worker was not initialized")
        pseudo_reaction_smiles = compose_pseudo_reaction_smiles(
            application.target_smiles,
            application.composite_rule,
        )
        extracted = _ALCHEMICAL_WORKER_EXTRACTOR.extract(pseudo_reaction_smiles)
        if extracted is None:
            return {
                "status": "skipped_rule_extraction",
                "application": application,
                "pseudo_reaction_smiles": pseudo_reaction_smiles,
            }
        return {
            "status": "ok",
            "application": application,
            "pseudo_reaction_smiles": pseudo_reaction_smiles,
            "rule_smarts": extracted.rule_smarts,
        }
    except Exception as exc:
        if isinstance(exc, RuleApplicationError):
            return {
                "status": "skipped_unwrap",
                "application": application,
                "error": collection_error_row(
                    application,
                    stage="skipped_unwrap",
                    exc=exc,
                ),
            }
        if is_standardization_error(exc):
            return {
                "status": "skipped_rule_extraction_error",
                "application": application,
            }
        return {
            "status": "error",
            "application": application,
            "error": {
                **collection_error_row(application, stage="error", exc=exc),
                "error_type": type(exc).__qualname__,
                "message": str(exc) or traceback.format_exc(limit=1).strip(),
            },
        }


def select_composite_rule_applications(
    composite_rule_tsvs: list[Path],
    *,
    limit_rows: int | None = None,
    limit_applications: int | None = None,
) -> tuple[list[Any], int]:
    """Select composite rule applications for the next processing step.

    The aggregation links each alchemical rule back to the composite rules, targets, and
    route IDs that produced it.
    """
    applications = []
    rows_seen: set[tuple[Path, int]] = set()
    composite_rows_seen = 0
    for application in iter_composite_rule_applications(composite_rule_tsvs):
        if limit_applications is not None and len(applications) >= limit_applications:
            break
        row_key = (application.source_tsv, application.row_index)
        if row_key not in rows_seen:
            if limit_rows is not None and composite_rows_seen >= limit_rows:
                break
            rows_seen.add(row_key)
            composite_rows_seen += 1
        applications.append(application)
    return applications, composite_rows_seen


def merge_alchemical_success(
    *,
    application: Any,
    pseudo_reaction_smiles: str,
    rule_smarts: str,
    aggregates: dict[str, AlchemicalRuleAggregate],
    aggregate_buckets: dict[tuple[Any, ...], list[AlchemicalRuleAggregate]],
    pseudo_reactions: list[PseudoReactionRecord],
) -> None:
    """Merge alchemical success into aggregate results.

    The aggregation links each alchemical rule back to the composite rules, targets, and
    route IDs that produced it.
    """
    query_cgr = rule_query_cgr(rule_smarts)
    cgr_key = rule_cgr_key(rule_smarts)
    aggregate = matching_aggregate(aggregate_buckets, query_cgr)
    if aggregate is None:
        aggregate = AlchemicalRuleAggregate(
            rule_smarts=rule_smarts,
            cgr_key=cgr_key,
            query_cgr=query_cgr,
        )
        aggregates[aggregate.cgr_key] = aggregate
        aggregate_buckets.setdefault(
            query_cgr_coarse_signature(query_cgr),
            [],
        ).append(aggregate)

    pseudo_reaction_id = f"p{len(pseudo_reactions)}"
    pseudo_reactions.append(
        PseudoReactionRecord(
            pseudo_reaction_id=pseudo_reaction_id,
            alchemical_cgr=aggregate.cgr_key,
            reaction_smiles=pseudo_reaction_smiles,
            source_tsv=str(application.source_tsv),
            source_row=application.row_index,
            route_ids=application.route_ids,
            target_smiles=application.target_smiles,
            composite_size=application.composite_size,
            composite_rule=application.composite_rule,
        )
    )

    aggregate.route_ids.update(application.route_ids)
    aggregate.target_molecules.add(application.target_smiles)
    aggregate.composite_rules.add(application.composite_rule)
    aggregate.composite_sizes.add(application.composite_size)
    aggregate.source_rows.add(f"{application.source_tsv.name}:{application.row_index}")
    aggregate.pseudo_reaction_ids.append(pseudo_reaction_id)


def collect_alchemical_rules(
    composite_rule_tsvs: list[Path],
    extractor: AlchemicalRuleExtractor,
    *,
    limit_rows: int | None = None,
    limit_applications: int | None = None,
    ignore_errors: bool = False,
    progress_interval: int = 250,
    n_cpu: int = 1,
    extractor_args: argparse.Namespace | None = None,
) -> tuple[
    dict[str, AlchemicalRuleAggregate],
    list[PseudoReactionRecord],
    AlchemicalCollectionStats,
    list[dict[str, Any]],
]:
    """Collect alchemical rules from route-analysis input data.

    The aggregation links each alchemical rule back to the composite rules, targets, and
    route IDs that produced it.
    """
    aggregates: dict[str, AlchemicalRuleAggregate] = {}
    aggregate_buckets: dict[tuple[Any, ...], list[AlchemicalRuleAggregate]] = {}
    pseudo_reactions: list[PseudoReactionRecord] = []
    errors: list[dict[str, Any]] = []
    stats = AlchemicalCollectionStats()
    applications, stats.composite_rows_seen = select_composite_rule_applications(
        composite_rule_tsvs,
        limit_rows=limit_rows,
        limit_applications=limit_applications,
    )
    n_cpu = normalize_n_cpu(n_cpu)

    def consume_result(result: dict[str, Any]) -> None:
        """Merge one worker result into the aggregate state.

        This step supports alchemical rule extraction from composite-rule applications
        and keeps duplicate rules merged by structural identity.
        """
        application = result["application"]
        stats.applications_seen += 1
        status = result["status"]

        if status == "ok":
            stats.pseudo_reactions_built += 1
            stats.alchemical_rules_extracted += 1
            merge_alchemical_success(
                application=application,
                pseudo_reaction_smiles=result["pseudo_reaction_smiles"],
                rule_smarts=result["rule_smarts"],
                aggregates=aggregates,
                aggregate_buckets=aggregate_buckets,
                pseudo_reactions=pseudo_reactions,
            )
        elif status == "skipped_rule_extraction":
            stats.pseudo_reactions_built += 1
            stats.skipped_rule_extractions += 1
        elif status == "skipped_rule_extraction_error":
            stats.skipped_rule_extraction_errors += 1
        elif status == "skipped_unwrap":
            stats.skipped_unwrap_applications += 1
            errors.append(result["error"])
        else:
            stats.errors += 1
            errors.append(result["error"])
            if not ignore_errors:
                raise RuntimeError(
                    "alchemical rule extraction failed for "
                    f"{application.source_tsv}:{application.row_index}: "
                    f"{result['error'].get('error_type', 'Error')}: "
                    f"{result['error'].get('message', '')}"
                )

        if progress_interval and stats.applications_seen % progress_interval == 0:
            print(
                "processed applications="
                f"{stats.applications_seen} alchemical_rules={len(aggregates)} "
                f"skipped_unwrap={stats.skipped_unwrap_applications} "
                f"errors={stats.errors}",
                flush=True,
            )

    if n_cpu > 1 and applications:
        if extractor_args is None:
            raise ValueError("extractor_args is required when n_cpu > 1")
        with ProcessPoolExecutor(
            max_workers=n_cpu,
            initializer=_init_alchemical_worker,
            initargs=(alchemical_extractor_args_dict(extractor_args),),
        ) as executor:
            for result in executor.map(
                _alchemical_application_worker,
                enumerate(applications),
            ):
                consume_result(result)
    else:
        for application_index, application in enumerate(applications):
            try:
                pseudo_reaction_smiles = compose_pseudo_reaction_smiles(
                    application.target_smiles,
                    application.composite_rule,
                )
                extracted = extractor.extract(pseudo_reaction_smiles)
                if extracted is None:
                    result = {
                        "status": "skipped_rule_extraction",
                        "application": application,
                        "pseudo_reaction_smiles": pseudo_reaction_smiles,
                    }
                else:
                    result = {
                        "status": "ok",
                        "application": application,
                        "pseudo_reaction_smiles": pseudo_reaction_smiles,
                        "rule_smarts": extracted.rule_smarts,
                    }
            except Exception as exc:
                if isinstance(exc, RuleApplicationError):
                    result = {
                        "status": "skipped_unwrap",
                        "application": application,
                        "error": collection_error_row(
                            application,
                            stage="skipped_unwrap",
                            exc=exc,
                        ),
                    }
                elif is_standardization_error(exc):
                    result = {
                        "status": "skipped_rule_extraction_error",
                        "application": application,
                    }
                else:
                    result = {
                        "status": "error",
                        "application": application,
                        "error": {
                            **collection_error_row(
                                application,
                                stage="error",
                                exc=exc,
                            ),
                            "error_type": type(exc).__qualname__,
                            "message": str(exc)
                            or traceback.format_exc(limit=1).strip(),
                        },
                    }
            consume_result(result)

    if progress_interval and stats.applications_seen % progress_interval:
        print(
            "processed applications="
            f"{stats.applications_seen} alchemical_rules={len(aggregates)} "
            f"skipped_unwrap={stats.skipped_unwrap_applications} "
            f"errors={stats.errors}",
            flush=True,
        )

    return aggregates, pseudo_reactions, stats, errors


def run(args: argparse.Namespace) -> int:
    """Run this module command with parsed CLI arguments.

    This step supports alchemical rule extraction from composite-rule applications and
    keeps duplicate rules merged by structural identity.
    """
    setup_runtime_cache_dirs()
    if getattr(args, "output", None) is None:
        if getattr(args, "output_dir", None) is None:
            raise ValueError("either --output or --output-dir is required")
        args.output = args.output_dir
    elif getattr(args, "output_dir", None) is not None:
        raise ValueError("--output and --output-dir cannot be used together")
    composite_rule_tsvs = expand_composite_rule_tsv_paths(args.composite_rule_tsv)
    extractor = AlchemicalRuleExtractor.from_args(args)
    aggregates, pseudo_reactions, stats, errors = collect_alchemical_rules(
        composite_rule_tsvs,
        extractor,
        limit_rows=args.limit_rows,
        limit_applications=args.limit_applications,
        ignore_errors=args.ignore_errors,
        progress_interval=args.progress_interval,
        n_cpu=args.n_cpu,
        extractor_args=args,
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
        "n_cpu": normalize_n_cpu(args.n_cpu),
        **output_stats,
    }
    write_json(summary_path, summary)
    summary["summary_file"] = str(summary_path)
    write_json(summary_path, summary)
    sidecars = write_standard_sidecars(
        rules_path.parent,
        command_name="extract-alchemical-rules",
        summary=summary,
        errors=errors,
        input_files=composite_rule_tsvs,
        output_files={
            "alchemical_rules": rules_path,
            "pseudo_reactions_smi": smi_path,
            "summary": summary_path,
            "errors": error_path,
        },
        config_path=getattr(args, "config", None),
        cli_args=args,
    )
    summary["sidecar_files"] = sidecars
    write_json(summary_path, summary)
    print(json.dumps(summary, indent=2), flush=True)
    return 0
