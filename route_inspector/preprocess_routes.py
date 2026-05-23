from __future__ import annotations

import argparse
import copy
import json
import sys
import traceback
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from dataclasses import replace
from itertools import islice, permutations
from pathlib import Path
from typing import Any

from route_inspector.composite_rules.extract import (
    reaction_smiles_from_node,
    route_items,
)
from route_inspector.composite_rules.extract import (
    normalize_route_tree,
)
from route_inspector.composite_rules.unwrap import unwrap_rule_sequence
from route_inspector.io import (
    dataset_prefix_from_path,
    normalize_n_cpu,
    read_json,
    resolve_existing_path,
    setup_runtime_cache_dirs,
    stage_output_dir,
    write_json,
    write_standard_sidecars,
)
from route_inspector.protection.analysis import (
    ProtectionAnalysisConfig,
    ReactionRecord,
    RouteIndex,
    analyze_route_protection,
    build_route_index,
    detect_deprotections,
    molecule_node_smiles,
    parse_molecule,
    same_molecule,
)
from route_inspector.protection.chython_rules import (
    ProtectionRule,
    load_chython_protection_rules,
)


@dataclass(frozen=True)
class SingleCenterRule:
    rule_smarts: str
    center_atoms: frozenset[int]
    forward_change_kind: str = "unknown"
    forward_bonds_formed: tuple[tuple[int, int], ...] = ()
    forward_bonds_broken: tuple[tuple[int, int], ...] = ()
    forward_bonds_changed: tuple[tuple[int, int], ...] = ()


@dataclass(frozen=True)
class ReactionGranularity:
    reaction_smiles: str
    multicenter_rule_smarts: str
    single_center_rules: tuple[SingleCenterRule, ...]
    center_components: tuple[frozenset[int], ...]
    skipped: bool = False

    @property
    def is_multicenter(self) -> bool:
        return len(self.single_center_rules) > 1


@dataclass(frozen=True)
class SplitPlan:
    route_id: Any
    reaction_id: str
    parent_mol_id: str
    reaction_node: dict[str, Any]
    parent_mol_node: dict[str, Any]
    original_children: tuple[dict[str, Any], ...]
    original_reaction_smiles: str
    parent_smiles: str
    extraction: ReactionGranularity
    protection_matches: tuple[Any, ...]
    protection_atom_ids: frozenset[int]


@dataclass
class RoutePreprocessResult:
    route_id: Any
    route: dict[str, Any]
    normalized: bool
    modified: bool
    multicenter_reactions: int
    protection_multicenter_reactions: int
    split_reactions: int
    protection_split_reactions: int
    changes: list[dict[str, Any]]
    unresolved_reactions: list[dict[str, Any]]
    errors: list[dict[str, Any]]


class RouteSplitError(ValueError):
    """Raised when a multicenter reaction cannot be split."""


def reaction_side_bonds(
    molecules: Iterable[Any],
    center_atoms: frozenset[int],
) -> dict[tuple[int, int], str]:
    bonds: dict[tuple[int, int], str] = {}
    for molecule in molecules:
        for atom_1, atom_2, bond in molecule.bonds():
            atom_1 = int(atom_1)
            atom_2 = int(atom_2)
            if atom_1 not in center_atoms and atom_2 not in center_atoms:
                continue
            key = (min(atom_1, atom_2), max(atom_1, atom_2))
            try:
                value = str(int(bond))
            except Exception:
                value = str(bond)
            bonds[key] = value
    return bonds


def component_forward_change(
    reaction: Any,
    center_atoms: frozenset[int],
) -> dict[str, Any]:
    reactant_bonds = reaction_side_bonds(reaction.reactants, center_atoms)
    product_bonds = reaction_side_bonds(reaction.products, center_atoms)
    reactant_keys = set(reactant_bonds)
    product_keys = set(product_bonds)

    formed = tuple(sorted(product_keys - reactant_keys))
    broken = tuple(sorted(reactant_keys - product_keys))
    changed = tuple(
        sorted(
            key
            for key in reactant_keys & product_keys
            if reactant_bonds[key] != product_bonds[key]
        )
    )

    if formed and not broken:
        kind = "bond_forming"
    elif broken and not formed:
        kind = "bond_breaking"
    elif formed and broken:
        kind = "bond_forming_and_breaking"
    elif changed:
        kind = "bond_order_change"
    else:
        kind = "other"

    return {
        "forward_change_kind": kind,
        "forward_bonds_formed": formed,
        "forward_bonds_broken": broken,
        "forward_bonds_changed": changed,
    }


class SynPlannerGranularityExtractor:
    """Extract multicenter and single-center SynPlanner rules for one reaction."""

    def __init__(self, config: Any):
        from synplan.chem.data.standardizing import RemoveReagentsStandardizer

        self.config = config
        self.standardizer = RemoveReagentsStandardizer()
        self.cache: dict[str, ReactionGranularity] = {}

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "SynPlannerGranularityExtractor":
        from synplan.utils.config import RuleExtractionConfig

        if getattr(args, "config", None):
            config = RuleExtractionConfig.from_yaml(
                str(resolve_existing_path(args.config))
            )
        else:
            config = RuleExtractionConfig(
                min_popularity=1,
                single_product_only=True,
                environment_atom_count=getattr(args, "environment_atom_count", 1),
                multicenter_rules=True,
                include_rings=getattr(args, "include_rings", False),
                include_func_groups=False,
                keep_leaving_groups=getattr(args, "keep_leaving_groups", True),
                keep_incoming_groups=getattr(args, "keep_incoming_groups", False),
                keep_reagents=False,
                reactor_validation=getattr(args, "reactor_validation", False),
            )
        return cls(config)

    def _config_with(self, **updates: Any) -> Any:
        if hasattr(self.config, "model_copy"):
            return self.config.model_copy(update=updates)
        if hasattr(self.config, "copy"):
            return self.config.copy(update=updates)
        values = dict(getattr(self.config, "__dict__", {}))
        values.update(updates)
        return type(self.config)(**values)

    def extract(self, reaction_smiles: str) -> ReactionGranularity:
        if reaction_smiles in self.cache:
            return self.cache[reaction_smiles]

        from chython import smiles as parse_smiles
        from synplan.chem.reaction_rules.extraction import (
            _rule_to_reactor_smarts,
            create_rule,
        )

        reaction = parse_smiles(reaction_smiles)
        if getattr(self.config, "ignore_stereo", True):
            reaction = reaction.copy()
            reaction.clean_stereo()
        standardized = self.standardizer(reaction)

        if getattr(self.config, "single_product_only", True) and (
            len(standardized.products) != 1
        ):
            result = ReactionGranularity(
                reaction_smiles=reaction_smiles,
                multicenter_rule_smarts="",
                single_center_rules=(),
                center_components=(),
                skipped=True,
            )
            self.cache[reaction_smiles] = result
            return result

        cgr = ~standardized
        center_components = tuple(
            frozenset(int(atom_id) for atom_id in component)
            for component in islice(cgr.centers_list, 15)
        )
        skip_full_validation = len(center_components) > 1
        multicenter_config = self._config_with(multicenter_rules=True)
        single_center_config = self._config_with(multicenter_rules=False)

        multicenter_rule = create_rule(multicenter_config, standardized)
        multicenter_rule_smarts = _rule_to_reactor_smarts(multicenter_rule)

        seen_cgrs: dict[Any, SingleCenterRule] = {}
        for component in center_components:
            rule = create_rule(
                single_center_config,
                standardized,
                _restrict_center_atoms=set(component),
                _skip_full_reaction_validation=skip_full_validation,
            )
            rule_cgr = ~rule
            if rule_cgr in seen_cgrs:
                continue
            seen_cgrs[rule_cgr] = SingleCenterRule(
                rule_smarts=_rule_to_reactor_smarts(rule),
                center_atoms=component,
                **component_forward_change(standardized, component),
            )

        result = ReactionGranularity(
            reaction_smiles=reaction_smiles,
            multicenter_rule_smarts=multicenter_rule_smarts,
            single_center_rules=tuple(seen_cgrs.values()),
            center_components=center_components,
        )
        self.cache[reaction_smiles] = result
        return result


def same_molecule_smiles(left: str, right: str) -> bool:
    if left == right:
        return True
    if not left or not right:
        return False
    try:
        return same_molecule(parse_molecule(left), parse_molecule(right))
    except Exception:
        return False


def protection_atom_ids_from_matches(matches: Iterable[Any]) -> frozenset[int]:
    atom_ids: set[int] = set()
    for match in matches:
        atom_ids.update(
            int(atom_id) for atom_id in getattr(match, "protected_atom_ids", ())
        )
        for _query_atom, route_atom in getattr(match, "raw_mapping", ()):
            atom_ids.add(int(route_atom))
    return frozenset(atom_ids)


NON_PROTECTION_SPLIT_PRIORITY = {
    "bond_breaking": 0,
    "bond_forming_and_breaking": 1,
    "bond_order_change": 1,
    "other": 1,
    "unknown": 1,
    "bond_forming": 2,
}


def non_protection_rule_groups(
    rules: tuple[SingleCenterRule, ...],
) -> tuple[tuple[SingleCenterRule, ...], ...]:
    grouped: dict[int, list[SingleCenterRule]] = {}
    for rule in rules:
        priority = NON_PROTECTION_SPLIT_PRIORITY.get(rule.forward_change_kind, 1)
        grouped.setdefault(priority, []).append(rule)
    return tuple(tuple(grouped[priority]) for priority in sorted(grouped))


def grouped_rule_permutations(
    groups: tuple[tuple[SingleCenterRule, ...], ...],
) -> Iterable[tuple[SingleCenterRule, ...]]:
    if not groups:
        yield ()
        return
    first, *rest = groups
    for left in permutations(first):
        for right in grouped_rule_permutations(tuple(rest)):
            yield tuple(left + right)


def split_candidate_rule_orders(
    rules: tuple[SingleCenterRule, ...],
    protection_atom_ids: frozenset[int],
    *,
    deprotection_first: bool,
) -> Iterable[tuple[SingleCenterRule, ...]]:
    if protection_atom_ids:
        protection_rules = tuple(
            rule for rule in rules if rule.center_atoms & protection_atom_ids
        )
        other_rules = tuple(
            rule for rule in rules if not (rule.center_atoms & protection_atom_ids)
        )
        if protection_rules:
            groups = (
                (protection_rules, other_rules)
                if deprotection_first
                else (other_rules, protection_rules)
            )
        else:
            groups = non_protection_rule_groups(rules)
    else:
        groups = non_protection_rule_groups(rules)
    groups = tuple(group for group in groups if group)

    primary = tuple(rule for group in groups for rule in group)
    if not primary:
        return
    yielded = {tuple(rule.rule_smarts for rule in primary)}
    yield primary

    if len(rules) > 6:
        return

    for candidate in grouped_rule_permutations(groups):
        key = tuple(rule.rule_smarts for rule in candidate)
        if key in yielded:
            continue
        yielded.add(key)
        yield candidate


def molecule_leaf_slots(
    node: dict[str, Any],
) -> list[tuple[list[dict[str, Any]], int, dict[str, Any]]]:
    slots: list[tuple[list[dict[str, Any]], int, dict[str, Any]]] = []

    def visit(current: dict[str, Any]) -> None:
        children = current.get("children", []) or []
        if current.get("type") == "mol":
            reaction_children = [
                child
                for child in children
                if isinstance(child, dict) and child.get("type") == "reaction"
            ]
            if not reaction_children:
                return
        for index, child in enumerate(children):
            if not isinstance(child, dict):
                continue
            if child.get("type") == "mol":
                child_reactions = [
                    grandchild
                    for grandchild in child.get("children", []) or []
                    if isinstance(grandchild, dict)
                    and grandchild.get("type") == "reaction"
                ]
                if not child_reactions:
                    slots.append((children, index, child))
                else:
                    visit(child)
            else:
                visit(child)

    visit(node)
    return slots


def collect_reaction_nodes(node: dict[str, Any]) -> list[dict[str, Any]]:
    reactions: list[dict[str, Any]] = []

    def visit(current: dict[str, Any]) -> None:
        if current.get("type") == "reaction":
            reactions.append(current)
        for child in current.get("children", []) or []:
            if isinstance(child, dict):
                visit(child)

    visit(node)
    return reactions


def reattach_original_children(
    generated_root: dict[str, Any],
    original_children: tuple[dict[str, Any], ...],
) -> dict[str, Any]:
    leaf_slots = molecule_leaf_slots(generated_root)
    used_original_indexes: set[int] = set()

    for parent_children, leaf_index, leaf_node in leaf_slots:
        leaf_smiles = molecule_node_smiles(leaf_node)
        matched_index = None
        for original_index, original_child in enumerate(original_children):
            if original_index in used_original_indexes:
                continue
            if same_molecule_smiles(leaf_smiles, molecule_node_smiles(original_child)):
                matched_index = original_index
                break
        if matched_index is None:
            raise RouteSplitError(
                f"generated leaf {leaf_smiles!r} does not match any original reactant"
            )
        used_original_indexes.add(matched_index)
        parent_children[leaf_index] = original_children[matched_index]

    if len(used_original_indexes) != len(original_children):
        unmatched = [
            molecule_node_smiles(child)
            for index, child in enumerate(original_children)
            if index not in used_original_indexes
        ]
        raise RouteSplitError(
            "split route did not regenerate all original reactants: "
            + ", ".join(unmatched)
        )
    return generated_root


def unwrap_split_route(
    plan: SplitPlan,
    ordered_rules: tuple[SingleCenterRule, ...],
    *,
    unwrapper: Callable[..., Any],
) -> dict[str, Any]:
    rule_smarts = [rule.rule_smarts for rule in ordered_rules]
    raw_result = unwrapper(
        plan.parent_smiles,
        rule_smarts,
        route_id=0,
        rule_key_prefix="preprocess_split",
        mark_leaves_in_stock=False,
    )
    if hasattr(raw_result, "routes_json"):
        generated_root = raw_result.routes_json[0]
    elif isinstance(raw_result, dict) and 0 in raw_result:
        generated_root = raw_result[0]
    elif isinstance(raw_result, dict):
        generated_root = raw_result
    else:
        raise TypeError(f"unsupported unwrap result: {type(raw_result)!r}")
    return reattach_original_children(generated_root, plan.original_children)


def original_reactants_smiles(reaction_smiles: str) -> str:
    if ">>" not in reaction_smiles:
        raise RouteSplitError("reaction SMILES does not contain >>")
    reactants, _products = reaction_smiles.split(">>", 1)
    return reactants


def molecule_has_bond(molecule: Any, atom_1: int, atom_2: int) -> bool:
    try:
        return bool(molecule.has_bond(atom_1, atom_2))
    except Exception:
        try:
            molecule.bond(atom_1, atom_2)
        except Exception:
            return False
        return True


def source_molecule_for_broken_rule(reaction: Any, rule: SingleCenterRule) -> Any:
    for molecule in reaction.reactants:
        if all(
            molecule_has_bond(molecule, atom_1, atom_2)
            for atom_1, atom_2 in rule.forward_bonds_broken
        ):
            return molecule
    raise RouteSplitError(
        f"could not find reactant source for broken center {sorted(rule.center_atoms)}"
    )


def restore_broken_rule_on_product(
    molecule: Any,
    source_molecule: Any,
    rule: SingleCenterRule,
) -> Any:
    restored = molecule.copy()
    atoms_to_restore = set(rule.center_atoms)
    for atom_1, atom_2 in rule.forward_bonds_broken:
        atoms_to_restore.add(atom_1)
        atoms_to_restore.add(atom_2)

    for atom_id in sorted(atoms_to_restore):
        try:
            source_atom = source_molecule.atom(atom_id)
        except Exception as exc:
            raise RouteSplitError(
                f"source molecule is missing atom {atom_id} for mapped split"
            ) from exc
        if not restored.has_atom(atom_id):
            restored.add_atom(
                source_atom.copy(),
                atom_id,
                charge=getattr(source_atom, "charge", 0),
                is_radical=getattr(source_atom, "is_radical", False),
                _skip_calculation=True,
            )
        else:
            restored_atom = restored.atom(atom_id)
            restored_atom.charge = getattr(source_atom, "charge", 0)
            restored_atom.is_radical = getattr(source_atom, "is_radical", False)

    for atom_1, atom_2 in rule.forward_bonds_broken:
        if molecule_has_bond(restored, atom_1, atom_2):
            continue
        restored.add_bond(
            atom_1,
            atom_2,
            int(source_molecule.bond(atom_1, atom_2)),
            _skip_calculation=True,
        )

    try:
        restored.fix_structure()
    except Exception:
        pass
    try:
        restored.fix_stereo()
    except Exception:
        pass
    return restored


def split_metadata(
    plan: SplitPlan,
    *,
    step_index: int,
    step_count: int,
    rule_smarts: str,
    forward_change_kind: str,
) -> dict[str, Any]:
    return {
        "route_preprocessing_split": {
            "route_id": str(plan.route_id),
            "original_reaction_id": plan.reaction_id,
            "original_reaction_smiles": plan.original_reaction_smiles,
            "split_step": step_index,
            "split_steps": step_count,
            "rule_smarts": rule_smarts,
            "forward_change_kind": forward_change_kind,
            "split_method": "mapped_intermediate",
        }
    }


def rule_detail(rule: SingleCenterRule) -> dict[str, Any]:
    return {
        "rule_smarts": rule.rule_smarts,
        "center_atoms": sorted(rule.center_atoms),
        "forward_change_kind": rule.forward_change_kind,
        "forward_bonds_formed": [list(bond) for bond in rule.forward_bonds_formed],
        "forward_bonds_broken": [list(bond) for bond in rule.forward_bonds_broken],
        "forward_bonds_changed": [list(bond) for bond in rule.forward_bonds_changed],
    }


def split_reaction_node_via_mapped_intermediates(
    plan: SplitPlan,
    *,
    restore_rules: tuple[SingleCenterRule, ...],
    remaining_rules: tuple[SingleCenterRule, ...],
) -> dict[str, Any]:
    if not restore_rules:
        raise RouteSplitError("mapped fallback has no bond-breaking rule to restore")
    if not remaining_rules and len(plan.original_children) != 1:
        raise RouteSplitError(
            "mapped fallback without a final reaction requires one original reactant"
        )

    from chython import smiles as parse_smiles

    reaction = parse_smiles(plan.original_reaction_smiles)
    if len(reaction.products) != 1:
        raise RouteSplitError("mapped fallback requires one reaction product")

    original_reactants = original_reactants_smiles(plan.original_reaction_smiles)
    active_molecule = reaction.products[0]
    active_smiles = format(active_molecule, "m")
    generated_reaction_smiles: list[str] = []
    ordered_rules = restore_rules + remaining_rules
    step_count = len(restore_rules) + int(bool(remaining_rules))
    first_reaction_node: dict[str, Any] | None = None
    current_mol_node: dict[str, Any] | None = None
    last_child_holder: list[dict[str, Any]] | None = None

    for step_index, rule in enumerate(restore_rules, start=1):
        source_molecule = source_molecule_for_broken_rule(reaction, rule)
        intermediate_molecule = restore_broken_rule_on_product(
            active_molecule,
            source_molecule,
            rule,
        )
        intermediate_smiles = format(intermediate_molecule, "m")
        if same_molecule_smiles(intermediate_smiles, active_smiles):
            raise RouteSplitError("mapped fallback did not change the active molecule")

        child_mol_node = {
            "type": "mol",
            "smiles": intermediate_smiles,
            "mapped_smiles": intermediate_smiles,
            "in_stock": False,
            "children": [],
        }
        reaction_smiles = f"{intermediate_smiles}>>{active_smiles}"
        reaction_node = {
            "type": "reaction",
            "smiles": reaction_smiles,
            "metadata": {
                "smiles": reaction_smiles,
                **split_metadata(
                    plan,
                    step_index=step_index,
                    step_count=step_count,
                    rule_smarts=rule.rule_smarts,
                    forward_change_kind=rule.forward_change_kind,
                ),
            },
            "children": [child_mol_node],
        }
        generated_reaction_smiles.append(reaction_smiles)

        if current_mol_node is None:
            first_reaction_node = reaction_node
        else:
            current_mol_node["children"] = [reaction_node]
        current_mol_node = child_mol_node
        last_child_holder = reaction_node["children"]
        active_molecule = intermediate_molecule
        active_smiles = intermediate_smiles

    if first_reaction_node is None or current_mol_node is None:
        raise RouteSplitError("mapped fallback did not create a split route")

    if remaining_rules:
        rule_smarts = "$".join(rule.rule_smarts for rule in remaining_rules)
        forward_change_kind = "+".join(
            rule.forward_change_kind for rule in remaining_rules
        )
        reaction_smiles = f"{original_reactants}>>{active_smiles}"
        final_reaction_node = {
            "type": "reaction",
            "smiles": reaction_smiles,
            "metadata": {
                "smiles": reaction_smiles,
                **split_metadata(
                    plan,
                    step_index=len(restore_rules) + 1,
                    step_count=step_count,
                    rule_smarts=rule_smarts,
                    forward_change_kind=forward_change_kind,
                ),
            },
            "children": list(plan.original_children),
        }
        generated_reaction_smiles.append(reaction_smiles)
        current_mol_node["children"] = [final_reaction_node]
    else:
        assert last_child_holder is not None
        original_child = plan.original_children[0]
        if not same_molecule_smiles(active_smiles, molecule_node_smiles(original_child)):
            raise RouteSplitError(
                "mapped fallback final intermediate does not match original reactant"
            )
        last_child_holder[0] = original_child

    replace_reaction_node(
        plan.parent_mol_node,
        plan.reaction_node,
        first_reaction_node,
    )
    return {
        "reaction_id": plan.reaction_id,
        "parent_mol_id": plan.parent_mol_id,
        "original_reaction_smiles": plan.original_reaction_smiles,
        "multicenter_rule_smarts": plan.extraction.multicenter_rule_smarts,
        "single_center_rules": [rule.rule_smarts for rule in ordered_rules],
        "single_center_rule_details": [rule_detail(rule) for rule in ordered_rules],
        "single_center_count": len(ordered_rules),
        "protection_related": bool(plan.protection_matches),
        "protection_match_count": len(plan.protection_matches),
        "protection_atom_ids": sorted(plan.protection_atom_ids),
        "generated_reaction_smiles": generated_reaction_smiles,
        "split_method": "mapped_intermediate",
    }


def mapped_intermediate_rule_orders(
    rules: tuple[SingleCenterRule, ...],
) -> Iterable[tuple[tuple[SingleCenterRule, ...], tuple[SingleCenterRule, ...]]]:
    restore_rules = tuple(
        rule
        for rule in rules
        if rule.forward_change_kind == "bond_breaking"
        and rule.forward_bonds_broken
        and not rule.forward_bonds_formed
    )
    remaining_rules = tuple(rule for rule in rules if rule not in restore_rules)
    if not restore_rules:
        return
    if (
        not any(rule.forward_bonds_formed for rule in remaining_rules)
        and remaining_rules
    ):
        return

    yielded: set[tuple[str, ...]] = set()
    restore_orders: Iterable[tuple[SingleCenterRule, ...]]
    if len(restore_rules) <= 6:
        restore_orders = permutations(restore_rules)
    else:
        restore_orders = (restore_rules,)
    for restore_order in restore_orders:
        key = tuple(rule.rule_smarts for rule in restore_order + remaining_rules)
        if key in yielded:
            continue
        yielded.add(key)
        yield tuple(restore_order), remaining_rules


def replace_reaction_node(
    parent_mol_node: dict[str, Any],
    old_reaction_node: dict[str, Any],
    new_reaction_node: dict[str, Any],
) -> None:
    children = parent_mol_node.get("children", []) or []
    for index, child in enumerate(children):
        if child is old_reaction_node:
            children[index] = new_reaction_node
            return
    raise RouteSplitError("original reaction node is no longer attached to its parent")


def split_reaction_node(
    plan: SplitPlan,
    *,
    protection_config: ProtectionAnalysisConfig,
    unwrapper: Callable[..., Any] = unwrap_rule_sequence,
) -> dict[str, Any]:
    if plan.protection_matches and not plan.protection_atom_ids:
        raise RouteSplitError("protection match did not expose mapped atoms")

    last_error: Exception | None = None
    for ordered_rules in split_candidate_rule_orders(
        plan.extraction.single_center_rules,
        plan.protection_atom_ids,
        deprotection_first=protection_config.deprotection_first,
    ):
        try:
            generated_root = unwrap_split_route(
                plan,
                ordered_rules,
                unwrapper=unwrapper,
            )
            first_reactions = [
                child
                for child in generated_root.get("children", []) or []
                if isinstance(child, dict) and child.get("type") == "reaction"
            ]
            if len(first_reactions) != 1:
                raise RouteSplitError("split route did not produce one root reaction")

            generated_reactions = collect_reaction_nodes(generated_root)
            if len(generated_reactions) != len(ordered_rules):
                raise RouteSplitError(
                    "split route reaction count does not match single-center rules"
                )
            for step_index, reaction_node in enumerate(generated_reactions, start=1):
                metadata = reaction_node.setdefault("metadata", {})
                metadata["route_preprocessing_split"] = {
                    "route_id": str(plan.route_id),
                    "original_reaction_id": plan.reaction_id,
                    "original_reaction_smiles": plan.original_reaction_smiles,
                    "split_step": step_index,
                    "split_steps": len(generated_reactions),
                    "rule_smarts": (
                        ordered_rules[step_index - 1].rule_smarts
                        if step_index <= len(ordered_rules)
                        else ""
                    ),
                    "forward_change_kind": (
                        ordered_rules[step_index - 1].forward_change_kind
                        if step_index <= len(ordered_rules)
                        else ""
                    ),
                }

            replace_reaction_node(
                plan.parent_mol_node,
                plan.reaction_node,
                first_reactions[0],
            )
            return {
                "reaction_id": plan.reaction_id,
                "parent_mol_id": plan.parent_mol_id,
                "original_reaction_smiles": plan.original_reaction_smiles,
                "multicenter_rule_smarts": plan.extraction.multicenter_rule_smarts,
                "single_center_rules": [rule.rule_smarts for rule in ordered_rules],
                "single_center_rule_details": [
                    rule_detail(rule) for rule in ordered_rules
                ],
                "single_center_count": len(ordered_rules),
                "protection_related": bool(plan.protection_matches),
                "protection_match_count": len(plan.protection_matches),
                "protection_atom_ids": sorted(plan.protection_atom_ids),
                "generated_reaction_smiles": [
                    reaction_smiles_from_node(reaction)
                    for reaction in generated_reactions
                ],
            }
        except Exception as exc:
            last_error = exc
            continue

    if not plan.protection_matches:
        for restore_rules, remaining_rules in mapped_intermediate_rule_orders(
            plan.extraction.single_center_rules
        ):
            try:
                return split_reaction_node_via_mapped_intermediates(
                    plan,
                    restore_rules=restore_rules,
                    remaining_rules=remaining_rules,
                )
            except Exception as exc:
                last_error = exc
                continue

    raise RouteSplitError(str(last_error) if last_error else "no split order generated")


def unresolved_reaction_record(
    route_id: Any,
    reaction_id: str,
    reaction_smiles: str,
    reason: str,
    extraction: ReactionGranularity | None = None,
    message: str = "",
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "route_id": str(route_id),
        "reaction_id": reaction_id,
        "reaction_smiles": reaction_smiles,
        "reason": reason,
    }
    if extraction is not None:
        record["multicenter_rule_smarts"] = extraction.multicenter_rule_smarts
        record["single_center_rules"] = [
            rule.rule_smarts for rule in extraction.single_center_rules
        ]
        record["single_center_rule_details"] = [
            {
                "rule_smarts": rule.rule_smarts,
                "center_atoms": sorted(rule.center_atoms),
                "forward_change_kind": rule.forward_change_kind,
                "forward_bonds_formed": [
                    list(bond) for bond in rule.forward_bonds_formed
                ],
                "forward_bonds_broken": [
                    list(bond) for bond in rule.forward_bonds_broken
                ],
                "forward_bonds_changed": [
                    list(bond) for bond in rule.forward_bonds_changed
                ],
            }
            for rule in extraction.single_center_rules
        ]
        record["center_components"] = [
            sorted(component) for component in extraction.center_components
        ]
    if message:
        record["message"] = message
    return record


def error_record(route_id: Any, stage: str, exc: Exception) -> dict[str, Any]:
    return {
        "route_id": str(route_id),
        "stage": stage,
        "error_type": type(exc).__qualname__,
        "message": str(exc) or traceback.format_exc(limit=1).strip(),
    }


def annotate_route(
    route: dict[str, Any],
    *,
    route_id: Any,
    status: str,
    changes: list[dict[str, Any]],
    unresolved_reactions: list[dict[str, Any]],
    errors: list[dict[str, Any]],
) -> None:
    metadata = route.setdefault("metadata", {})
    metadata["route_preprocessing"] = {
        "route_id": str(route_id),
        "status": status,
        "changes": changes,
        "unresolved_reactions": unresolved_reactions,
        "errors": errors,
    }


def process_route(
    route_id: Any,
    route: dict[str, Any],
    *,
    extractor: Any,
    protection_rules: dict[str, ProtectionRule],
    protection_config: ProtectionAnalysisConfig,
    normalizer: Callable[[dict[str, Any]], dict[str, Any]] = normalize_route_tree,
    protection_detector: Callable[
        [RouteIndex, ReactionRecord, dict[str, ProtectionRule], ProtectionAnalysisConfig],
        list[Any],
    ] = detect_deprotections,
    protection_event_analyzer: Callable[..., tuple[list[Any], list[Any], RouteIndex]] = (
        analyze_route_protection
    ),
    unwrapper: Callable[..., Any] = unwrap_rule_sequence,
    ignore_errors: bool = False,
) -> RoutePreprocessResult:
    errors: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    split_plans: list[SplitPlan] = []
    multicenter_reactions = 0
    protection_multicenter_reactions = 0

    try:
        normalized_route = normalizer(route)
    except Exception as exc:
        if not ignore_errors:
            raise
        errors.append(error_record(route_id, "normalize_route", exc))
        route_copy = copy.deepcopy(route)
        annotate_route(
            route_copy,
            route_id=route_id,
            status="error",
            changes=[],
            unresolved_reactions=[],
            errors=errors,
        )
        return RoutePreprocessResult(
            route_id=route_id,
            route=route_copy,
            normalized=False,
            modified=False,
            multicenter_reactions=0,
            protection_multicenter_reactions=0,
            split_reactions=0,
            protection_split_reactions=0,
            changes=[],
            unresolved_reactions=[],
            errors=errors,
        )

    try:
        index = build_route_index(normalized_route)
    except Exception as exc:
        if not ignore_errors:
            raise
        errors.append(error_record(route_id, "build_route_index", exc))
        annotate_route(
            normalized_route,
            route_id=route_id,
            status="error",
            changes=[],
            unresolved_reactions=[],
            errors=errors,
        )
        return RoutePreprocessResult(
            route_id=route_id,
            route=normalized_route,
            normalized=True,
            modified=False,
            multicenter_reactions=0,
            protection_multicenter_reactions=0,
            split_reactions=0,
            protection_split_reactions=0,
            changes=[],
            unresolved_reactions=[],
            errors=errors,
        )

    protection_events_by_protection_node: dict[str, tuple[Any, ...]] | None = None

    def protection_events_for_node(reaction_id: str) -> tuple[Any, ...]:
        nonlocal protection_events_by_protection_node
        if protection_events_by_protection_node is None:
            event_config = replace(protection_config, collect_interval_rules=False)
            events, _interval_rules, _event_index = protection_event_analyzer(
                index.route,
                route_id,
                protection_rules,
                config=event_config,
            )
            grouped: dict[str, list[Any]] = {}
            for event in events:
                protection_node_id = getattr(event, "protection_node_id", "")
                if not protection_node_id:
                    continue
                grouped.setdefault(str(protection_node_id), []).append(event)
            protection_events_by_protection_node = {
                node_id: tuple(node_events)
                for node_id, node_events in grouped.items()
            }
        return protection_events_by_protection_node.get(reaction_id, ())

    for reaction_id in index.reaction_order:
        rxn_record = index.reaction_records[reaction_id]
        try:
            extraction = extractor.extract(rxn_record.reaction_smiles)
        except Exception as exc:
            if not ignore_errors:
                raise
            errors.append(error_record(route_id, "extract_reaction_rules", exc))
            unresolved.append(
                unresolved_reaction_record(
                    route_id,
                    reaction_id,
                    rxn_record.reaction_smiles,
                    "rule_extraction_error",
                    message=str(exc),
                )
            )
            continue

        if extraction.skipped or not extraction.is_multicenter:
            continue

        multicenter_reactions += 1
        try:
            matches = tuple(
                protection_detector(
                    index,
                    rxn_record,
                    protection_rules,
                    protection_config,
                )
            )
        except Exception as exc:
            if not ignore_errors:
                raise
            errors.append(error_record(route_id, "detect_deprotections", exc))
            matches = ()

        if not matches:
            try:
                matches = protection_events_for_node(reaction_id)
            except Exception as exc:
                if not ignore_errors:
                    raise
                errors.append(error_record(route_id, "analyze_route_protection", exc))
                matches = ()

        protection_atom_ids = frozenset()
        if matches:
            protection_atom_ids = protection_atom_ids_from_matches(matches)
            if not any(
                rule.center_atoms & protection_atom_ids
                for rule in extraction.single_center_rules
            ):
                unresolved.append(
                    unresolved_reaction_record(
                        route_id,
                        reaction_id,
                        rxn_record.reaction_smiles,
                        "protection_center_not_matched_to_single_center_rule",
                        extraction,
                    )
                )
                continue
            protection_multicenter_reactions += 1

        parent_mol_id = index.parent_mol_by_reaction[reaction_id]
        parent_mol_node = index.molecule_records[parent_mol_id].node
        original_children = tuple(
            child
            for child in rxn_record.node.get("children", []) or []
            if isinstance(child, dict) and child.get("type") == "mol"
        )
        split_plans.append(
            SplitPlan(
                route_id=route_id,
                reaction_id=reaction_id,
                parent_mol_id=parent_mol_id,
                reaction_node=rxn_record.node,
                parent_mol_node=parent_mol_node,
                original_children=original_children,
                original_reaction_smiles=rxn_record.reaction_smiles,
                parent_smiles=index.molecule_records[parent_mol_id].smiles,
                extraction=extraction,
                protection_matches=matches,
                protection_atom_ids=protection_atom_ids,
            )
        )

    if unresolved:
        annotate_route(
            index.route,
            route_id=route_id,
            status="unresolved",
            changes=[],
            unresolved_reactions=unresolved,
            errors=errors,
        )
        return RoutePreprocessResult(
            route_id=route_id,
            route=index.route,
            normalized=True,
            modified=False,
            multicenter_reactions=multicenter_reactions,
            protection_multicenter_reactions=protection_multicenter_reactions,
            split_reactions=0,
            protection_split_reactions=0,
            changes=[],
            unresolved_reactions=unresolved,
            errors=errors,
        )

    changes: list[dict[str, Any]] = []
    protection_split_reactions = 0
    for plan in split_plans:
        try:
            change = split_reaction_node(
                plan,
                protection_config=protection_config,
                unwrapper=unwrapper,
            )
            changes.append(change)
            protection_split_reactions += int(bool(change.get("protection_related")))
        except Exception as exc:
            if not ignore_errors:
                raise
            errors.append(error_record(route_id, "split_multicenter_reaction", exc))
            reason = (
                "protection_related_split_failed"
                if plan.protection_matches
                else "multicenter_split_failed"
            )
            unresolved.append(
                unresolved_reaction_record(
                    route_id,
                    plan.reaction_id,
                    plan.original_reaction_smiles,
                    reason,
                    plan.extraction,
                    message=str(exc),
                )
            )

    if unresolved:
        annotate_route(
            normalized_route,
            route_id=route_id,
            status="unresolved",
            changes=[],
            unresolved_reactions=unresolved,
            errors=errors,
        )
        return RoutePreprocessResult(
            route_id=route_id,
            route=normalized_route,
            normalized=True,
            modified=False,
            multicenter_reactions=multicenter_reactions,
            protection_multicenter_reactions=protection_multicenter_reactions,
            split_reactions=0,
            protection_split_reactions=0,
            changes=[],
            unresolved_reactions=unresolved,
            errors=errors,
        )

    if changes:
        try:
            final_route = normalizer(index.route)
        except Exception as exc:
            if not ignore_errors:
                raise
            errors.append(error_record(route_id, "normalize_split_route", exc))
            unresolved = [
                unresolved_reaction_record(
                    route_id,
                    plan.reaction_id,
                    plan.original_reaction_smiles,
                    "split_route_normalization_failed",
                    plan.extraction,
                    message=str(exc),
                )
                for plan in split_plans
            ]
            annotate_route(
                normalized_route,
                route_id=route_id,
                status="unresolved",
                changes=[],
                unresolved_reactions=unresolved,
                errors=errors,
            )
            return RoutePreprocessResult(
                route_id=route_id,
                route=normalized_route,
                normalized=True,
                modified=False,
                multicenter_reactions=multicenter_reactions,
                protection_multicenter_reactions=protection_multicenter_reactions,
                split_reactions=0,
                protection_split_reactions=0,
                changes=[],
                unresolved_reactions=unresolved,
                errors=errors,
            )
        annotate_route(
            final_route,
            route_id=route_id,
            status="modified",
            changes=changes,
            unresolved_reactions=[],
            errors=errors,
        )
        return RoutePreprocessResult(
            route_id=route_id,
            route=final_route,
            normalized=True,
            modified=True,
            multicenter_reactions=multicenter_reactions,
            protection_multicenter_reactions=protection_multicenter_reactions,
            split_reactions=len(changes),
            protection_split_reactions=protection_split_reactions,
            changes=changes,
            unresolved_reactions=[],
            errors=errors,
        )

    return RoutePreprocessResult(
        route_id=route_id,
        route=index.route,
        normalized=True,
        modified=False,
        multicenter_reactions=multicenter_reactions,
        protection_multicenter_reactions=protection_multicenter_reactions,
        split_reactions=0,
        protection_split_reactions=0,
        changes=[],
        unresolved_reactions=[],
        errors=errors,
    )


def empty_collection_like(routes_json: Any) -> Any:
    if isinstance(routes_json, list):
        return []
    if isinstance(routes_json, dict):
        return {}
    raise TypeError(f"unsupported routes JSON root: {type(routes_json)!r}")


def add_route_to_collection(collection: Any, route_id: Any, route: dict[str, Any]) -> None:
    if isinstance(collection, list):
        collection.append(route)
        return
    collection[str(route_id)] = route


def route_id_sort_key(value: Any) -> tuple[int, Any]:
    if isinstance(value, int):
        return (0, value)
    if isinstance(value, str) and value.isdigit():
        return (0, int(value))
    return (1, str(value))


def preprocess_routes_json(
    routes_json: Any,
    *,
    extractor: Any,
    protection_rules: dict[str, ProtectionRule],
    protection_config: ProtectionAnalysisConfig,
    normalizer: Callable[[dict[str, Any]], dict[str, Any]] = normalize_route_tree,
    protection_detector: Callable[
        [RouteIndex, ReactionRecord, dict[str, ProtectionRule], ProtectionAnalysisConfig],
        list[Any],
    ] = detect_deprotections,
    protection_event_analyzer: Callable[..., tuple[list[Any], list[Any], RouteIndex]] = (
        analyze_route_protection
    ),
    unwrapper: Callable[..., Any] = unwrap_rule_sequence,
    ignore_errors: bool = False,
    limit: int | None = None,
    progress_interval: int = 0,
) -> tuple[Any, dict[str, dict[str, Any]], dict[str, dict[str, Any]], dict[str, Any]]:
    cleaned_routes = empty_collection_like(routes_json)
    resolved_routes: dict[str, dict[str, Any]] = {}
    unresolved_routes: dict[str, dict[str, Any]] = {}
    errors: list[dict[str, Any]] = []
    modified_route_ids: list[str] = []
    unresolved_route_ids: list[str] = []
    all_changes: list[dict[str, Any]] = []
    total_routes = 0
    normalized_routes = 0
    multicenter_reactions = 0
    protection_multicenter_reactions = 0
    split_reactions = 0
    protection_split_reactions = 0

    for index, (route_id, route) in enumerate(route_items(routes_json), start=1):
        if limit is not None and index > limit:
            break
        result = process_route(
            route_id,
            route,
            extractor=extractor,
            protection_rules=protection_rules,
            protection_config=protection_config,
            normalizer=normalizer,
            protection_detector=protection_detector,
            protection_event_analyzer=protection_event_analyzer,
            unwrapper=unwrapper,
            ignore_errors=ignore_errors,
        )
        total_routes += 1
        normalized_routes += int(result.normalized)
        multicenter_reactions += result.multicenter_reactions
        protection_multicenter_reactions += result.protection_multicenter_reactions
        split_reactions += result.split_reactions
        protection_split_reactions += result.protection_split_reactions
        errors.extend(result.errors)
        add_route_to_collection(cleaned_routes, route_id, result.route)

        route_id_text = str(route_id)
        if result.modified:
            modified_route_ids.append(route_id_text)
            resolved_routes[route_id_text] = result.route
            all_changes.extend(
                {"route_id": route_id_text, **change} for change in result.changes
            )
        if result.unresolved_reactions:
            unresolved_route_ids.append(route_id_text)
            unresolved_routes[route_id_text] = result.route

        if progress_interval and total_routes % progress_interval == 0:
            print(
                "[preprocess-routes] processed "
                f"{total_routes} routes; modified={len(modified_route_ids)}; "
                f"unresolved={len(unresolved_route_ids)}; errors={len(errors)}",
                file=sys.stderr,
                flush=True,
            )

    summary = {
        "total_routes_processed": total_routes,
        "number_of_normalized_routes": normalized_routes,
        "number_of_routes_modified": len(modified_route_ids),
        "number_of_multicenter_reactions_found": multicenter_reactions,
        "number_of_multicenter_reactions_split": split_reactions,
        "number_of_protection_related_multicenter_reactions_split": (
            protection_split_reactions
        ),
        "number_of_non_protection_multicenter_reactions_split": (
            split_reactions - protection_split_reactions
        ),
        "number_of_protection_related_multicenter_reactions_found": (
            protection_multicenter_reactions
        ),
        "number_of_unresolved_multicenter_routes": len(unresolved_route_ids),
        "route_ids_for_modified_routes": sorted(
            modified_route_ids,
            key=route_id_sort_key,
        ),
        "route_ids_for_unresolved_routes": sorted(
            unresolved_route_ids,
            key=route_id_sort_key,
        ),
        "changes": all_changes,
        "errors": errors,
        "number_of_errors": len(errors),
    }
    return cleaned_routes, resolved_routes, unresolved_routes, summary


def dataset_output_paths(
    output_dir: Path,
    dataset_name: str,
    *,
    sidecar_dir: Path | None = None,
) -> dict[str, Path]:
    output_path = output_dir / dataset_name
    report_dir = sidecar_dir or output_dir
    stem = output_path.stem
    return {
        "cleaned": output_path,
        "resolved": report_dir / f"{stem}_multicenter_resolved.json",
        "unresolved": report_dir / f"{stem}_multicenter_unresolved.json",
        "summary": report_dir / f"{stem}_preprocess_summary.json",
    }


def resolve_dataset_path(input_dir: Path, dataset: str | Path) -> Path:
    dataset_path = Path(dataset)
    candidates = []
    if dataset_path.is_absolute():
        candidates.append(dataset_path)
    else:
        candidates.append(input_dir / dataset_path)
        name = dataset_path.name
        candidates.append(input_dir / name.replace("_", "-"))
        candidates.append(input_dir / name.replace("-", "_"))

    for candidate in candidates:
        resolved = resolve_existing_path(candidate)
        if resolved.exists():
            return resolved
    return candidates[0]


def preprocess_routes_file(
    input_path: Path,
    output_paths: dict[str, Path],
    *,
    extractor: Any,
    protection_rules: dict[str, ProtectionRule],
    protection_config: ProtectionAnalysisConfig,
    ignore_errors: bool,
    limit: int | None = None,
    progress_interval: int = 0,
    normalizer: Callable[[dict[str, Any]], dict[str, Any]] = normalize_route_tree,
    protection_detector: Callable[
        [RouteIndex, ReactionRecord, dict[str, ProtectionRule], ProtectionAnalysisConfig],
        list[Any],
    ] = detect_deprotections,
    protection_event_analyzer: Callable[..., tuple[list[Any], list[Any], RouteIndex]] = (
        analyze_route_protection
    ),
    unwrapper: Callable[..., Any] = unwrap_rule_sequence,
) -> dict[str, Any]:
    routes_json = read_json(input_path)
    cleaned, resolved, unresolved, summary = preprocess_routes_json(
        routes_json,
        extractor=extractor,
        protection_rules=protection_rules,
        protection_config=protection_config,
        ignore_errors=ignore_errors,
        limit=limit,
        progress_interval=progress_interval,
        normalizer=normalizer,
        protection_detector=protection_detector,
        protection_event_analyzer=protection_event_analyzer,
        unwrapper=unwrapper,
    )
    summary = {
        "input_file": str(input_path),
        "output_files": {name: str(path) for name, path in output_paths.items()},
        **summary,
    }

    write_json(output_paths["cleaned"], cleaned)
    write_json(output_paths["resolved"], resolved)
    write_json(output_paths["unresolved"], unresolved)
    write_json(output_paths["summary"], summary)
    return summary


def preprocess_datasets(args: argparse.Namespace) -> dict[str, Any]:
    setup_runtime_cache_dirs()
    input_dir = resolve_existing_path(args.input_dir)
    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_root = (
        Path(args.summary_dir).expanduser()
        if getattr(args, "summary_dir", None)
        else None
    )
    protection_config_path = (
        resolve_existing_path(args.protection_config)
        if getattr(args, "protection_config", None)
        else None
    )
    protection_config = ProtectionAnalysisConfig.from_yaml(protection_config_path)
    protection_config.collect_interval_rules = False
    if getattr(args, "ignore_errors", False):
        protection_config.ignore_errors = True
    extractor = SynPlannerGranularityExtractor.from_args(args)
    protection_rules = load_chython_protection_rules()

    dataset_summaries = {}
    for dataset in args.datasets:
        input_path = resolve_dataset_path(input_dir, dataset)
        stage_dir = (
            stage_output_dir(
                summary_root,
                dataset_prefix_from_path(dataset),
                "preprocess",
            )
            if summary_root is not None
            else None
        )
        output_paths = dataset_output_paths(
            output_dir,
            Path(dataset).name,
            sidecar_dir=stage_dir,
        )
        print(
            f"[preprocess-routes] processing {input_path} -> {output_paths['cleaned']}",
            file=sys.stderr,
            flush=True,
        )
        dataset_summary = preprocess_routes_file(
            input_path,
            output_paths,
            extractor=extractor,
            protection_rules=protection_rules,
            protection_config=protection_config,
            ignore_errors=getattr(args, "ignore_errors", False),
            limit=getattr(args, "limit", None),
            progress_interval=getattr(args, "progress_interval", 0),
        )
        write_standard_sidecars(
            output_paths["summary"].parent,
            command_name="preprocess-routes",
            summary=dataset_summary,
            errors=dataset_summary.get("errors", []),
            input_files=[input_path],
            output_files=dataset_summary["output_files"],
            config_path=getattr(args, "config", None),
            cli_args=args,
        )
        dataset_summaries[Path(dataset).name] = dataset_summary

    aggregate = {
        "datasets": dataset_summaries,
        "total_routes_processed": sum(
            summary["total_routes_processed"]
            for summary in dataset_summaries.values()
        ),
        "number_of_normalized_routes": sum(
            summary["number_of_normalized_routes"]
            for summary in dataset_summaries.values()
        ),
        "number_of_routes_modified": sum(
            summary["number_of_routes_modified"]
            for summary in dataset_summaries.values()
        ),
        "number_of_multicenter_reactions_found": sum(
            summary["number_of_multicenter_reactions_found"]
            for summary in dataset_summaries.values()
        ),
        "number_of_multicenter_reactions_split": sum(
            summary["number_of_multicenter_reactions_split"]
            for summary in dataset_summaries.values()
        ),
        "number_of_protection_related_multicenter_reactions_split": sum(
            summary["number_of_protection_related_multicenter_reactions_split"]
            for summary in dataset_summaries.values()
        ),
        "number_of_non_protection_multicenter_reactions_split": sum(
            summary["number_of_non_protection_multicenter_reactions_split"]
            for summary in dataset_summaries.values()
        ),
        "number_of_unresolved_multicenter_routes": sum(
            summary["number_of_unresolved_multicenter_routes"]
            for summary in dataset_summaries.values()
        ),
        "number_of_errors": sum(
            summary["number_of_errors"] for summary in dataset_summaries.values()
        ),
    }
    aggregate_path = (summary_root or output_dir) / "preprocess_routes_summary.json"
    aggregate["summary_file"] = str(aggregate_path)
    write_json(aggregate_path, aggregate)
    return aggregate


def run(args: argparse.Namespace) -> int:
    if normalize_n_cpu(getattr(args, "n_cpu", 1)) != 1:
        print(
            "[preprocess-routes] --n-cpu is accepted for CLI consistency, "
            "but preprocessing currently runs sequentially.",
            file=sys.stderr,
            flush=True,
        )
    summary = preprocess_datasets(args)
    print(json.dumps(summary, indent=2), flush=True)
    return 0
