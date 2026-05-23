import copy
import json
from dataclasses import dataclass

from route_inspector.preprocess_routes import (
    ReactionGranularity,
    SingleCenterRule,
    dataset_output_paths,
    preprocess_routes_file,
    preprocess_routes_json,
)
from route_inspector.protection.analysis import ProtectionAnalysisConfig


PARENT = "[CH2:1]=[O:2]"
CHILD = "[CH3:1][OH:2]"
REACTION = f"{CHILD}>>{PARENT}"


class FakeExtractor:
    def __init__(self, extraction):
        self.extraction = extraction
        self.reactions = []

    def extract(self, reaction_smiles):
        self.reactions.append(reaction_smiles)
        return self.extraction


@dataclass(frozen=True)
class FakeProtectionMatch:
    protected_atom_ids: tuple[int, ...] = (2,)
    raw_mapping: tuple[tuple[int, int], ...] = ((1, 2),)


@dataclass(frozen=True)
class FakeProtectionEvent:
    protection_node_id: str = "r0"
    protected_atom_ids: tuple[int, ...] = (2,)


def route_with_one_reaction():
    return {
        "smiles": PARENT,
        "type": "mol",
        "in_stock": False,
        "children": [
            {
                "type": "reaction",
                "metadata": {"smiles": REACTION},
                "children": [
                    {
                        "smiles": CHILD,
                        "type": "mol",
                        "in_stock": True,
                        "children": [],
                    }
                ],
            }
        ],
    }


def multicenter_extraction():
    return ReactionGranularity(
        reaction_smiles=REACTION,
        multicenter_rule_smarts="deprotect$transform",
        single_center_rules=(
            SingleCenterRule(
                "deprotect",
                frozenset({2}),
                forward_change_kind="bond_breaking",
                forward_bonds_broken=((2, 3),),
            ),
            SingleCenterRule(
                "transform",
                frozenset({1}),
                forward_change_kind="bond_forming",
                forward_bonds_formed=((1, 2),),
            ),
        ),
        center_components=(frozenset({2}), frozenset({1})),
    )


def non_protection_multicenter_extraction():
    return ReactionGranularity(
        reaction_smiles=REACTION,
        multicenter_rule_smarts="form$break",
        single_center_rules=(
            SingleCenterRule(
                "form",
                frozenset({1}),
                forward_change_kind="bond_forming",
                forward_bonds_formed=((1, 2),),
            ),
            SingleCenterRule(
                "break",
                frozenset({2}),
                forward_change_kind="bond_breaking",
                forward_bonds_broken=((2, 3),),
            ),
        ),
        center_components=(frozenset({1}), frozenset({2})),
    )


def identity_normalizer(route):
    return copy.deepcopy(route)


def test_route_normalization_is_called():
    calls = []

    def normalizer(route):
        calls.append(route)
        return copy.deepcopy(route)

    cleaned, _resolved, _unresolved, summary = preprocess_routes_json(
        [{"smiles": "CCO", "type": "mol", "in_stock": False, "children": []}],
        extractor=FakeExtractor(
            ReactionGranularity("", "", (), (), skipped=True)
        ),
        protection_rules={},
        protection_config=ProtectionAnalysisConfig(collect_interval_rules=False),
        normalizer=normalizer,
    )

    assert len(calls) == 1
    assert summary["number_of_normalized_routes"] == 1
    assert cleaned[0]["smiles"] == "CCO"


def test_non_protection_multicenter_reaction_is_split_by_forward_bond_order():
    seen_rule_orders = []

    def fake_unwrapper(target_smiles, rule_smarts, **_kwargs):
        seen_rule_orders.append(tuple(rule_smarts))
        return {
            0: {
                "type": "mol",
                "smiles": target_smiles,
                "children": [
                    {
                        "type": "reaction",
                        "smiles": "intermediate>>target",
                        "children": [
                            {
                                "type": "mol",
                                "smiles": "intermediate",
                                "children": [
                                    {
                                        "type": "reaction",
                                        "smiles": "child>>intermediate",
                                        "children": [
                                            {
                                                "type": "mol",
                                                "smiles": CHILD,
                                                "children": [],
                                            }
                                        ],
                                    }
                                ],
                            }
                        ],
                    }
                ],
            }
        }

    cleaned, resolved, unresolved, summary = preprocess_routes_json(
        [route_with_one_reaction()],
        extractor=FakeExtractor(non_protection_multicenter_extraction()),
        protection_rules={},
        protection_config=ProtectionAnalysisConfig(collect_interval_rules=False),
        normalizer=identity_normalizer,
        protection_detector=lambda *_args: [],
        unwrapper=fake_unwrapper,
    )

    root_reaction = cleaned[0]["children"][0]
    nested_reaction = root_reaction["children"][0]["children"][0]

    assert seen_rule_orders[0] == ("break", "form")
    assert root_reaction["smiles"] == "intermediate>>target"
    assert nested_reaction["smiles"] == "child>>intermediate"
    assert list(resolved) == ["0"]
    assert unresolved == {}
    assert summary["number_of_multicenter_reactions_found"] == 1
    assert summary["number_of_multicenter_reactions_split"] == 1
    assert summary["number_of_protection_related_multicenter_reactions_split"] == 0
    assert summary["number_of_non_protection_multicenter_reactions_split"] == 1


def test_non_protection_mapped_intermediate_fallback_splits_n_oxide_chlorination():
    product = (
        "[c:14]1([cH:31][cH:32][c:6]([Br:69])[c:7]([Cl:42])[n:13]1)"
        "[C:15](=[O:38])[O:68][CH3:67]"
    )
    n_oxide = (
        "[cH:32]1[cH:31][c:14]([n+:13]([cH:7][c:6]1[Br:69])[O-:77])"
        "[C:15](=[O:38])[O:68][CH3:67]"
    )
    reagent = "[Cl:42][P:74]([Cl:75])([Cl:76])=[O:73]"
    reaction = f"{reagent}.{n_oxide}>>{product}"
    route = {
        "smiles": product,
        "type": "mol",
        "in_stock": False,
        "children": [
            {
                "type": "reaction",
                "metadata": {"smiles": reaction},
                "children": [
                    {"smiles": reagent, "type": "mol", "in_stock": True},
                    {"smiles": n_oxide, "type": "mol", "in_stock": True},
                ],
            }
        ],
    }
    extraction = ReactionGranularity(
        reaction_smiles=reaction,
        multicenter_rule_smarts="chlorination$n_oxide_cleavage",
        single_center_rules=(
            SingleCenterRule(
                "chlorination",
                frozenset({7, 42, 74}),
                forward_change_kind="bond_forming_and_breaking",
                forward_bonds_formed=((7, 42),),
                forward_bonds_broken=((42, 74),),
            ),
            SingleCenterRule(
                "n_oxide_cleavage",
                frozenset({13, 77}),
                forward_change_kind="bond_breaking",
                forward_bonds_broken=((13, 77),),
            ),
        ),
        center_components=(frozenset({7, 42, 74}), frozenset({13, 77})),
    )

    def failing_unwrapper(*_args, **_kwargs):
        raise RuntimeError("force mapped fallback")

    cleaned, resolved, unresolved, summary = preprocess_routes_json(
        [route],
        extractor=FakeExtractor(extraction),
        protection_rules={},
        protection_config=ProtectionAnalysisConfig(collect_interval_rules=False),
        normalizer=identity_normalizer,
        protection_detector=lambda *_args: [],
        unwrapper=failing_unwrapper,
        ignore_errors=True,
    )

    first_reaction = cleaned[0]["children"][0]
    intermediate = first_reaction["children"][0]
    second_reaction = intermediate["children"][0]
    change = summary["changes"][0]

    assert "[n+:13]([O-:77])" in first_reaction["smiles"].split(">>", 1)[0]
    assert first_reaction["smiles"].endswith(f">>{product}")
    assert second_reaction["smiles"].startswith(f"{reagent}.{n_oxide}>>")
    assert "[n+:13]([O-:77])" in second_reaction["smiles"].split(">>", 1)[1]
    assert (
        change["single_center_rule_details"][0]["forward_change_kind"]
        == "bond_breaking"
    )
    assert (
        change["single_center_rule_details"][1]["forward_change_kind"]
        == "bond_forming_and_breaking"
    )
    assert change["split_method"] == "mapped_intermediate"
    assert list(resolved) == ["0"]
    assert unresolved == {}
    assert summary["number_of_non_protection_multicenter_reactions_split"] == 1


def test_protection_related_multicenter_reaction_is_split():
    seen_rule_orders = []

    def fake_unwrapper(target_smiles, rule_smarts, **_kwargs):
        seen_rule_orders.append(tuple(rule_smarts))
        return {
            0: {
                "type": "mol",
                "smiles": target_smiles,
                "children": [
                    {
                        "type": "reaction",
                        "smiles": "intermediate>>target",
                        "children": [
                            {
                                "type": "mol",
                                "smiles": "intermediate",
                                "children": [
                                    {
                                        "type": "reaction",
                                        "smiles": "child>>intermediate",
                                        "children": [
                                            {
                                                "type": "mol",
                                                "smiles": CHILD,
                                                "children": [],
                                            }
                                        ],
                                    }
                                ],
                            }
                        ],
                    }
                ],
            }
        }

    cleaned, resolved, unresolved, summary = preprocess_routes_json(
        [route_with_one_reaction()],
        extractor=FakeExtractor(multicenter_extraction()),
        protection_rules={},
        protection_config=ProtectionAnalysisConfig(collect_interval_rules=False),
        normalizer=identity_normalizer,
        protection_detector=lambda *_args: [FakeProtectionMatch()],
        unwrapper=fake_unwrapper,
    )

    root_reaction = cleaned[0]["children"][0]
    nested_reaction = root_reaction["children"][0]["children"][0]

    assert seen_rule_orders[0] == ("deprotect", "transform")
    assert root_reaction["smiles"] == "intermediate>>target"
    assert nested_reaction["smiles"] == "child>>intermediate"
    assert nested_reaction["children"][0]["smiles"] == CHILD
    assert list(resolved) == ["0"]
    assert unresolved == {}
    assert summary["number_of_routes_modified"] == 1
    assert (
        summary["number_of_protection_related_multicenter_reactions_split"] == 1
    )


def test_protection_introduction_multicenter_reaction_is_split_from_trace_event():
    def fake_unwrapper(target_smiles, rule_smarts, **_kwargs):
        assert tuple(rule_smarts) == ("deprotect", "transform")
        return {
            0: {
                "type": "mol",
                "smiles": target_smiles,
                "children": [
                    {
                        "type": "reaction",
                        "smiles": "protected>>target",
                        "children": [
                            {
                                "type": "mol",
                                "smiles": "deprotected",
                                "children": [
                                    {
                                        "type": "reaction",
                                        "smiles": "child>>protected",
                                        "children": [
                                            {
                                                "type": "mol",
                                                "smiles": CHILD,
                                                "children": [],
                                            }
                                        ],
                                    }
                                ],
                            }
                        ],
                    }
                ],
            }
        }

    _cleaned, resolved, unresolved, summary = preprocess_routes_json(
        [route_with_one_reaction()],
        extractor=FakeExtractor(multicenter_extraction()),
        protection_rules={},
        protection_config=ProtectionAnalysisConfig(collect_interval_rules=False),
        normalizer=identity_normalizer,
        protection_detector=lambda *_args: [],
        protection_event_analyzer=lambda *_args, **_kwargs: (
            [FakeProtectionEvent()],
            [],
            None,
        ),
        unwrapper=fake_unwrapper,
    )

    assert list(resolved) == ["0"]
    assert unresolved == {}
    assert summary["number_of_routes_modified"] == 1


def test_unresolved_routes_and_summary_are_written(tmp_path):
    def failing_unwrapper(*_args, **_kwargs):
        raise RuntimeError("cannot split")

    input_path = tmp_path / "n1_routes.json"
    input_path.write_text(json.dumps([route_with_one_reaction()]), encoding="utf-8")
    output_paths = dataset_output_paths(tmp_path / "clean", "n1_routes.json")

    summary = preprocess_routes_file(
        input_path,
        output_paths,
        extractor=FakeExtractor(multicenter_extraction()),
        protection_rules={},
        protection_config=ProtectionAnalysisConfig(collect_interval_rules=False),
        ignore_errors=True,
        normalizer=identity_normalizer,
        protection_detector=lambda *_args: [],
        unwrapper=failing_unwrapper,
    )

    cleaned = json.loads(output_paths["cleaned"].read_text())
    unresolved = json.loads(output_paths["unresolved"].read_text())
    resolved = json.loads(output_paths["resolved"].read_text())
    written_summary = json.loads(output_paths["summary"].read_text())

    assert cleaned[0]["metadata"]["route_preprocessing"]["status"] == "unresolved"
    assert list(unresolved) == ["0"]
    assert resolved == {}
    assert summary["number_of_unresolved_multicenter_routes"] == 1
    assert written_summary["total_routes_processed"] == 1
    assert output_paths["summary"].exists()
