from chython import smarts

from route_inspector.protection.analysis import (
    ProtectionAnalysisResult,
    ProtectionSingleRuleObservation,
    ProtectionAnalysisConfig,
    analyze_protection_in_routes,
    analyze_route_protection,
    summarize_aggregate_single_rules,
    summarize_single_rules,
)
from route_inspector.protection.chython_rules import ProtectionRule
from route_inspector.protection.outputs import write_protection_outputs


PROTECTED = "[CH3:1][CH2:2][O:3][Si:4]([CH3:5])([CH3:6])[CH3:7]"
UNPROTECTED = "[CH3:1][CH2:2][OH:3]"
DEP_RXN = f"{PROTECTED}>>{UNPROTECTED}"
PROT_RXN = f"{UNPROTECTED}>>{PROTECTED}"


def tms_rule() -> ProtectionRule:
    return ProtectionRule(
        rule_id="hydroxyl_tms",
        name="hydroxyl_tms",
        query=smarts("[O;D2:1]-[Si:2]([C:3])([C:4])[C:5]"),
        product_query=smarts("[A:1]"),
        product_smarts="[A:1]",
        atoms_to_keep=(1,),
        atoms_to_add=(),
        protected_example="CCO[Si](C)(C)C",
        cleaved_example="CCO",
        decoys=(),
        source="test",
    )


def route_with_child(child):
    return {
        "smiles": "CCO",
        "type": "mol",
        "in_stock": False,
        "children": [
            {
                "type": "reaction",
                "metadata": {"smiles": DEP_RXN},
                "children": [child],
            }
        ],
    }


def protected_mol(*, in_stock=False, children=None):
    return {
        "smiles": "CCO[Si](C)(C)C",
        "type": "mol",
        "in_stock": in_stock,
        "children": children or [],
    }


def unprotected_stock():
    return {"smiles": "CCO", "type": "mol", "in_stock": True, "children": []}


def test_detects_stock_derived_deprotection():
    events, interval_rules, _index = analyze_route_protection(
        route_with_child(protected_mol(in_stock=True)),
        "r0",
        {"hydroxyl_tms": tms_rule()},
        config=ProtectionAnalysisConfig(collect_interval_rules=False),
    )

    assert len(events) == 1
    assert events[0].pg_type == "hydroxyl_tms"
    assert events[0].trace_status == "stock"
    assert events[0].source_type == "stock"
    assert events[0].stock_node_id
    assert interval_rules == []


def test_detects_route_introduced_protecting_group():
    route = route_with_child(
        protected_mol(
            children=[
                {
                    "type": "reaction",
                    "metadata": {"smiles": PROT_RXN},
                    "children": [unprotected_stock()],
                }
            ]
        )
    )

    events, _interval_rules, _index = analyze_route_protection(
        route,
        "r1",
        {"hydroxyl_tms": tms_rule()},
        config=ProtectionAnalysisConfig(collect_interval_rules=False),
    )

    assert len(events) == 1
    assert events[0].trace_status == "introduced"
    assert events[0].source_type == "introduced"
    assert events[0].protection_node_id == "r1"
    assert events[0].lifetime_steps == 0


def test_ambiguous_trace_is_reported():
    protected_child = protected_mol(in_stock=True)
    route = route_with_child(
        protected_mol(
            children=[
                {
                    "type": "reaction",
                    "metadata": {"smiles": f"{PROTECTED}>>{PROTECTED}"},
                    "children": [protected_child],
                },
                {
                    "type": "reaction",
                    "metadata": {"smiles": f"{PROTECTED}>>{PROTECTED}"},
                    "children": [protected_mol(in_stock=True)],
                },
            ]
        )
    )

    events, _interval_rules, _index = analyze_route_protection(
        route,
        "r2",
        {"hydroxyl_tms": tms_rule()},
        config=ProtectionAnalysisConfig(collect_interval_rules=False),
    )

    assert len(events) == 1
    assert events[0].trace_status == "ambiguous"
    assert events[0].failure_reason == "multiple_candidate_ancestors"


def test_normalizes_step_local_maps_before_tracing():
    route = route_with_child(
        protected_mol(
            children=[
                {
                    "type": "reaction",
                    "metadata": {"smiles": "[CH4:8]>>" + PROTECTED},
                    "children": [{"smiles": "C", "type": "mol", "in_stock": True}],
                }
            ]
        )
    )

    events, _interval_rules, _index = analyze_route_protection(
        route,
        "r3",
        {"hydroxyl_tms": tms_rule()},
        config=ProtectionAnalysisConfig(collect_interval_rules=False),
    )

    assert len(events) == 1
    assert events[0].trace_status == "introduced"
    assert events[0].source_type == "introduced"


def test_multicenter_deprotection_is_kept(monkeypatch):
    route = route_with_child(protected_mol(in_stock=True))
    monkeypatch.setattr(
        "route_inspector.protection.analysis.reaction_center_atoms",
        lambda _reaction_smiles: frozenset({3, 4, 8}),
    )

    events, _interval_rules, _index = analyze_route_protection(
        route,
        "r4",
        {"hydroxyl_tms": tms_rule()},
        config=ProtectionAnalysisConfig(
            collect_interval_rules=False,
            include_multicenter=True,
        ),
    )

    assert len(events) == 1
    assert events[0].multicenter_status == "deprotection_plus_other"
    assert events[0].n_other_reaction_centers == 1


def test_single_rule_summary_merges_by_querycgr():
    rows = summarize_single_rules(
        [
            ProtectionSingleRuleObservation(
                route_id="r1",
                pg_type="hydroxyl_tms",
                rule_smarts="[C:1]>>[O:1]",
                reaction_smiles="[CH4:1]>>[OH2:1]",
                event_id="r1:pg1",
            ),
            ProtectionSingleRuleObservation(
                route_id="r1",
                pg_type="hydroxyl_tms",
                rule_smarts="[C:1]>>[O:1]",
                reaction_smiles="[CH4:1]>>[OH2:1]",
                event_id="r1:pg1",
            ),
            ProtectionSingleRuleObservation(
                route_id="r2",
                pg_type="hydroxyl_tms",
                rule_smarts="[C:9]>>[O:9]",
                reaction_smiles="[CH4:9]>>[OH2:9]",
                event_id="r2:pg1",
            ),
        ]
    )

    assert rows == [
        {
            "source_pg_type": "hydroxyl_tms",
            "rule": "[C:1]>>[O:1]",
            "route_count": 2,
            "rule_count": 3,
            "route_ids": "r1,r2",
        }
    ]


def test_aggregate_single_rule_summary_merges_across_pg_types():
    rows = summarize_aggregate_single_rules(
        [
            ProtectionSingleRuleObservation(
                route_id="r1",
                pg_type="amine_boc",
                rule_smarts="[C:1]>>[O:1]",
                reaction_smiles="[CH4:1]>>[OH2:1]",
                event_id="r1:pg1",
            ),
            ProtectionSingleRuleObservation(
                route_id="r2",
                pg_type="carboxyl_methyl",
                rule_smarts="[C:9]>>[O:9]",
                reaction_smiles="[CH4:9]>>[OH2:9]",
                event_id="r2:pg1",
            ),
        ]
    )

    assert rows == [
        {
            "rule": "[C:1]>>[O:1]",
            "pg_types": "amine_boc,carboxyl_methyl",
            "route_count": 2,
            "rulec_count": 2,
            "route_ids": "r1,r2",
        }
    ]


def test_protection_free_routes_are_reported():
    route = {"smiles": "CCO", "type": "mol", "in_stock": False, "children": []}
    result = analyze_protection_in_routes(
        [route],
        protection_rules={"hydroxyl_tms": tms_rule()},
        config=ProtectionAnalysisConfig(collect_interval_rules=False),
    )

    assert result.summary["n_protection_free_routes"] == 1


def test_protection_outputs_write_single_rules_and_drop_network_edges(tmp_path):
    stale_network = tmp_path / "n1_protection_network_edges.tsv"
    stale_network.write_text("old\n", encoding="utf-8")
    stale_free_routes = tmp_path / "n1_protection_free_routes.tsv"
    stale_free_routes.write_text("old\n", encoding="utf-8")
    result = ProtectionAnalysisResult(
        route_stats_rows=[],
        event_rows=[],
        interval_rule_rows=[],
        single_rule_rows=[
            {
                "source_pg_type": "hydroxyl_tms",
                "rule": "[C:1]>>[O:1]",
                "route_count": 1,
                "rule_count": 2,
                "route_ids": "1",
            }
        ],
        aggregate_single_rule_rows=[
            {
                "rule": "[C:1]>>[O:1]",
                "pg_types": "hydroxyl_tms",
                "route_count": 1,
                "rulec_count": 2,
                "route_ids": "1",
            }
        ],
        group_summary_rows=[],
        rule_family_rows=[],
        trace_failure_rows=[],
        summary={},
    )

    output_info = write_protection_outputs(result, tmp_path, dataset_prefix="n1")

    assert not stale_network.exists()
    assert not stale_free_routes.exists()
    assert "network_edges" not in output_info["output_files"]
    assert "protection_free_routes" not in output_info["output_files"]
    assert (tmp_path / "n1_protection_single_rules.tsv").exists()
    assert (tmp_path / "n1_protection_agg_single_rule.tsv").exists()
