import argparse
import json

from route_inspector.composite_rules import extract


def test_extract_run_writes_routes_without_composites_file(tmp_path, monkeypatch):
    routes_json = tmp_path / "routes.json"
    routes_json.write_text(
        json.dumps(
            [
                {
                    "smiles": "CCO",
                    "type": "mol",
                    "in_stock": False,
                    "children": [],
                }
            ]
        ),
        encoding="utf-8",
    )
    output = tmp_path / "n1.tsv"

    monkeypatch.setattr(
        extract.SynPlannerRuleExtractor,
        "from_args",
        classmethod(lambda cls, args: object()),
    )

    exit_code = extract.run(
        argparse.Namespace(
            routes_json=routes_json,
            output=output,
            config=None,
            environment_atom_count=1,
            include_rings=False,
            keep_leaving_groups=True,
            keep_incoming_groups=False,
            reactor_validation=False,
            limit=None,
            min_length=2,
            max_length=5,
            routes_without_composites_output=None,
            skip_routes_without_composites_output=False,
            ignore_errors=False,
            progress_interval=0,
            n_cpu=1,
            worker_chunksize=16,
            max_pending_chunks=None,
        )
    )

    assert exit_code == 0

    no_rules_path = tmp_path / "n1_routes_without_composite_rules.json"
    routes_without_rules = json.loads(no_rules_path.read_text())

    assert list(routes_without_rules) == ["0"]
    assert routes_without_rules["0"]["smiles"] == "CCO"
    assert routes_without_rules["0"]["metadata"]["composite_rule_extraction"] == {
        "route_id": 0,
        "target_smiles": "CCO",
        "reactions_seen": 0,
        "extracted_reaction_rules": 0,
        "skipped_reactions": 0,
        "reason": "no_reactions",
    }

    summary = json.loads(
        (tmp_path / "n1_composite_rule_extraction_summary.json").read_text()
    )
    assert summary["routes_without_composite_rules"] == 1
    assert summary["routes_without_composite_rules_file"] == str(no_rules_path)
    assert summary["routes_without_composite_rules_output_skipped"] is False
    assert summary["routes_without_composite_rules_by_reason"] == {"no_reactions": 1}


def test_extract_run_accepts_output_dir_and_infers_dataset_prefix(tmp_path, monkeypatch):
    routes_json = tmp_path / "data" / "clean" / "n1_routes.json"
    routes_json.parent.mkdir(parents=True)
    routes_json.write_text(
        json.dumps(
            [
                {
                    "smiles": "CCO",
                    "type": "mol",
                    "in_stock": False,
                    "children": [],
                }
            ]
        ),
        encoding="utf-8",
    )
    output_dir = tmp_path / "outputs" / "n1" / "10_composite_rules"

    monkeypatch.setattr(
        extract.SynPlannerRuleExtractor,
        "from_args",
        classmethod(lambda cls, args: object()),
    )

    exit_code = extract.run(
        argparse.Namespace(
            routes_json=routes_json,
            output=None,
            output_dir=output_dir,
            config=None,
            environment_atom_count=1,
            include_rings=False,
            keep_leaving_groups=True,
            keep_incoming_groups=False,
            reactor_validation=False,
            limit=None,
            min_length=2,
            max_length=5,
            routes_without_composites_output=None,
            skip_routes_without_composites_output=False,
            unique_reactions_first=False,
            ignore_errors=False,
            progress_interval=0,
            n_cpu=1,
            worker_chunksize=16,
            max_pending_chunks=None,
        )
    )

    assert exit_code == 0
    assert (output_dir / "n1_composite_rule_extraction_summary.json").exists()
    assert (output_dir / "n1_routes_without_composite_rules.json").exists()
    assert (output_dir / "summary.json").exists()
    assert (output_dir / "errors.tsv").exists()
    manifest = json.loads((output_dir / "manifest.json").read_text())
    assert manifest["command_name"] == "extract-composite-rules"
    assert manifest["input_files"] == [str(routes_json)]


def test_extract_run_can_skip_routes_without_composites_file(tmp_path, monkeypatch):
    routes_json = tmp_path / "routes.json"
    routes_json.write_text(
        json.dumps(
            [
                {
                    "smiles": "CCO",
                    "type": "mol",
                    "in_stock": False,
                    "children": [],
                }
            ]
        ),
        encoding="utf-8",
    )
    output = tmp_path / "n1.tsv"

    monkeypatch.setattr(
        extract.SynPlannerRuleExtractor,
        "from_args",
        classmethod(lambda cls, args: object()),
    )

    exit_code = extract.run(
        argparse.Namespace(
            routes_json=routes_json,
            output=output,
            config=None,
            environment_atom_count=1,
            include_rings=False,
            keep_leaving_groups=True,
            keep_incoming_groups=False,
            reactor_validation=False,
            limit=None,
            min_length=2,
            max_length=5,
            routes_without_composites_output=None,
            skip_routes_without_composites_output=True,
            ignore_errors=False,
            progress_interval=0,
            n_cpu=1,
            worker_chunksize=16,
            max_pending_chunks=None,
        )
    )

    assert exit_code == 0

    no_rules_path = tmp_path / "n1_routes_without_composite_rules.json"
    assert not no_rules_path.exists()

    summary = json.loads(
        (tmp_path / "n1_composite_rule_extraction_summary.json").read_text()
    )
    assert summary["routes_without_composite_rules"] == 1
    assert summary["routes_without_composite_rules_file"] is None
    assert summary["routes_without_composite_rules_output_skipped"] is True
    assert summary["routes_without_composite_rules_by_reason"] == {"no_reactions": 1}


def test_extract_run_unique_reactions_first_mode(tmp_path, monkeypatch):
    routes_json = tmp_path / "routes.json"
    routes_json.write_text(
        json.dumps(
            [
                {
                    "smiles": "CCO",
                    "type": "mol",
                    "in_stock": False,
                    "children": [],
                }
            ]
        ),
        encoding="utf-8",
    )
    output = tmp_path / "n1.tsv"

    monkeypatch.setattr(
        extract.SynPlannerRuleExtractor,
        "from_args",
        classmethod(lambda cls, args: object()),
    )

    exit_code = extract.run(
        argparse.Namespace(
            routes_json=routes_json,
            output=output,
            config=None,
            environment_atom_count=1,
            include_rings=False,
            keep_leaving_groups=True,
            keep_incoming_groups=False,
            reactor_validation=False,
            limit=None,
            min_length=2,
            max_length=5,
            routes_without_composites_output=None,
            skip_routes_without_composites_output=True,
            unique_reactions_first=True,
            ignore_errors=False,
            progress_interval=0,
            n_cpu=1,
            worker_chunksize=16,
            max_pending_chunks=None,
        )
    )

    assert exit_code == 0

    summary = json.loads(
        (tmp_path / "n1_composite_rule_extraction_summary.json").read_text()
    )
    assert summary["extraction_mode"] == "unique_reactions_first"
    assert summary["unique_reactions_seen"] == 0
    assert summary["routes_without_composite_rules"] == 1
