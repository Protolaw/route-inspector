import csv
import json

from route_analysis.io import (
    write_composite_routes_without_rules,
    write_composite_rules,
    write_composite_summary,
)


def test_write_composite_rules_splits_by_size_and_adds_popularity(tmp_path):
    output = tmp_path / "n1.tsv"
    summary = write_composite_rules(
        output,
        {
            ("z",): {4},
            ("a", "b"): {2, 1},
            ("b", "c"): {3},
            ("a", "b", "c"): {1},
        },
        target_molecules_by_sequence={
            ("z",): {4: {"CCCl"}},
            ("a", "b"): {2: {"CCO"}, 1: {"CCN"}},
            ("b", "c"): {3: {"CCC"}},
            ("a", "b", "c"): {1: {"CCN"}},
        },
    )

    t1 = tmp_path / "n1_t1_single_rules.tsv"
    t2 = tmp_path / "n1_t2_composite_rules.tsv"
    t3 = tmp_path / "n1_t3_composite_rules.tsv"
    assert set(summary["output_files"]) == {"1", "2", "3"}
    assert summary["unique_single_rules"] == 1
    assert summary["unique_composite_rules_by_size"] == {"2": 2, "3": 1}
    assert t1.exists()
    assert t2.exists()
    assert t3.exists()

    with t1.open() as file:
        rows = list(csv.DictReader(file, delimiter="\t"))

    assert rows[0] == {
        "Rule": "z",
        "popularity": "1",
        "Reference": "4",
        "Target_molecules": "CCCl",
    }

    with t2.open() as file:
        rows = list(csv.DictReader(file, delimiter="\t"))

    assert rows[0] == {
        "Composite_rule": "a$b",
        "output_reactants_num": "1",
        "popularity": "2",
        "Reference": "1,2",
        "Target_molecules": "CCN,CCO",
    }
    assert rows[1] == {
        "Composite_rule": "b$c",
        "output_reactants_num": "1",
        "popularity": "1",
        "Reference": "3",
        "Target_molecules": "CCC",
    }


def test_summary_file_name(tmp_path):
    output = tmp_path / "n1.tsv"
    path = write_composite_summary(output, {"ok": True})

    assert path == tmp_path / "n1_composite_rule_extraction_summary.json"
    assert json.loads(path.read_text()) == {"ok": True}


def test_write_composite_routes_without_rules(tmp_path):
    output = tmp_path / "n1.tsv"
    path = write_composite_routes_without_rules(
        output,
        {
            7: {
                "target_smiles": "CCO",
            }
        },
    )

    assert path == tmp_path / "n1_routes_without_composite_rules.json"
    assert json.loads(path.read_text()) == {
        "7": {
            "target_smiles": "CCO",
        }
    }
