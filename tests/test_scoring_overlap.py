import csv
import json

from alchems.composite_rules.extract import RouteProcessingStats
from alchems.scoring import overlap


def test_score_composite_rule_overlap_reports_order_sensitive_matches(
    tmp_path,
    monkeypatch,
):
    extracted = tmp_path / "extracted.tsv"
    reference_routes = tmp_path / "routes.json"
    output_dir = tmp_path / "scores"

    extracted.write_text(
        "Composite_rule\tpopularity\tReference\n"
        "a$b\t3\t1,2,3\n"
        "b$a\t2\t4,5\n",
        encoding="utf-8",
    )
    reference_routes.write_text("[]", encoding="utf-8")

    def fake_reference(*_args, **_kwargs):
        return (
            overlap.CompositeRuleSet(
                source=str(reference_routes),
                popularity_by_rule={"a$b": 1, "x$y": 1},
                references_by_rule={"a$b": {"10"}, "x$y": {"11"}},
                rows_seen=2,
            ),
            RouteProcessingStats(routes_seen=2, routes_with_composites=2),
            [],
        )

    monkeypatch.setattr(
        overlap,
        "reference_composite_rules_from_routes",
        fake_reference,
    )

    summary = overlap.score_composite_rule_overlap(
        [extracted],
        reference_routes,
        output_dir,
        rule_extractor=object(),
    )

    output_path = output_dir / "composite_rule_overlap_scores.tsv"
    with output_path.open() as file:
        rows = list(csv.DictReader(file, delimiter="\t"))

    assert json.loads((output_dir / "composite_rule_overlap_summary.json").read_text())
    assert rows[0]["overlap_unique_composite_rules"] == "1"
    assert rows[0]["extracted_overlap_ratio"] == "0.5"
    assert rows[0]["reference_coverage_ratio"] == "0.5"
    assert rows[0]["popularity_overlap_ratio"] == "0.6"

    with (output_dir / "composite_rule_overlap_matches.tsv").open() as file:
        match_rows = list(csv.DictReader(file, delimiter="\t"))

    assert match_rows[0]["Composite_rule"] == "a$b"
    assert match_rows[0]["present_in_reference"] == "True"
    assert match_rows[1]["Composite_rule"] == "b$a"
    assert match_rows[1]["present_in_reference"] == "False"
