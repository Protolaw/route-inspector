from __future__ import annotations

import json
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]

if __package__ in (None, "") and str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from route_analysis.alchemical_rules.alchemical import rule_cgr_key
from route_analysis.io import (
    default_classification_output_path,
    default_classification_summary_path,
    normalize_n_cpu,
    read_tsv_rows,
    resolve_classification_output_paths,
    resolve_existing_path,
    rule_column,
    setup_runtime_cache_dirs,
    write_json,
    write_tsv,
)


def _rule_cgr_worker(item: tuple[int, str]) -> tuple[int, str, str, str]:
    index, smarts = item
    try:
        return index, smarts, rule_cgr_key(smarts), ""
    except Exception as exc:
        return index, smarts, "", str(exc)


def load_default_rule_cgrs(
    rules_tsv: Path,
    *,
    n_cpu: int = 1,
) -> tuple[dict[str, list[tuple[int, str]]], int, int]:
    default_rules: dict[str, list[tuple[int, str]]] = defaultdict(list)
    parsed = 0
    errors = 0
    fieldnames, rows = read_tsv_rows(rules_tsv)
    column = rule_column(fieldnames, ("rule_smarts", "Rule", "SMARTS"))
    work_items = [
        (index, row.get(column, "").strip())
        for index, row in enumerate(rows)
        if row.get(column, "").strip()
    ]
    n_cpu = normalize_n_cpu(n_cpu)
    if n_cpu > 1 and work_items:
        with ProcessPoolExecutor(max_workers=n_cpu) as executor:
            results = executor.map(_rule_cgr_worker, work_items)
            for index, smarts, cgr_key, error in results:
                if error:
                    errors += 1
                    continue
                default_rules[cgr_key].append((index, smarts))
                parsed += 1
        return default_rules, parsed, errors

    for index, row in enumerate(rows):
        smarts = row.get(column, "").strip()
        if not smarts:
            continue
        try:
            default_rules[rule_cgr_key(smarts)].append((index, smarts))
            parsed += 1
        except Exception:
            errors += 1
    return default_rules, parsed, errors


def classify_alchemical_rules(
    alchemical_rules_tsv: Path,
    default_rules_tsv: Path,
    output: Path,
    summary_path: Path,
    *,
    n_cpu: int = 1,
) -> dict[str, Any]:
    alchemical_rules_tsv = resolve_existing_path(alchemical_rules_tsv)
    default_rules_tsv = resolve_existing_path(default_rules_tsv)
    output, summary_path = resolve_classification_output_paths(
        alchemical_rules_tsv,
        output,
        summary_path,
    )
    default_cgrs, default_rules_parsed, default_rule_errors = load_default_rule_cgrs(
        default_rules_tsv,
        n_cpu=n_cpu,
    )
    positive = 0
    negative = 0
    rows_seen = 0
    cgr_errors = 0

    input_fieldnames, input_rows = read_tsv_rows(alchemical_rules_tsv)
    alchemical_column = rule_column(
        input_fieldnames,
        ("Alchemical_rule", "Alchemical_rules", "rule_smarts"),
    )
    fieldnames = list(input_fieldnames)
    for column in (
        "classification",
        "Matched_default_rule_ids",
        "Matched_default_rules",
    ):
        if column not in fieldnames:
            fieldnames.append(column)

    output_rows = []
    positive_rows = []
    negative_rows = []
    n_cpu = normalize_n_cpu(n_cpu)
    cgr_keys_by_row: dict[int, tuple[str, bool]] = {}
    work_items = [
        (row_index, row.get(alchemical_column, ""))
        for row_index, row in enumerate(input_rows)
    ]
    if n_cpu > 1 and work_items:
        with ProcessPoolExecutor(max_workers=n_cpu) as executor:
            for row_index, _rule, cgr_key, error in executor.map(
                _rule_cgr_worker,
                work_items,
            ):
                cgr_keys_by_row[row_index] = (cgr_key, bool(error))

    for row_index, row in enumerate(input_rows):
        rows_seen += 1
        if n_cpu > 1:
            cgr_key, had_error = cgr_keys_by_row.get(row_index, ("", True))
            if had_error:
                cgr_errors += 1
                cgr_key = row.get("Alchemical_cgr", "").strip()
        else:
            try:
                cgr_key = rule_cgr_key(row[alchemical_column])
            except Exception:
                cgr_errors += 1
                cgr_key = row.get("Alchemical_cgr", "").strip()

        matches = default_cgrs.get(cgr_key, []) if cgr_key else []
        if matches:
            classification = "negative"
            negative += 1
        else:
            classification = "positive"
            positive += 1

        row["classification"] = classification
        row["Matched_default_rule_ids"] = ",".join(
            str(index) for index, _smarts in matches
        )
        row["Matched_default_rules"] = " || ".join(
            smarts for _index, smarts in matches
        )
        output_rows.append(row)
        if classification == "positive":
            positive_rows.append(row.copy())
        else:
            negative_rows.append(row.copy())

    write_tsv(output, fieldnames, output_rows)
    positive_output = output.with_name(f"{output.stem}_pos{output.suffix or '.tsv'}")
    negative_output = output.with_name(f"{output.stem}_neg{output.suffix or '.tsv'}")
    write_tsv(positive_output, fieldnames, positive_rows)
    write_tsv(negative_output, fieldnames, negative_rows)

    summary = {
        "alchemical_rules_tsv": str(alchemical_rules_tsv),
        "default_rules_tsv": str(default_rules_tsv),
        "output": str(output),
        "positive_output": str(positive_output),
        "negative_output": str(negative_output),
        "rows_seen": rows_seen,
        "positive": positive,
        "negative": negative,
        "default_rules_parsed": default_rules_parsed,
        "default_rule_parse_errors": default_rule_errors,
        "alchemical_cgr_errors": cgr_errors,
        "n_cpu": n_cpu,
    }
    write_json(summary_path, summary)
    summary["summary_file"] = str(summary_path)
    write_json(summary_path, summary)
    return summary


def run(args: argparse.Namespace) -> int:
    setup_runtime_cache_dirs()
    output = args.output or default_classification_output_path(args.alchemical_rules_tsv)
    summary_path = args.summary or (
        args.output if args.output else default_classification_summary_path(output)
    )
    summary = classify_alchemical_rules(
        resolve_existing_path(args.alchemical_rules_tsv),
        resolve_existing_path(args.default_rules_tsv),
        output,
        summary_path,
        n_cpu=args.n_cpu,
    )
    print(json.dumps(summary, indent=2), flush=True)
    return 0
