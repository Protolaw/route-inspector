from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]

if __package__ in (None, "") and str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from route_analysis.composite_rules.unwrap import unwrap_rule_sequence
from route_analysis.io import (
    read_alchemical_rule_from_tsv,
    resolve_existing_path,
    setup_runtime_cache_dirs,
    write_json,
)


def unwrap_alchemical_rule(
    target_smiles: str,
    alchemical_rule: str,
    *,
    route_id: int = 0,
    mark_leaves_in_stock: bool = True,
) -> dict[int, dict[str, Any]]:
    return unwrap_rule_sequence(
        target_smiles,
        [alchemical_rule],
        route_id=route_id,
        rule_key_prefix="alchemical",
        mark_leaves_in_stock=mark_leaves_in_stock,
    ).routes_json


def run(args: argparse.Namespace) -> int:
    setup_runtime_cache_dirs()

    alchemical_rule = args.alchemical_rule
    if alchemical_rule is None:
        alchemical_rule = read_alchemical_rule_from_tsv(
            resolve_existing_path(args.alchemical_rule_tsv),
            args.row,
        )

    routes_json = unwrap_alchemical_rule(
        args.smiles,
        alchemical_rule,
        route_id=args.route_id,
        mark_leaves_in_stock=not args.do_not_mark_leaves_in_stock,
    )

    if args.output_json:
        write_json(args.output_json, routes_json)
    else:
        print(json.dumps(routes_json, indent=2))

    if args.output_svg:
        from synplan.utils.visualisation import get_route_svg_from_json

        svg = get_route_svg_from_json(routes_json, args.route_id, labeled=args.labeled)
        args.output_svg.parent.mkdir(parents=True, exist_ok=True)
        args.output_svg.write_text(svg, encoding="utf-8")

    return 0

