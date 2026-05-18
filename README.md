# Route Inspector

Route Inspector is a standalone SynPlanner/chython toolkit for route-level rule
analysis. It works with PaRoutes-style `mol`/`reaction` route trees and uses
SynPlanner rule extraction instead of RDKit/AiZynthFinder templates.

The package can:

- extract one-step rules and ordered composite rules from synthesis routes;
- keep only contiguous rule sequences whose adjacent reaction centers are shared
  or chemically connected by the route molecule;
- unwrap composite or alchemical rules back into route JSON/SVG depictions;
- collect alchemical rules by collapsing unwrapped routes into pseudo-reactions;
- classify alchemical rules against default one-step rules by QueryCGR identity;
- score extracted composite-rule overlap against reference route JSON files;
- analyze protection/deprotection strategies in mapped routes.

## Installation

Use an environment where `synplan` and `chython-synplan` are
available. From the repository root:

```bash
conda activate synplan
python -m pip install -e .
route-inspector --help
```

For development without installation, run the module directly:

```bash
python -m route_analysis.cli --help
```

## Quick Start

Extract composite rules from a PaRoutes-style route JSON:

```bash
python -m route_analysis.cli extract-composite-rules \
  --routes-json data/n1-routes.json \
  --output comp_output/n1/n1.tsv \
  --config configs/rule_extraction_functional_groups.yaml \
  --ignore-errors
```

The command writes `*_t1_single_rules.tsv`, `*_t2_composite_rules.tsv`, larger
composite-rule TSVs, and a JSON extraction summary.

## Documentation

- `docs/composite_rules.rst`: extraction principles, output schema, large-run
  options, and composite unwrapping.
- `docs/alchemical_rules.rst`: alchemical extraction, pseudo-reaction output,
  classification, and alchemical unwrapping.
- `docs/protection_analysis.rst`: protection/deprotection tracing outputs and
  options.
- `docs/scoring.rst`: Retro-BLEU-style overlap scoring and positive/negative
  classification-aware metrics.
- `docs/tutorials.rst`: notebook map and recommended tutorial order.

Interactive examples live in `tutorials/`.
