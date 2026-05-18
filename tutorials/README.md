# Tutorials

Run examples from the Route Inspector repository root after installing the
package in editable mode, or use `python -m route_analysis.cli` directly from
the root. SynPlanner is expected to be importable from the active `synplan`
environment.

## 1. Extract Composite Rules

```bash
python -m route_analysis.cli extract-composite-rules \
  --routes-json data/n1-routes.json \
  --output comp_output/n1/n1.tsv \
  --config configs/rule_extraction_functional_groups.yaml \
  --ignore-errors
```

Expected outputs:

```text
comp_output/n1/n1_t1_single_rules.tsv
comp_output/n1/n1_t2_composite_rules.tsv
comp_output/n1/n1_t3_composite_rules.tsv
comp_output/n1/n1_t4_composite_rules.tsv
comp_output/n1/n1_t5_composite_rules.tsv
comp_output/n1/n1_composite_rule_extraction_summary.json
```

## 2. Unwrap One Composite Rule

```bash
python -m route_analysis.cli unwrap-composite-rule \
  --smiles 'CCO' \
  --composite-rule-tsv comp_output/n1/n1_t2_composite_rules.tsv \
  --row 0 \
  --output-json /private/tmp/composite_unwrapped_route.json \
  --output-svg /private/tmp/composite_unwrapped_route.svg
```

## 3. Collect Alchemical Rules

```bash
python -m route_analysis.cli extract-alchemical-rules \
  --composite-rule-tsv comp_output/n1 \
  --output res_alchem/n1 \
  --config configs/rule_extraction_functional_groups.yaml \
  --ignore-errors
```

When `--output` is a directory, the collector derives all sidecar filenames:

```text
res_alchem/n1/n1_alchemical_rules.tsv
res_alchem/n1/n1_alchemical_reactions.smi
res_alchem/n1/n1_alchemical_rule_collection_summary.json
```

## 4. Classify Positive/Negative Rules

```bash
python -m route_analysis.cli classify-alchemical-rules \
  --alchemical-rules-tsv res_alchem/n1/n1_alchemical_rules.tsv \
  --default-rules-tsv reference/pop_3_rules.tsv \
  --output reference/n1_classified_alchemical_rules.tsv
```

## 5. Analyze Protection Strategies

```bash
python -m route_analysis.cli analyze-protection \
  --routes-json data/n1-routes.json \
  --composite-rule-tsv comp_output/n1 \
  --output-dir protection_out/n1 \
  --config configs/protection_analysis.yaml \
  --include-multicenter \
  --deprotection-first \
  --querycgr-compare \
  --ignore-errors
```

See `analyze_protection.ipynb` for a visual single-route example and the
pool-level summary tables.

## 6. Score Overlap

```bash
python -m route_analysis.cli score-composite-overlap \
  --extracted-tsv comp_output/n1 \
  --reference-routes-json data/n1-routes.json \
  --classification-tsv reference/n1_classified_alchemical_rules_pos.tsv reference/n1_classified_alchemical_rules_neg.tsv \
  --output /private/tmp/composite_rule_overlap_scores \
  --config configs/rule_extraction_functional_groups.yaml \
  --ignore-errors
```

The scoring output reports unique-rule overlap, reference coverage, Jaccard,
popularity-weighted extracted overlap, and optional classification-aware
`pos_overlap`/`neg_overlap` reward and penalty ratios.
