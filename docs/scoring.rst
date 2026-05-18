Scoring
========

Composite-rule overlap scoring compares extracted composite-rule TSV files with
composite rules extracted from reference route JSON files.

Command
-------

.. code-block:: bash

   python -m route_analysis.cli score-composite-overlap \
     --extracted-tsv comp_output/n1 \
     --reference-routes-json data/n1-routes.json \
     --output /private/tmp/composite_rule_scores \
     --config configs/rule_extraction_functional_groups.yaml \
     --ignore-errors

The scoring command extracts reference rules on demand, compares
``$``-separated ordered rule sequences, and reports unique-rule overlap,
reference coverage, Jaccard overlap, and popularity-weighted overlap.

Classification-Aware Overlap
----------------------------

If classified alchemical rules are provided, the scorer also reports positive
and negative overlap components:

.. code-block:: bash

   python -m route_analysis.cli score-composite-overlap \
     --extracted-tsv comp_output/n1 \
     --reference-routes-json data/n1-routes.json \
     --classification-tsv reference/n1_classified_alchemical_rules_pos.tsv reference/n1_classified_alchemical_rules_neg.tsv \
     --output /private/tmp/composite_rule_scores \
     --config configs/rule_extraction_functional_groups.yaml \
     --ignore-errors

``pos_overlap`` rewards overlap linked to positive alchemical rules. A positive
alchemical rule is one whose QueryCGR is not equivalent to a default one-step
rule. ``neg_overlap`` is the corresponding penalty component for negative
alchemical rules. If one composite rule is linked to both positive and negative
alchemical rows, its contribution is split by the positive and negative weights
in the classification TSVs.

Relation to Retro-BLEU
----------------------

The scoring follows Retro-BLEU's high-level idea: extract ordered route
fragments, convert reactions into reusable rule identifiers, and score how many
extracted fragments occur in a reference vocabulary. This project differs by
using SynPlanner/chython rules, requiring reaction-center-sharing fragments, and
reading references directly from route JSON instead of a precomputed pickle.

Tutorial
--------

See ``tutorials/test_scoring.ipynb`` for a five-example explanation with route
depictions, overlap reactions, and positive/negative alchemical classifications.
