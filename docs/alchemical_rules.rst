Alchemical Rules
================

Concept
-------

Alchemical rules are derived from composite-rule applications. The unwrapped
route leaves that are in stock are collected as pseudo-reactants, the original
target molecule is used as pseudo-product, and SynPlanner rule extraction is
run on the resulting pseudo-reaction.

Extraction
----------

.. code-block:: bash

   python -m route_analysis.cli extract-alchemical-rules \
     --composite-rule-tsv comp_output/n1 \
     --output res_alchem/n1 \
     --config configs/rule_extraction_functional_groups.yaml \
     --ignore-errors

``--composite-rule-tsv`` accepts individual TSV files or directories containing
``*_composite_rules.tsv`` files. When ``--output`` is a directory, filenames are
derived from the composite-rule prefix.

Main outputs:

.. code-block:: text

   n1_alchemical_rules.tsv
   n1_alchemical_reactions.smi
   n1_alchemical_rule_collection_summary.json
   n1_alchemical_rule_collection_errors.tsv

The alchemical TSV links each alchemical rule back to all composite rules and
route IDs that generated it:

.. code-block:: text

   Alchemical_rule	output_reactants_num	popularity	route_ids_size	Reference	Target_molecules	composite_rules_size	Composite_rule_sizes	Composite_rules	Source_composite_rows	pseudo_reactions_size	Pseudo_reaction_ids

Duplicates are merged by comparing composed QueryCGR objects, not raw SMARTS
strings, so equivalent rules with different atom-map renderings are grouped.

Unwrapping
----------

Apply one alchemical rule directly to a target:

.. code-block:: bash

   python -m route_analysis.cli unwrap-alchemical-rule \
     --smiles 'CCO' \
     --alchemical-rule-tsv res_alchem/n1/n1_alchemical_rules.tsv \
     --row 0 \
     --output-json /private/tmp/alchemical_unwrapped_route.json \
     --output-svg /private/tmp/alchemical_unwrapped_route.svg

Classification
--------------

Classify alchemical rules against a default SynPlanner one-step rule TSV:

.. code-block:: bash

   python -m route_analysis.cli classify-alchemical-rules \
     --alchemical-rules-tsv res_alchem/n1/n1_alchemical_rules.tsv \
     --default-rules-tsv reference/pop_3_rules.tsv \
     --output reference/n1_classified_alchemical_rules.tsv

A rule is negative when its QueryCGR matches any default one-step rule.
Otherwise it is positive. The classifier also writes positive and negative
split TSVs.
