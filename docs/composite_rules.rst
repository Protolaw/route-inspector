Composite Rule Extraction
=========================

Principle
---------

The extractor reads PaRoutes-style ``mol``/``reaction`` route trees, normalizes
route atom maps, extracts individual SynPlanner reaction rules, and emits
ordered rule sequences.

A sequence is considered a composite rule only when adjacent route reactions
share the same reaction-center region on the route molecule. For a valid
ordered chain ``t1, t2, t3``, the extractor emits ``t1$t2``, ``t2$t3``, and
``t1$t2$t3``. It does not emit non-contiguous combinations such as ``t1$t3``.

Basic Command
-------------

.. code-block:: bash

   python -m route_analysis.cli extract-composite-rules \
     --routes-json data/n1-routes.json \
     --output comp_output/n1/n1.tsv \
     --config configs/rule_extraction_functional_groups.yaml \
     --ignore-errors

Large Runs
----------

For large collections, use worker processes and chunked dispatch:

.. code-block:: bash

   python -m route_analysis.cli extract-composite-rules \
     --routes-json data/all_routes.json \
     --output comp_output/all/all \
     --config configs/rule_extraction_functional_groups.yaml \
     --n_cpu 4 \
     --worker-chunksize 64 \
     --skip-routes-without-composites-output

``--skip-routes-without-composites-output`` keeps summary counts but skips the
large SVG-ready JSON sidecar of routes without composite rules.

When normalized mapped reactions repeat across many routes, add
``--unique-reactions-first``:

.. code-block:: bash

   python -m route_analysis.cli extract-composite-rules \
     --routes-json data/all_routes.json \
     --output comp_output/all/all \
     --config configs/rule_extraction_functional_groups.yaml \
     --n_cpu 4 \
     --worker-chunksize 64 \
     --unique-reactions-first \
     --skip-routes-without-composites-output

This mode scans routes once for unique normalized reaction SMILES, extracts
each unique reaction rule once, then scans routes again to compose rule
sequences from the precomputed cache. Check ``unique_reactions_seen`` versus
``reactions_seen`` in the summary to see whether this mode helped.

Output Files
------------

For output prefix ``n1``, extraction writes:

.. code-block:: text

   n1_t1_single_rules.tsv
   n1_t2_composite_rules.tsv
   n1_t3_composite_rules.tsv
   ...
   n1_composite_rule_extraction_summary.json
   n1_routes_without_composite_rules.json

The single-rule TSV schema is:

.. code-block:: text

   Rule	popularity	Reference	Target_molecules

The composite-rule TSV schema is:

.. code-block:: text

   Composite_rule	output_reactants_num	popularity	Reference	Target_molecules

``Composite_rule`` is a ``$``-separated sequence of SynPlanner rule SMARTS.
``output_reactants_num`` is the estimated final leaf count after applying the
linear rule sequence. ``Reference`` is a comma-separated list of route IDs.
``popularity`` is ``len(Reference)``. ``Target_molecules`` lists unique
molecules where the rule sequence starts.

Routes Without Composite Rules
------------------------------

Successfully processed routes that yield no composite rules are written as a
JSON lookup compatible with ``get_route_svg_from_json``:

.. code-block:: text

   {route_id: route_tree}

Each route root receives ``metadata.composite_rule_extraction`` with target
SMILES, route-level reaction counts, and a reason such as ``no_reactions``,
``fewer_than_min_length_extracted_reactions``, or
``no_reaction_center_sharing_sequence``.

Composite Unwrapping
--------------------

Apply a composite rule sequence to a target molecule:

.. code-block:: bash

   python -m route_analysis.cli unwrap-composite-rule \
     --smiles 'CCO' \
     --composite-rule-tsv comp_output/n1/n1_t2_composite_rules.tsv \
     --row 0 \
     --output-json /private/tmp/composite_unwrapped_route.json \
     --output-svg /private/tmp/composite_unwrapped_route.svg
