Protection Analysis
===================

Concept
-------

Protection analysis detects curated chython deprotection events in mapped route
steps, traces the protected atom or group backward through the normalized route,
and reports whether the protecting group was introduced during synthesis or
came from an in-stock precursor.

Command
-------

.. code-block:: bash

   python -m route_analysis.cli analyze-protection \
     --routes-json data/n1-routes.json \
     --composite-rule-tsv comp_output/n1 \
     --output-dir protection_out/n1 \
     --config configs/protection_analysis.yaml \
     --include-multicenter \
     --deprotection-first \
     --querycgr-compare \
     --ignore-errors

Use ``--n_cpu`` for route-level parallelism. If running through ``conda run``,
add ``--no-capture-output`` to see progress live.

Outputs
-------

The analysis writes route statistics, protection events, interval
composite-rule rows, group summaries, rule-family summaries, network edges,
trace failures, and a JSON summary.

Important event-level fields include:

- protecting-group type;
- source type: introduced, stock, or unresolved;
- deprotection route depth;
- source route depth;
- lifetime steps between source and deprotection;
- multicenter/deprotective-combo status;
- selected interval composite rules.

Tutorial
--------

See ``tutorials/analyze_protection.ipynb`` for a visual route-level walkthrough
using ``get_route_svg_from_json`` and aggregate protection summary tables.
