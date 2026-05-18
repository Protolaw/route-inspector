Tutorials
=========

The ``tutorials/`` directory contains exploratory notebooks and small visual
checks. Recommended order:

1. ``test_composite.ipynb``: inspect extracted composite rules and route-level
   examples.
2. ``routes_without_composite.ipynb``: inspect routes that did or did not yield
   composite rules, including stacked histograms.
3. ``test_alchemical.ipynb``: inspect alchemical rules, pseudo-reactions, and
   unwrapping.
4. ``classifiy_alchemy.ipynb``: inspect positive/negative alchemical
   classification outputs.
5. ``test_scoring.ipynb``: understand overlap scoring with concrete route
   examples and classified alchemical rules.
6. ``test_diversity.ipynb``: compute rule fingerprints and pairwise
   similarities.
7. ``alchem_overlap_venn.ipynb``: compare alchemical rule sets and inspect
   unique rules between route collections.
8. ``analyze_protection.ipynb`` and ``test_protecting.ipynb``: inspect
   protection/deprotection event tracing and summary tables.
9. ``check_mol_rule.ipynb``: ad hoc rule and molecule inspection.

``tutorials/README.md`` keeps a short command-oriented walkthrough for running
the main CLIs. The notebooks are intentionally more detailed than the package
README and may include local paths from the development environment.
