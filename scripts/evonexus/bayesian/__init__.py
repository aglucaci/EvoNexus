"""Pipeline-side driver for the Bayesian modelling stage.

The model itself lives in the standalone ``evonexus_bayes`` distribution
(``packages/evonexus-bayes``). What remains here is the part that is specific
to this repository rather than to the statistics: resolving ``config/``, the
``RunRecord`` manifest, the sensitivity grid, and the layout of the output
directory that ``evonexus.reporting`` reads back.
"""
