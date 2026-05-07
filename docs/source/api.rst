Renoir API
==========

This page documents every public function in Renoir.
Internal helper functions (``entropy``, ``ISM``, ``PEM``, ``min_max``, etc.) are
intentionally omitted — they are implementation details not intended for direct use.

----

Core functions (``Renoir.renoir``)
----------------------------------

These functions cover the full pre-processing and score-computation pipeline,
from building the neighborhood graph through to producing the final
``AnnData`` of neighborhood scores.

.. autofunction:: Renoir.renoir.compute_neighborhood_scores

.. autofunction:: Renoir.renoir.neighborhood

.. autofunction:: Renoir.renoir.get_ligand_target

.. autofunction:: Renoir.renoir.compute_celltype_expression

.. autofunction:: Renoir.renoir.register_table

----

Downstream analysis (``Renoir.downstream``)
--------------------------------------------

These functions operate on the neighborhood-score ``AnnData`` produced by
:func:`Renoir.renoir.compute_neighborhood_scores` and cover the full
downstream workflow: clustering, domain identification, pathway activity,
ligand ranking, and visualization.

.. autofunction:: Renoir.downstream.get_msig

.. autofunction:: Renoir.downstream.create_cluster

.. autofunction:: Renoir.downstream.downstream_analysis

.. autofunction:: Renoir.downstream.spot_v_spot

.. autofunction:: Renoir.downstream.pcs_v_neighbscore

.. autofunction:: Renoir.downstream.ligand_ranking

.. autofunction:: Renoir.downstream.sankeyPlot
