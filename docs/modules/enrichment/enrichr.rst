enrichR module
--------------

.. autoclass:: magine.enrichment.enrichr.Enrichr
   :members:
   :undoc-members:
   :show-inheritance:


Using a ExperimentalData instance, we can run enrichR for all the databases using a simple wrapper around.

.. autofunction:: magine.enrichment.enrichr.run_enrichment_for_project


We also provide some tools to clean up and standardize enrichRs output.


Functions to cleanup enrichR term names
+++++++++++++++++++++++++++++++++++++++
*Note: this are in progress and not fully tested! Warning!*

.. autofunction:: magine.enrichment.enrichr.clean_term_names

.. autofunction:: magine.enrichment.enrichr.clean_lincs

.. autofunction:: magine.enrichment.enrichr.clean_drug_pert_geo

.. autofunction:: magine.enrichment.enrichr.clean_tf_names

Download reference databases
++++++++++++++++++++++++++++
.. autofunction:: magine.enrichment.enrichr.get_background_list
