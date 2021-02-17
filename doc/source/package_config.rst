Custom Protein Lists
=====================

Pathway and GO analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

| Pathway and GO analysis protein files are lists of proteins that are somewhat associated with the respective pathway
  or GO term. Importantly, these lists have been manually created and should merely fulfill the purpose of giving an
  idea or providing an example for a potential analysis (see :ref:`gallery`).

  * GO protein lists are required for the GO analysis plot using :meth:`~mspypeline.BasePlotter.plot_go_analysis`, in
    which an enrichment analysis for each selected GO Term file is created.
  * Pathway protein lists serve multiple plotting options including the rank plot using
    :meth:`~mspypeline.BasePlotter.plot_rank`, the pathway analysis plot using
    :meth:`~mspypeline.BasePlotter.plot_pathway_analysis` and the volcano plot using
    :meth:`~mspypeline.BasePlotter.plot_r_volcano`.

| Pathway (and GO term) lists are configured globally for a data analysis and not individually for each plot. If one or
  more pathway protein lists are selected in the GUI or :ref:`configs file <default-yaml>`, these pathways will
  be used for any of the three listed plots if they are being created.
| To change the choice of pathways or GO terms for an analysis, the previously selected pathways have to be un-checked
  in the corresponding selection box in the GUI or the *"pathways:"* or *"go_term:"* arguments in the configs file have
  to be edited manually.

.. tip::
    | Any desired pathway and GO analysis protein file can be manually provided to the ``mspypeline`` by the user. The
      file simply has to:
    | 1. follow the *one-column-txt-format* that can be seen in the exemplary files listed below,
    | 2. be saved in the *.../mspypeline/config/pathway or go_term* directory, where all the other files are stored.

.. attention::
    All protein and GO analysis protein files provided here are based on the **HUMAN** proteome.

.. _pathway-proteins:

Pathways
~~~~~~~~
.. toctree::
   :glob:

   pathways/*

.. _go-term-proteins:

GO Terms
~~~~~~~~
.. toctree::
   :glob:

   go_terms/*