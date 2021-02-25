Custom Gene Lists
=====================

Pathway and GO analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

| Pathway and GO analysis gene set files are lists of genes that are somewhat associated with the respective pathway
  or GO term. Upon data preparation by the :class:`MQReader` measured proteins are indexed by their gene name
  originating from the FASTA file (header) to which the protein was mapped. Thus, measured protein intensities can be
  analysed using functional gene sets like those that are incorporated by the ``mspypeline``.
| Since some analysis methods deploy such pathway and GO term lists the provided files can be used to get an idea and
  generate an example for a potential analysis (see :ref:`gallery`).

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
    | 2. store the file on one of these two locations:
         - saved in a *pathways* and *go_terms* dir where the experiment data is stored
         - saved in the *.../mspypeline/config/pathway or go_term* directory, where all the other files are stored.
           (files saved here are available for all experiments)

.. attention::
    * All pathway and GO analysis protein files provided here are based on the **HUMAN** proteome.
    * All pathway and GO analysis protein files are retrieved from the open source
      `GSEA Molecular Signature Data Base <https://www.gsea-msigdb.org/gsea/msigdb/index.jsp>`__ (22. Feb. 2021)

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