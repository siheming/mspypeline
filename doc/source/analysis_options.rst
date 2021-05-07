.. _Plot-Options:

Analysis options
=================

Data Analysis
~~~~~~~~~~~~~
| The ``mspypeline`` package was developed to consolidate the proteomics data analysis workflow while providing both
  simplicity and flexibility.
| The :ref:`usage of the GUI <gui-quickstart>` provides the opportunity to explore and visualize data in a straightforward manner,
  without demanding any interaction with the code. The GUI offers the majority of plotting options currently available
  in ``mspypeline``.
| Additional plots and more advanced analysis may be performed by specifying optional arguments in the
  :ref:`configuration <default-yaml>` file or by accessing the package as a :ref:`python module <python-quickstart>`.


Plots
*****

The GUI allows to create the following plots:

.. _norm-plots:

**Normalization plots:**

* Normalization overview: :meth:`~mspypeline.BasePlotter.plot_normalization_overview` (exemplary plot in :ref:`gallery <norm-overview>`)
* Heatmap overview: :meth:`~mspypeline.BasePlotter.plot_heatmap_overview_all_normalizers` (exemplary plot in :ref:`gallery <heatmap-overview>`)

.. _detection-plots:

**Outlier detection and comparison plots:**

* Detection counts: :meth:`~mspypeline.BasePlotter.plot_detection_counts` (exemplary plot in :ref:`gallery <detection-counts>`)
* Number of detected proteins :meth:`~mspypeline.BasePlotter.plot_detected_proteins_per_replicate` (exemplary plot in :ref:`gallery <detected-proteins>`)
* Venn diagrams: :meth:`~mspypeline.BasePlotter.plot_venn_results` (exemplary plot in :ref:`gallery <venn-rep>`)
* Group diagrams: :meth:`~mspypeline.BasePlotter.plot_venn_groups` (exemplary plot in :ref:`gallery <venn-group>`)
* PCA overview: :meth:`~mspypeline.BasePlotter.plot_pca_overview` (exemplary plot in :ref:`gallery <pca>`)
* Intensity histogram: :meth:`~mspypeline.BasePlotter.plot_intensity_histograms` (exemplary plot in :ref:`gallery <int-hist>`)
* Relative std: :meth:`~mspypeline.BasePlotter.plot_relative_std` (exemplary plot in :ref:`gallery <rel-std>`)
* Scatter replicates: :meth:`~mspypeline.BasePlotter.plot_scatter_replicates` (exemplary plot in :ref:`gallery <scatter-rep>`)
* Experiment comparison: :meth:`~mspypeline.BasePlotter.plot_experiment_comparison` (exemplary plot in :ref:`gallery <scatter-group>`)
* Rank: :meth:`~mspypeline.BasePlotter.plot_rank` (exemplary plot in :ref:`gallery <rank>`)

.. _statistic-plots:

**Statistical inference plots:**

* Pathway analysis: :meth:`~mspypeline.BasePlotter.plot_pathway_analysis` (exemplary plot in :ref:`gallery <pathway-analysis>`)
* Go analysis: :meth:`~mspypeline.BasePlotter.plot_go_analysis` (exemplary plot in :ref:`gallery <go-analysis>`)
* Volcano plot (R): :meth:`~mspypeline.BasePlotter.plot_r_volcano` (exemplary plot in :ref:`gallery <volcano>`)

.. _add-python-plots:

**Additionally via python:**

* :meth:`~mspypeline.BasePlotter.plot_kde` (exemplary plot in :ref:`gallery <kde>`)
* :meth:`~mspypeline.BasePlotter.plot_boxplot` (exemplary plot in :ref:`gallery <boxplot>`)
* :meth:`~mspypeline.BasePlotter.plot_n_proteins_vs_quantile` (exemplary plot in :ref:`gallery <proteins-vs-quantiles>`)
* :meth:`~mspypeline.BasePlotter.plot_intensity_heatmap` (exemplary plot in :ref:`gallery <int-heatmap>`)



.. _plotters:

Plotters
~~~~~~~~~
To perform data analysis and visualisation the Plotters from the MSPlots module are used. The MSPypeline currently
contains two Plotters, the ``mspypeline.BasePlotter`` and the ``mspypeline.MaxQuantPlotter``.

Base Plotter
*************
The Base Plotter provides all plots listed above. No quality control report is provided.

MaxQuant Plotter
*****************
The MaxQuant Plotter is a child class of the Base Plotter and inherits all functionality and plotting options listed
above. Additionally, the MaxQuant Plotter provides a quality control report based on supplementary MaxQuant data.

MaxQuant Report
***********************
| A quality control report for combining information from output tables of MaxQuant.
  See :meth:`~mspypeline.MaxQuantPlotter.create_report` for a description of the output and the
  :ref:`gallery <mqreport>` for an example of such a report.
| This quality control report is specifically designed to process available supplementary MaxQuant output tables
  and generate a multi-page pdf document. Here, the quality of the samples measured by mass spectrometry can be assessed,
  for instance the shape of a chromatogram or retention time and retention length of peptides. Additionally, information
  on sample quality such as missed cleavages, number of peptides measured and sequenced, and the proportion
  of contaminants among protein intensities is provided.
