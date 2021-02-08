.. _Plot-Options:

Analysis options
=================

Data Analysis
~~~~~~~~~~~~~
The ``mspypeline`` was developed to consolidate the proteomics data analysis workflow while providing both simplicity
and flexibility. The :ref:`usage of the GUI <gui-quickstart>` provides the opportunity to explore and visualize data
straight forward, without interaction of the code. The GUI offers the majority of data analysis options currently
available in the ``mspypeline``. Additional plots and more advanced analysis may be performed by specifying optional
arguments in the :ref:`configuration <default-yaml>` file or by accessing the package as a
:ref:`python module <python-quickstart>`.


Plots
*****

The GUI allows for the following plots:

Normalization plots:

* Normalization overview: :meth:`~mspypeline.BasePlotter.plot_normalization_overview`
* Heatmap overview: :meth:`~mspypeline.BasePlotter.plot_heatmap_overview_all_normalizers`

Outlier detection / Comparison plots:

* Detection counts: :meth:`~mspypeline.BasePlotter.plot_detection_counts`
* Number of detected proteins :meth:`~mspypeline.BasePlotter.plot_detected_proteins_per_replicate`
* Venn diagrams: :meth:`~mspypeline.BasePlotter.plot_venn_results`
* Group diagrams: :meth:`~mspypeline.BasePlotter.plot_venn_groups`
* PCA overview: :meth:`~mspypeline.BasePlotter.plot_pca_overview`
* Intensity histogram: :meth:`~mspypeline.BasePlotter.plot_intensity_histograms`
* Relative std: :meth:`~mspypeline.BasePlotter.plot_relative_std`
* Scatter replicates: :meth:`~mspypeline.BasePlotter.plot_scatter_replicates`
* Experiment comparison: :meth:`~mspypeline.BasePlotter.plot_experiment_comparison`
* Rank: :meth:`~mspypeline.BasePlotter.plot_rank`

Statistical inference plots:

* Pathway analysis: :meth:`~mspypeline.BasePlotter.plot_pathway_analysis`
* Go analysis: :meth:`~mspypeline.BasePlotter.plot_go_analysis`
* Volcano plot (R): :meth:`~mspypeline.BasePlotter.plot_r_volcano`

Additionally via python:

* :meth:`~mspypeline.BasePlotter.plot_kde`
* :meth:`~mspypeline.BasePlotter.plot_boxplot`
* :meth:`~mspypeline.BasePlotter.plot_n_proteins_vs_quantile`
* :meth:`~mspypeline.BasePlotter.plot_intensity_heatmap`



.. _plotters:

Plotters
~~~~~~~~~
To perform data analysis and visualisation the Plotters from the MSPlots module are used. The MSPypeline currently
contains two Plotters, the ``mspypeline.BasePlotter`` and the ``mspypeline.MaxQuantPlotter``.

Base Plotter
*************
The BasePlotter provides all plots listed above. No quality control report is provided.

MaxQuant Plotter
*****************
The MaxQuant Plotter is a child class of the Base Plotter and inherits all functionality and plotting options listed
above. Additionally, the MaxQuant Plotter provides a quality control report based on supplementary MaxQuant data.

Quality Control Report
***********************
A quality control report for the output of MaxQuant. See :meth:`~mspypeline.MaxQuantPlotter.create_report`
for a description of the output.
This quality control report is specifically designed to process supplementary MaxQuant files
available to generate a multi-page pdf document. Here, the quality of the raw data can be investigated, for instance,
the influence of experimental parameters such as protein digestion, technical information like retention time and
length of individual samples as well as the corresponding chromatograms or the number of peptides measured and
sequenced and the proportion of contamination of protein intensities.
