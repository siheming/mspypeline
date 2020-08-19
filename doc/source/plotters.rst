.. _plotters:

Plotters
===============

Base Plotter
~~~~~~~~~~~~

Report
******
The Base Plotter does not provide a report.

Plots
*****

The GUI allows for these plots:

Normalization plots:

* Normalization overview: :meth:`~mspypeline.BasePlotter.plot_normalization_overview`
* Heatmap overview: :meth:`~mspypeline.BasePlotter.plot_intensity_heatmap`

Comparison plots:

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
* Pathway timecourse: :meth:`~mspypeline.BasePlotter.plot_pathway_timecourse`
* Go analysis: :meth:`~mspypeline.BasePlotter.plot_go_analysis`
* Volcano plot (R): :meth:`~mspypeline.BasePlotter.plot_r_volcano`

Additionally via python:

* :meth:`~mspypeline.BasePlotter.plot_kde`
* :meth:`~mspypeline.BasePlotter.plot_boxplot`
* :meth:`~mspypeline.BasePlotter.plot_n_proteins_vs_quantile`

Max Quant Plotter
~~~~~~~~~~~~~~~~~

Report
*******
A quality control report for the output of MaxQuant. See :meth:`~mspypeline.MaxQuantPlotter.create_report`
for a description of the output.

Plots
*****
The MaxQuant Plotter provides the same plots as the Base Plotter.
