.. currentmodule:: mspypeline

.. _gallery:

Gallery
========

In the following, each plot that can be created with the ``mspypeline`` will be shown and explained giving a minimal
code example to create this plot when using the package as a :ref:`python module <python-quickstart>`.

Plotter creation
^^^^^^^^^^^^^^^^^
Firstly, a plotter object has to be created too make the plots. Here, the :class:`~MaxQuantPlotter` is build from
the :class:`~MSPInitializer` class which creates and reads in the :ref:`configuration file <settings>` and initiates the
:class:`~MQReader` that loads the exemplary data set provided.

.. ipython:: python

    import pandas as pd
    import os
    from mspypeline import load_example_dataset, MaxQuantPlotter
    # load the data that is provided in a submodule
    init = load_example_dataset(configs={
        "pathways": ["BIOCARTA_EGF_PATHWAY.txt", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY.txt"],
        "go_terms": ["GO_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION.txt", "GO_INFLAMMATORY_RESPONSE.txt"]
        })
    plotter = MaxQuantPlotter.from_MSPInitializer(init)

    # create a second plotter without collapsed technical replicates
    init = load_example_dataset(configs={"has_techrep": False, "pathways":[]})
    plotter_with_tech_reps = MaxQuantPlotter.from_MSPInitializer(init)

define some helper functions

.. ipython:: python

    def select_fig(plts, idx):
        # sets active figure which can then be saved by the @savefig decorator
        plt.figure(plts[idx][0].number)


MaxQuant Report
^^^^^^^^^^^^^^^^^^^^^^^
Created using: :meth:`~mspypeline.MaxQuantPlotter.create_report`

    The MaxQuant report was built with the intention to offer a broad insight into the different sources of information
    from a MaxQuant output. Besides the protein intensities (from the *proteinGroups.txt* file) which are the only
    source of data for all other parts of the analysis with the :ref:`MaxQuant Plotter <plotters>`, further information
    about experimental and technical parameters of the experiment are taken into account. The MaxQuant report can
    function as quality control of the data and will output a multi-page pdf document composed of a variety of
    information and graphics.
    Make sure that :ref:`all MaxQuant files <file-readers>` are provided, which are used to create the report.

.. ipython:: python

    #plotter.create_report("./source/_static");
    print("skipping report")

The resulting `MaxQuant Report <./_static/MaxQuantReport.pdf>`_.


Normalization Plots
^^^^^^^^^^^^^^^^^^^^

The helper function :meth:`~mspypeline.BasePlotter.plot_all_normalizer_overview` is used to generate the same plot
multiple times with different normalizations methods of the base data.

.. _norm-overview:

Normalization overview
***********************
Created using: :meth:`~mspypeline.BasePlotter.plot_normalization_overview_all_normalizers` by calling
:meth:`~mspypeline.BasePlotter.plot_normalization_overview`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_normalization_overview

.. ipython:: python

    #plotter.plot_normalization_overview_all_normalizers("raw_log2", 0, save_path="./source/_static");
    print("skipping norm overview")

View `this normalization overview example <./_static/normalization_overview_all_normalizers_raw_log2.pdf>`_.

.. _heatmap-overview:

Heatmap overview
******************
Created using: :meth:`~mspypeline.BasePlotter.plot_heatmap_overview_all_normalizers` by calling
:meth:`~mspypeline.BasePlotter.plot_intensity_heatmap`.

    The Heatmap overview offers the opportunity to visually inspect how the distribution of protein intensities and missing
    values for each sample. As in the normalization overview, with this method, a separate plot is generated for each
    normalizer and attached to the document as another page. The heatmap overview can help to understand the differences
    between the distinct :ref:`protein intensity options <hyperparameter>` and :ref:`normalization methods <hyperparameter>`
    as it allows for instance to spot patterns between the different plot. The heatmap overview is based on the
    :meth:`~mspypeline.BasePlotter.plot_intensity_heatmap` method.

.. ipython:: python

    #plotter.plot_heatmap_overview_all_normalizers("raw_log2", 0, save_path="./source/_static");
    print("skipping heatmap overview")

View `this heatmap overview example <./_static/heatmap_overview_all_normalizers_raw_log2.pdf>`_.


Outlier detection and comparison plots
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Detection counts
*****************
Created using: :meth:`~mspypeline.BasePlotter.plot_detection_counts`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_detection_counts

.. ipython:: python

    @savefig detection_counts.png width=400 align=center
    plotter.plot_detection_counts("lfq_log2", 0, save_path=None);


Number of detected proteins
****************************
Created using: :meth:`~mspypeline.BasePlotter.plot_detected_proteins_per_replicate`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_detected_proteins_per_replicate

    Depending on whether :ref:`technical replicates <tech-reps>` should be recognized/averaged (top graph) or not
    (bottom graph) the data and resulting plot will have different outcomes. The number of detected proteins in total and
    per sample changes as 0 values are handled as missing values ("nan") and neglected when calculating the mean of samples.

.. ipython:: python

    @savefig detected_proteins.png width=600 align=center
    plotter.plot_detected_proteins_per_replicate("lfq_log2", 1, save_path=None);

    @savefig detected_proteins_tech_rep.png width=600 align=center
    plotter_with_tech_reps.plot_detected_proteins_per_replicate("lfq_log2", 1, save_path=None);


Venn diagrams
**************
Created using: :meth:`~mspypeline.BasePlotter.plot_venn_results`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_venn_results

.. ipython:: python

    plots = plotter.plot_venn_results("lfq_log2", 1, close_plots=None, save_path=None)

    @savefig venn_plot1.png width=600 align=center
    select_fig(plots, 0);

    @savefig venn_plot2.png width=450 align=center
    select_fig(plots, 1);


Group diagrams
***************
Created using: :meth:`~mspypeline.BasePlotter.plot_venn_groups`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_venn_groups

.. ipython:: python

    plots = plotter.plot_venn_groups("lfq_log2", 0, close_plots=None, save_path=None);

    @savefig venn_group_plot1.png height=300 width=540 align=left
    select_fig(plots, 0);

    @savefig venn_group_plot2.png height=300 width=100 align=right
    select_fig(plots, 1);


PCA overview
*************
Created using: :meth:`~mspypeline.BasePlotter.plot_pca_overview`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_pca_overview

.. ipython:: python

    @savefig pca_overview.png width=450 align=center
    plotter.plot_pca_overview("lfq_log2", 1, save_path=None);


Intensity histogram
********************
Created using: :meth:`~mspypeline.BasePlotter.plot_intensity_histograms`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_intensity_histograms

.. ipython:: python

    @savefig intensity_hist.png width=600 align=center
    plotter.plot_intensity_histograms("lfq_log2", 1, save_path=None);


Relative std
*************
Created using: :meth:`~mspypeline.BasePlotter.plot_relative_std`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_relative_std

.. ipython:: python

    @savefig relative_std.png width=500 align=center
    plotter.plot_relative_std("lfq_log2", 0, save_path=None);


Scatter replicates
*******************
Created using: :meth:`~mspypeline.BasePlotter.plot_scatter_replicates`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_scatter_replicates

.. ipython:: python

    plots = plotter.plot_scatter_replicates("lfq_log2", 1, close_plots=None, save_path=None);

    @savefig scatter_replicates.png width=500 align=center
    select_fig(plots, 0);


Experiment comparison
**********************
Created using: :meth:`~mspypeline.BasePlotter.plot_experiment_comparison`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_experiment_comparison

.. ipython:: python

    @savefig experiment_comparison.png width=500 align=center
    plotter.plot_experiment_comparison("lfq_log2", 0, save_path=None);


Rank
*****
Created using: :meth:`~mspypeline.BasePlotter.plot_rank`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_rank

.. ipython:: python

    @savefig rank_plot.png width=500 align=center
    plotter.plot_rank("lfq_log2", 0, save_path=None);


Statistical inference plots
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pathway analysis
*****************
Created using: :meth:`~mspypeline.BasePlotter.plot_pathway_analysis`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_pathway_analysis

.. ipython:: python

    plots_level0 = plotter.plot_pathway_analysis("lfq_log2", 0, close_plots=None, save_path=None)
    plots_level1 = plotter.plot_pathway_analysis("lfq_log2", 1, close_plots=None, save_path=None)

    @savefig pathway_analysis_level0.png width=700 align=center
    select_fig(plots_level0, 0);

    @savefig pathway_analysis_level1.png width=700 align=center
    select_fig(plots_level1, 0);


Go analysis
************
Created using: :meth:`~mspypeline.BasePlotter.plot_go_analysis`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_go_analysis

.. ipython:: python

    @savefig go_analysis.png width=700 align=left
    plotter.plot_go_analysis("lfq_log2", 1, save_path=None);


Volcano plot (R)
*****************
Created using: :meth:`~mspypeline.BasePlotter.plot_r_volcano`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_r_volcano

.. ipython:: python

    @savefig volcano_plot.png width=500 align=center
    plotter_with_tech_reps.plot_r_volcano("lfq_log2", 0, sample1="H1975", sample2="H838", adj_pval=True, save_path=None);


Additionally via python
^^^^^^^^^^^^^^^^^^^^^^^

Kernel density estimate plot
*******************************
Created using: :meth:`~mspypeline.BasePlotter.plot_kde`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_kde

.. ipython:: python

    @savefig kde_lfq_plot.png width=320 align=left
    plotter.plot_kde("lfq_log2", 3, save_path=None);

    @savefig kde_raw_plot.png width=320 align=right
    plotter.plot_kde("raw_log2", 3, save_path=None);

Boxplot
********
Created using: :meth:`~mspypeline.BasePlotter.plot_boxplot`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_boxplot

.. ipython:: python

    @savefig boxplot.png width=700 align=center
    plotter.plot_boxplot("lfq_log2", 3, save_path=None);

Number of Proteins vs Quantiles
********************************
Created using: :meth:`~mspypeline.BasePlotter.plot_n_proteins_vs_quantile`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_n_proteins_vs_quantile


.. ipython:: python
    :okwarning:

    @savefig n_proteins_vs_quantile.png width=700 align=center
    plotter.plot_n_proteins_vs_quantile("lfq_log2", 3, save_path=None);

Intensity Heatmap
******************
Created using: :meth:`~mspypeline.BasePlotter.plot_intensity_heatmap`

.. autodescriptiononly:: mspypeline.BasePlotter.plot_intensity_heatmap

.. ipython:: python

    @savefig intensity_heatmap.png width=700 align=center
    plotter.plot_intensity_heatmap("lfq_log2", 3, save_path=None);
