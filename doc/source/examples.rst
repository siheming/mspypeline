.. _examples:

Examples
========

Create a plotter from custom dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When data is available is some other format or already preprocessed it is possible to create a plotter
by providing the data in a DataFrame.

.. ipython:: python

    import pandas as pd
    import numpy as np
    from mspypeline import BasePlotter
    from mspypeline.helpers import get_analysis_design
    samples = [f"Tumor_Stage{s}_Experiment{e}" for s in ["I", "II", "III", "IV"] for e in range(1, 4)] + [
        f"Control_Stage{s}_Experiment{e}" for s in ["I", "II", "III", "IV"] for e in range(1, 4)]

    data = pd.DataFrame(np.exp2(np.random.normal(26, 3, (100, 24))).astype(int), columns=samples)
    data.iloc[:, 12:] = data.iloc[:, 12:] + 1e8  # this is just for later, not required here
    analysis_design = get_analysis_design(samples)
    plotter = BasePlotter("result_dir", reader_data={"custom_reader": {"my_data": data}},
        intensity_df_name="my_data", configs={"analysis_design": analysis_design},
        required_reader="custom_reader", intensity_entries=[("raw", "", "Intensity")])
    plotter.all_tree_dict["raw_log2"].aggregate(None, None).head()
    plotter.all_tree_dict["raw_log2"].groupby(0).head()

Or alternatively this is also possible:

.. ipython:: python

    plotter = BasePlotter("result_dir", configs={"analysis_design": analysis_design})
    plotter.add_intensity_column("raw", "", "Intensity", df=data)
    plotter.all_tree_dict["raw_log2"].aggregate(None, None).head()
    plotter.all_tree_dict["raw_log2"].groupby(0).head()


Add normalization options
~~~~~~~~~~~~~~~~~~~~~~~~~
Normalizing data is easy by adding a normalized option. The source data needs to be provided, then a Normalizer
and lastly a name for the new option.

.. ipython:: python

    from mspypeline import MedianNormalizer
    plotter.add_normalized_option("raw", MedianNormalizer, "median_norm")
    plotter.all_tree_dict.keys()

This added two new options: 'raw_median_norm' and 'raw_median_norm_log2'.

Global protein level
~~~~~~~~~~~~~~~~~~~~
This is an example plot which compares the global protein level and additionally shows mean and standard deviation.

.. ipython:: python

    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    def create_plot(intensity):
        df_grouped_normal = plotter.all_tree_dict[intensity]["Control"].groupby()
        df_grouped_tumor = plotter.all_tree_dict[intensity]["Tumor"].groupby()
        segs_normal = [[(i, value) for i, value in enumerate(df_grouped_normal.loc[protein])]
            for protein in df_grouped_normal.index]
        linecoll_normal = LineCollection(segs_normal, color="gray", linewidth=0.05)
        segs_tumor = [[(i, value) for i, value in enumerate(df_grouped_tumor.loc[protein])]
            for protein in df_grouped_tumor.index]
        linecoll_tumor = LineCollection(segs_tumor, color="lightcoral", linewidth=0.05)
        fig, ax = plt.subplots(1, 1, figsize=(14, 7))
        ax.add_collection(linecoll_normal);
        ax.add_collection(linecoll_tumor);
        ax.errorbar([i for i, x in enumerate(df_grouped_normal.columns)], df_grouped_normal.mean(),
            yerr=df_grouped_normal.std(), color="black");
        ax.errorbar([i for i, x in enumerate(df_grouped_tumor.columns)], df_grouped_tumor.mean(),
            yerr=df_grouped_tumor.std(), color="red");
        ax.set_ylim(18, 40);
        ax.set_xlim(-0.5, 3.5);
        ax.set_xticks([0, 1, 2, 3]);
        ax.set_xticklabels(df_grouped_tumor.columns);
        ax.set_ylabel(intensity);
        return fig

    @savefig plot_global_protein_level_raw.png width=6in
    create_plot("raw_log2");

    @savefig plot_global_protein_level_median_norm.png width=6in
    create_plot("raw_median_norm_log2");



.. missing_tumor = plotter.all_tree_dict[intensity]["Tumor"].groupby(lambda x: x.isnull().sum())