.. _workflow:

Workflow
========

Quality control
~~~~~~~~~~~~~~~
Ensure that the quality of the results to analyze are of sufficient quality. This can be done with the report provided
by the :ref:`mspypeline plotters <plotters>`.

Determining normalization method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Decide for a normalization method. Normalization options are:

* No normalization
* Median Normalization via: :class:`~mspypeline.MedianNormalizer`
* Quantile Normalization with missing value handling via: :class:`~mspypeline.QuantileNormalizer`
  and :func:`~mspypeline.interpolate_data`
* Tail Robust Quantile Normalization (TRQN) via: :class:`~mspypeline.TailRobustNormalizer` and
  :class:`~mspypeline.QuantileNormalizer`
* TRQN with missing value handling via: same as above and :func:`~mspypeline.interpolate_data`
* Tail Robust Median Normalization via: :class:`~mspypeline.TailRobustNormalizer` and
  :class:`~mspypeline.MedianNormalizer`

To help to decide for a normalization method two plots can help:
:meth:`~mspypeline.BasePlotter.plot_normalization_overview` and :meth:`~mspypeline.BasePlotter.plot_intensity_heatmap`.
Please read the function description how normalized data should look like. Also, after deciding for a normalization
method it should not be changed, especially in hindsight after seeing the results of other plots.

(Optionally) Select pathways and GO-Terms of interest
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Select :ref:`pathway-proteins`. Selected pathways will have following effects:

* for the :meth:`~mspypeline.BasePlotter.plot_pathway_analysis` one plot per pathway will be created
* in the :meth:`~mspypeline.BasePlotter.plot_rank`, if a protein is found it will be marked on the plot
  and colored by the pathway

Select :ref:`go-term-proteins`. Selected GO-Terms will have following effects:

* for the :meth:`~mspypeline.BasePlotter.plot_go_analysis` one additional barplot is added

Create descriptive and comparison plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The descriptive and comparison plots can e.g. help to analyze how biological replicates compare to another or
how different conditions effect detected proteins.

Create statistical inference plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The statistical inference plots show how significant different proteins are expressed in different conditions.
The Volcano plot uses the R package limma, so additional R packages might be downloaded when this plot is created
for the first time.
