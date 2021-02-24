.. _workflow:

Workflow
========

Quality control
~~~~~~~~~~~~~~~
| Ensure that the quality of the results to analyze are of sufficient quality. This can be done with the report provided
  by the :ref:`MaxQuant plotter <plotters>`.

.. _hyperparameter:

Data Preprocessing
~~~~~~~~~~~~~~~~~~~
| Data may be processed in multiple ways and such modeling has the potential to shape the results of an analysis
  substantially.
| This kind of data preprocessing comprises the choice of protein intensities provided by MaxQuant, including raw,
  label-free quantification (LFQ) or intensity-based absolute quantification (iBAQ) intensities, averaging technical
  replicates, the removal of erroneous samples or the normalization and standardization of the data set.

Intensity options
******************

* LFQ Intensity ("lfq_log2")
* raw Intensity ("raw_log2")
* iBAQ Intensity ("ibaq_log2")

| Despite the choice of protein intensity, the GUI handles all data in log2 format. However, it is possible to analyze
  the data without log2 scale ("lfq", "raw", "ibaq") if advanced data analysis is performed by interacting rather direct
  with the plotters.


Normalization options
*********************

* No normalization
* Median Normalization via: :class:`~mspypeline.MedianNormalizer`
* Quantile Normalization with missing value handling via: :class:`~mspypeline.QuantileNormalizer`
  and :func:`~mspypeline.interpolate_data`
* Tail Robust Quantile Normalization (TRQN) via: :class:`~mspypeline.TailRobustNormalizer` and
  :class:`~mspypeline.QuantileNormalizer`
* TRQN with missing value handling via: same as above and :func:`~mspypeline.interpolate_data`
* Tail Robust Median Normalization via: :class:`~mspypeline.TailRobustNormalizer` and
  :class:`~mspypeline.MedianNormalizer`

| To aid the determination of the best possible normalization method, two plots may be created:
  :meth:`~mspypeline.BasePlotter.plot_normalization_overview` and
  :meth:`~mspypeline.BasePlotter.plot_heatmap_overview_all_normalizers`.
| These methods will output a multipage PDF file in which the data is plotted repeatedly after applying the different
  normalization options. Thereby it is possible to get a better understanding of the effect of each normalization method
  on the data.
| Please read the function description explaining what normalized data should look like. Once a normalization method is
  chosen, it is highly recommended to perform all further analysis with the same normalized data.



Exploratory Analysis
~~~~~~~~~~~~~~~~~~~~~

Create outlier detection and comparison plots
**********************************************
| The descriptive and comparison plots can for example help to analyze how biological replicates compare to another or
  how different conditions effect detected proteins.

Create statistical inference plots
**********************************
| Statistical inference plots can inform about differential protein intensities between groups of the data set.
  The calculation of statistical significances of the protein intensity deviation of groups can help to exploit
  biological questions by incorporation the functional profile of proteins or protein sets.
| Statistics for each plot are calculated based on the intended proposition of the plot.

* for the :meth:`~mspypeline.BasePlotter.plot_pathway_analysis` an independent t-test is applied
* for the :meth:`~mspypeline.BasePlotter.plot_go_analysis` a fisher'S exact test is applied
* for the :meth:`~mspypeline.BasePlotter.plot_r_volcano` plot the moderated t-statistics is applied which is
  implemented by the R package limma. Additional R packages might
  be downloaded when this plot is created for the first time.


Select pathways and GO-Terms of interest
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Select :ref:`pathway-proteins`. Selected pathways will have following effects:

* for the :meth:`~mspypeline.BasePlotter.plot_pathway_analysis` one plot per pathway will be created
* in the :meth:`~mspypeline.BasePlotter.plot_rank`, if a protein is found it will be marked on the plot
  and colored by the pathway
* in the :meth:`~mspypeline.BasePlotter.plot_r_volcano`, if a pathway is selected, proteins of that pathway will be
  annotated in the plot instead of the most significant proteins that are annotated by default

Select :ref:`go-term-proteins`. Selected GO-Terms will have following effects:

* for the :meth:`~mspypeline.BasePlotter.plot_go_analysis` one additional barplot is added per GO term

