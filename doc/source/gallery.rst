.. _gallery:

Gallery
========

In the following, each plot that can be created with the ``mspypeline`` will be shown and explained giving a minimal
code example to create this plot when using the package as a :ref:`python module <python-quickstart>`.

Quality Control Report
^^^^^^^^^^^^^^^^^^^^^^^
Given that :ref:`all MaxQuant files <file-readers>` which are used to create the report are provided, the resulting
Quality Control Report will be a multi-page pdf document composed of a variety of information and graphics.
Several examples of diagrams from the report are given below, a complete report can be viewed in more detail here (ref).

--> include figure


Normalization Plots
^^^^^^^^^^^^^^^^^^^^

.. _norm-overview:

Normalization overview
***********************
Created using: :meth:`~mspypeline.BasePlotter.plot_normalization_overview`

The Normalization overview offers the opportunity to examine different aspects of the data in three distinct plots. For
each :ref:`normalization method <hyperparameter>` provided an additional page will be attached to the resulting pdf file starting with the
raw or not normalized data. That way it is possible to get a better understanding of the effects of the normalization
methods on the data, to inspect the different approaches and to find the best suitable normalization for the data. The
normalization overview combines the plots :meth:`~mspypeline.BasePlotter.plot_kde`,
:meth:`~mspypeline.BasePlotter.plot_n_proteins_vs_quantile` and :meth:`~mspypeline.BasePlotter.plot_boxplot`.



.. _heatmap-overview:

Heatmap overview
******************
Created using: :meth:`~mspypeline.BasePlotter.plot_heatmap_overview_all_normalizers`

The Heatmap overview offers the opportunity to visually inspect how the distribution of protein intensities and missing
values for each sample. As in the normalization overview, with this method, a separate plot is generated for each
normalizer and attached to the document as another page. The heatmap overview can help to understand the differences
between the distinct :ref:`protein intensity options <hyperparameter>` and :ref:`normalization methods <hyperparameter>`
as it allows for instance to spot patterns between the different plot. The heatmap overview is based on the
:meth:`~mspypeline.BasePlotter.plot_intensity_heatmap` method.





Outlier detection and comparison plots
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Detection counts
*****************
Created using: :meth:`~mspypeline.BasePlotter.plot_detection_counts`

This bar diagram shows how often proteins were detected in a number of replicates for each group.


Number of detected proteins
****************************
Created using: :meth:`~mspypeline.BasePlotter.plot_detected_proteins_per_replicate`

This bar plot shows the number of detected proteins per sample as well as the total number of detected proteins for each
group of a selected level. The average number of detected proteins is indicated as gray dashed line. Depending on
whether :ref:`technical replicates <tech-reps>` should be recognized/averaged (top graph) or not (bottom graph) the data
and resulting plot will have different outcomes. The number of detected proteins in total and per sample changes as 0
values are handled as missing values ("nan") and neglected when calculating the mean of samples.



Venn diagrams
**************
Created using: :meth:`~mspypeline.BasePlotter.plot_venn_results`

With this plotting method, both a venn diagram and a bar-venn diagram are created comparing the similarity of the
replicates of each group from the selected level (based on protein counts).
The venn diagram depicts the number of identified proteins per replicate/set as circle. Overlapping areas indicate the
number of detected proteins that are shared between the overlapping replicates. Non overlapping areas indicate the
number of proteins uniquely found in the replicate/set. A maximum of 3 replicates per group can be compared in the venn
diagram.



The bar-venn diagrams allow for the comparison of more than 3 replicates. The plot consists of two graphs, an upper bar
diagram, tha indicates the number of unique or shared proteins of a set or overlapping sets. The lower graph indicates
which set or sets are being compared, respectively, which protein count (upper graph) belongs to which comparison
(lower graph). An example of a bar-venn diagram is shown in the paragraph below (ref group diagrams).



Group diagrams
***************
Created using: :meth:`~mspypeline.BasePlotter.plot_venn_groups`

With this plotting method, both a venn diagram and a bar-venn diagram are created comparing the similarity of the
replicates of each group from the selected level (based on protein counts).
The venn diagram depicts the number of identified proteins per replicate/set as circle. Overlapping areas indicate the
number of detected proteins that are shared between the overlapping replicates. Non overlapping areas indicate the
number of proteins uniquely found in the replicate/set. A maximum of 3 replicates per group can be compared in the venn
diagram. An example of a venn diagram is shown in the paragraph above (see ref venn).
The bar-venn diagrams allow for the comparison of more than 3 replicates. The plot consists of two graphs, an upper bar
diagram, tha indicates the number of unique or shared proteins of a set or overlapping sets. The lower graph indicates
which set or sets are being compared, respectively, which protein count (upper graph) belongs to which comparison
(lower graph).

.. note::
    To determine which proteins can be compared between the groups and which are unique for one group an
    internal :ref:`threshold function <thresholding>` is applied.



PCA overview
*************
Created using: :meth:`~mspypeline.BasePlotter.plot_pca_overview`

This plotting method creates a PCA plot comparing all components against each other. The default is 2 components where
only PC 1 and PC 2 are compared. The PCA results do not change in dependence on the chosen level, however, determining
the level on which the data should be compared influences the coloring of the scatter elements. Each group of the
selected level is colored differently. Multiple different analysis options can be chosen to generate a PCA
(see: :ref:`multiple option config <default-yaml>`).




Intensity histogram
********************
Created using: :meth:`~mspypeline.BasePlotter.plot_intensity_histograms`

For each group of the selected level a histogram is created that counts the occurrence of the binned intensity values of
each sample. If *"show_mean"* is set to True in the :ref:`configs <default-yaml>` the mean intensity of the plotted
samples of a group will be shown as gray dashed line.



Relative std
*************
Created using: :meth:`~mspypeline.BasePlotter.plot_relative_std`

One plot per group of the selected level is created with the relative standard deviation of each protein between the
samples of a group. Low deviation shows that measured intensities are stable over multiple samples.

.. note::
    To determine which proteins can be compared between the two samples an internal :ref:`threshold function
    <thresholding>` is applied.



Scatter replicates
*******************
Created using: :meth:`~mspypeline.BasePlotter.plot_scatter_replicates`

With this plotting methos, for each group of the selected level, pairwise scatter comparisons of all replicates of a
group are plotted above each other in one graph (based on protein intensities). Unique proteins per replicate are shown
at the bottom and right side of the graph (substitution of na values by min value of data set). Pearsons's Correlation
Coefficient r² is given in the legend and calculated based on proteins of diagonal scatter/proteins that have a non na
value in both samples compared.



Experiment comparison
**********************
Created using: :meth:`~mspypeline.BasePlotter.plot_experiment_comparison`

Creates a pairwise scatter comparison between each combination of the groups of the selected level (based on protein
intensities). For each comparison a new plot is created. Unique proteins per replicate are shown at the bottom and right
side of the graph (substitution of na values by min*0.95 value of sample data set). Pearsons's Correlation Coefficient
r² is given in the legend and calculated based on proteins of diagonal scatter/proteins that have a non na value in both
samples.

.. note::
    To determine which proteins can be compared between the two groups and which are unique for one group an
    internal :ref:`threshold function <thresholding>` is applied.



Rank
*****
Created using: :meth:`~mspypeline.BasePlotter.plot_rank`

With this method one plot per group (mean intensity of samples) is created, where all proteins are sorted by intensity
value and plotted against their rank. The highest intensity accounts for rank 0 the lowest intensity for the number of
proteins - 1 whereby proteins with missing values are neglected. The median intensity of all proteins is given in the
legend. Additionally, if a protein is part of a selected :ref:`pathway <pathway-proteins>` it will be presented in color
and the median rank of all proteins of a given pathway is indicated. Multiple pathways can be selected and will be
represented in the same graph as distinct groups.




Statistical inference plots
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pathway analysis
*****************
Created using: :meth:`~mspypeline.BasePlotter.plot_pathway_analysis`

With the pathway analysis, two plots per selected pathway are created, one indicating significances and the other
without showing significances. For each protein of the pathway a subplot is created displaying the intensities of the
protein for all groups and significances are calculated for each pairwise comparison between groups with an independent
`t-test <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html>`__. For a group of multiple
samples, the protein intensity per sample is shown as a single scatter dot colored per group.

.. note::
    To determine which proteins can be compared between two groups an internal :ref:`threshold function
    <thresholding>` is applied.



Go analysis
************
Created using: :meth:`~mspypeline.BasePlotter.plot_go_analysis`

With this plotting method, an enrichment analysis for each selected :ref:`GO Term file <go-term-proteins>` is created
(based on protein counts).
First, for each GO term list a list *"pathway_genes"* is created by taking the intersection of the proteins from the GO
list and the total detected proteins. Secondly, a list of *"non_pathway_genes"* is created which comprises total
detected proteins but proteins in *"pathway_genes"*. Third, a list of *"experiment_genes"* and *"non_experiment_genes"*
is created in a similar fashion where an experiment references to a sample/group of samples of the data set.
Lastly, a `fisher exact test <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html>`__ is
calculated with the following contingency table and the "greater" alternative.

+------------------------+--------------------------------------+------------------------------------------+
|                        | in pathway                           | not in pathway                           |
+========================+======================================+==========================================+
| **in experiment**      | experiment_genes & pathway_genes     | experiment_genes & not_pathway_genes     |
+------------------------+--------------------------------------+------------------------------------------+
| **not in experiment**  | not_experiment_genes & pathway_genes | not_experiment_genes & not_pathway_genes |
+------------------------+--------------------------------------+------------------------------------------+

The resulting p-value is thus, also dependent on the overall protein count of the sample/group of samples.


Volcano plot (R)
*****************
Created using: :meth:`~mspypeline.BasePlotter.plot_r_volcano`

When using the volcano plot method, two plots for each pairwise comparison of the groups of the selected level (min 3
samples per group required) where one plot has a set of proteins annotated and the other does not. The volcano plot
shows the log2 fold change between the two different conditions against the -log10(p value) (based on protein
intensities). The p value is determined using the R limma package (`moderated t-statistic
<https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf>`__). A p value and fold change
cutoff are applied and all proteins below the cutoff are considered non significant. Additionally, the intensities of
unique proteins of both conditions are shown next to the volcano plot.

.. note::
    To determine which proteins can be compared between the two groups and which are unique for one group an internal
    :ref:`threshold function <thresholding>` is applied.




Additionally via python
^^^^^^^^^^^^^^^^^^^^^^^

Kernel density estimate plot
*******************************
Created using: :meth:`~mspypeline.BasePlotter.plot_kde`

In the Kernel density estimate (KDE) plot, one density graph per sample is plotted indicating the Intensity on the x
axis and the density on the y axis. The KDE is part of the :ref:`Normalization overview <norm-overview>`.


Boxplot
********
Created using: :meth:`~mspypeline.BasePlotter.plot_boxplot`

Creates one boxplot per group of the selected level sorted by median intensity. The boxplot is part of the
:ref:`Normalization overview <norm-overview>`.


Number of Proteins vs Quantiles
********************************
Created using: :meth:`~mspypeline.BasePlotter.plot_n_proteins_vs_quantile`

This plot shows the protein intensities against the number of identified proteins per sample. The samples are indicated
as a horizontal line of scatter dots where the color anf x position of a dot indicate the intensity value of the
respective quantile. The y position of the dots of a sample indicate the total number of detected proteins in that
sample. Solid, rather vertical lines indicate a linear fit of each quantile for all the samples. This plot is part of
the :ref:`Normalization overview <norm-overview>`.


Intensity Heatmap
******************
Created using: :meth:`~mspypeline.BasePlotter.plot_intensity_heatmap`

This heatmap is showing protein intensities, where samples are given in rows on the y axis and proteins on the x axis.
Missing values are colored in gray. The heatmap can be used to spot patterns in the different
:ref:`normalization methods <hyperparameter>` and to understand how different :ref:`intensity types <hyperparameter>`
affect the data. The :ref:`Heatmap-overview <heatmap-overview>` is created from a series of these intensity heatmap
plot.