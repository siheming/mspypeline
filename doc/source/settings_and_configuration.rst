.. _settings:

Settings and Configurations
============================

The configuration YAML file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
| The YAML file stores the main configurations that determine which results are created. This configuration file
  also offers the specification of multiple optional analysis options.
| A :ref:`default YAML file <default-yaml>`, containing the default analysis settings, is provided by ``mspypeline``
  at the start of the analysis. This file can be edited to further individualize the results.
  If the data analysis is performed via the GUI, no further interaction with the YAML file is necessary.
| Plots can be easily reproduced by reusing the same settings/YAML file.

.. _default-yaml:

Analysis settings
*************************
.. toctree::
    default_config
.. toctree::
    adjustable_options_config


.. _analysis-design:

Analysis Design
~~~~~~~~~~~~~~~
| To perform comparative data analysis, ``mspypeline`` assumes that data consists of samples that can be arranged
  into a tree structure resembling the experimental setup. Different samples of an experiment are arranged
  in groups and subgroups dependent on the sample's name. This naming convention is the key principle to draw comparisons
  between distinct samples of different groups/at different levels. The analysis design can be of any level of depth.

.. warning::
   * Different levels of the analysis design need to be separated by an underscore (**_**)
   * All samples must have the same number of levels (meaning same number of underscores)


**Example Analysis Design**

| Assume an experiment, where two cancer cell lines are compared to two non-cancer cell lines. Additionally, we have 3
  technical replicates per cell line. One sample could be called e.g. "Cancer_Line1_Rep1", or "Control_Line2_Rep3".
| The resulting analysis design could look like this:

.. ipython:: python

    from mspypeline.helpers import get_analysis_design
    from pprint import pprint
    samples = [f"Cancer_Line{l}_Rep{r}" for l in range(1, 3) for r in range(1, 4)] + [
        f"Control_Line{l}_Rep{r}" for l in range(1, 3) for r in range(1, 4)]
    analysis_design = get_analysis_design(samples)
    pprint(analysis_design)

| Level 0 corresponds to "Cancer" and "Control", level 1 corresponds to "Cancer_Line1", "Cancer_Line2", "Control_Line1"
  and "Control_Line2" and level 2 corresponds to the technical replicates "Cancer_Line1_Rep1", ... . Here, the "Cancer"
  group has a total of 6 children, which are then split in "Line1" and "Line2". "Line1" and "Line2" both have 3 children.
| If a comparison between the "Cancer" and "Control" groups should be performed, the lowest level (level 0) of the Data
  Tree must be chosen for the analysis. Consequently, the results show the comparison of both level 0 groups "Cancer"
  and "Control" with all their children.
| On the other hand, if the different cancer and control cell lines should be compared, the second level (level 1) of the
  Data Tree must be selected. The results will consequently show the comparison of the four different level 1 groups
  "Cancer_Line1", "Cancer_Line2", "Control_Line1 and "Control_Line2", each with their 3 replicates.
| If all the 12 replicates should be compared to each other, the last level (level 2) must be selected in the analysis.


.. _sample-mapping:

Sample Mapping
***************

| Should the naming convention deviate from the expected standard, it is possible to subsequently correct the sample
  naming with a provided sample mapping file so that samples translate to a proper analysis design. Some examples for
  potential reasons for naming convention violation are given in the table below.
| If the naming convention is violated a sample mapping can be provided manually or by using the default file.

* default **sample_mapping_template.txt** file: This file is created automatically if the naming convention is
  incorrect. The file already provides the general structure of the sample mapping consisting of two columns *old name*
  which is readily filled out and *new name* which needs to be filled with the new desired sample name. Then the file
  needs to be renamed to **sample_mapping.txt**.
* manual **sample_mapping.txt** file: The file needs to be tab-separated with two columns and saved on the same level
  as the config directory. The first column named *old name* should contain the sample name
  of the ms run. The second column named *new name* should follow the naming convention.

| A *sample-mapping.txt* file should have a simple structure as shown in the table below. All samples need to be mapped,
  non-conforming sample names need to be corrected in the *new name* column (e.g. row 1-3) and sample names that are
  already correct simply need to be copied in the *new name* column (e.g. row 4-5).
| Since this example is only intended to give a small insight into the sample-mapping file, the remaining samples of the
  analysis design are not listed here, which would be necessary for a real data analysis.

.. list-table:: sample_mapping.txt
    :widths: 25, 25
    :header-rows: 1

    * - old name
      - new name
    * - Cancer_Line1-Rep1
      - Cancer_Line1_Rep1
    * - Cance_Line1_Rep2
      - Cancer_Line1_Rep2
    * - Cancer_Line_1_Rep3
      - Cancer_Line1_Rep3
    * - Cancer_Line2_Rep1
      - Cancer_Line2_Rep1
    * - Cancer_Line2_Rep2
      - Cancer_Line2_Rep2
    * - Cancer...
      - Cancer..

.. _tech-reps:

Technical Replicates
********************
| If the :ref:`configuration setting <default-yaml>` `has_techrep` is set to True or the corresponding checkbox in the
  GUI is ticked, the highest level of the analysis design is considered technical replicates.
| Technical replicates are averaged and averaged to one sample of the next highest level in the analysis design.
  Respectively, the mean of all samples below a node is calculated and assigned to that node. The last level of the Data
  Tree is thus omitted. Values that are 0 are replaced with missing values, which are neglected when calculating the
  mean of samples (e.g. the average of the three values 32, 30, 0 would be replaced with 32, 30 and NaN resulting in an
  average of 31).

.. ipython:: python

    import numpy as np
    import pandas as pd
    from mspypeline import DataTree
    data = pd.DataFrame(np.exp2(np.random.normal(26, 2, (3, 12))).astype(int), columns=samples)
    tree_agg = DataTree.from_analysis_design(analysis_design, data, True)
    tree_no_agg = DataTree.from_analysis_design(analysis_design, data, False)
    tree_no_agg.aggregate(None, None)
    tree_agg.aggregate(None, None)

| Columns of the latter output are now named "Cancer_Line1", "Cancer_Line2", etc. and the values of the replicates
  "Rep1" to "Rep3" are averaged to a mean. This procedure can help to improve data reproducibility since measurement
  results can be quite noisy and/or proteins might be missing in some of the samples by random chance.



.. _thresholding:

Thresholds and Comparisons
~~~~~~~~~~~~~~~~~~~~~~~~~~
| Several analysis methods require the determination whether a detected protein can be compared between two groups A and B.
| For group comparisons, the protein counts or intensities for all samples of a group are cumulated and averaged whereby
  missing values are omitted in the calculation. In mass spectrometry data, missing values are frequently observed and
  can be found in multiple samples for one protein. To ensure appropriate data analysis a sufficient number of samples
  with non-missing values per group must be provided. The required number of samples per group, the threshold, is
  thereby individually determined per group and dynamically controlled based on the number of samples in the group. The
  MSPypeline provides an internal thresholding function (Fig. 2.6), however any other desired function may be applied.
| Besides the determination of comparable proteins, several parts of the data analysis further distinguish proteins that
  are unique for/ exlusively found in a group. There are four potential scenarios of categorizing the protein:

* Unique in A: above threshold in A and completely absent in B
* Unique in B: above threshold in B and completely absent in A
* Can be compared: above threshold in A and B
* Otherwise: not considered

| Thresholding is important for the venn group diagrams, the relative standard deviation graph, the group comparison
  scatter plot and for the volcano plot.

.. ipython:: python

    import matplotlib.pyplot as plt
    from mspypeline.helpers import get_number_of_non_na_values, get_non_na_percentage
    x_values = [x for x in range(1, 31)]
    y_values = [get_number_of_non_na_values(x) for x in x_values]
    fig, ax = plt.subplots(1,1, figsize=(7,5))
    ax.plot(x_values, y_values, marker=".");
    ax.set_xticks(x_values);
    ax.set_yticks(y_values);
    ax.set_xlabel("Number of samples");
    @savefig plot_non_na_samples.png width=6in
    ax.set_ylabel("Required number of non na (missing) values");

Starting with a minimum **number** of 3, the **number** of non missing values to function as threshold increases
steadily with rising numbers of samples per group.


An example: Group A has 7 samples, Group B has 8 Samples.

* Unique in A: Group A has equals or more than 5 non missing values and Group B has only missing values
* Unique in B: Group B has equals or more than 6 non missing values and Group A has only missing values
* Can be compared: Group A has equals or more than 5 non missing values and
  Group B has equals or more than 6 non missing values
* Not considered: In all other cases

This threshold criterion is quite harsh, but the results will be reliable.

The next plot shows the required **percentage** of non zero values as a function of the number of samples in a group.

.. ipython:: python

    y_values = [get_non_na_percentage(x) for x in x_values]
    fig, ax = plt.subplots(1,1, figsize=(7,5))
    ax.plot(x_values, y_values, marker=".");
    ax.set_ylim(0, 1);
    ax.set_xlabel("Number of samples");
    ax.set_xticks(x_values);
    @savefig plot_non_zero_values.png width=6in
    ax.set_ylabel("Required percentage/100 of non zero values");
