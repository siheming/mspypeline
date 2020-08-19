.. py:currentmodule:: mspypeline

Package design
==============

.. _analysis-design:

Analysis design
~~~~~~~~~~~~~~~

The project assumes that samples can be arranged into a tree structure which represents an experimental design.
The analysis design can have any level of depth.

.. warning::
   * Different levels of the analysis design need to be separated by underscore (**_**)
   * All samples must have the same number of levels (meaning same number of underscores)

**Example Analysis Design**

Assume an experiment where two cancer cell lines are compared to two non-cancer cell lines. Additionally, we have 5
technical replicates per cell line. One sample could be called e.g. "Cancer_Line1_Rep1", or "Control_Line2_Rep3".
The resulting analysis design could look like this:

.. ipython:: python

    from mspypeline.helpers import get_analysis_design
    from pprint import pprint
    samples = [f"Cancer_Line{l}_Rep{r}" for l in range(1, 3) for r in range(1, 6)] + [
        f"Control_Line{l}_Rep{r}" for l in range(1, 3) for r in range(1, 6)]
    analysis_design = get_analysis_design(samples)
    pprint(analysis_design)

Here, the "Cancer" group has a total of 10 children, which are then split in "Line1" and "Line2". "Line1" and "Line"
both have 5 children. Level 0 corresponds to "Cancer" and "Control", level 1 corresponds to "Cancer_Line1",
"Cancer_Line2", "Control_Line1" and "Control_Line2".

An analysis design could also be more complicated like: "Cancer_Human_Week10_Cell1_Rep1".


Sample mapping
^^^^^^^^^^^^^^
If the naming in the file does not translate to a proper analysis design an additional mapping can be provided.

* This file is created automatically if the naming convention is violated. If the file is created automatically it will
  be created with the name **sample_mapping_template.txt**. The general structure is already present, but the second
  column needs to be changed into the new desired name. Then the file needs to be renamed to **sample_mapping.txt**.
* Otherwise, it can be created manually with the name **sample_mapping.txt** on the same level as the config directory.
  The file needs to be tab-separated with two columns. The first column named old_name should contain the sample name
  of the ms run. The second column named new_name should follow the naming convention.

Technical Replicates
^^^^^^^^^^^^^^^^^^^^
If the configuration setting `has_techrep` is set to True or the corresponding checkbox in the GUI is ticked the
lowest level of the analysis design is considered as a technical replicate. Technical replicates are averaged to one
level higher in the analysis design by taking the mean of all samples below a node. Values that are 0 are replaced with
missing values, then for calculating the mean missing values are ignored (meaning the average of the three values
32, 30, 0 would be replaced with 32, 30 and NaN resulting in an average of 31).

.. ipython:: python

    import numpy as np
    import pandas as pd
    from mspypeline import DataTree
    data = pd.DataFrame(np.exp2(np.random.normal(26, 2, (5, 20))).astype(int), columns=samples)
    tree_agg = DataTree.from_analysis_design(analysis_design, data, True)
    tree_no_agg = DataTree.from_analysis_design(analysis_design, data, False)
    tree_no_agg.aggregate(None, None)
    tree_agg.aggregate(None, None)

The columns of the last output are now named "Cancer_Line1", "Cancer_Line2", etc. and the values of the replicates
"Rep1" to "Rep5" is averaged as mean. This can help to improve results since measurement results are noisy or proteins
might be missing in some of the samples by random chance.

Thresholds and Comparisons
^^^^^^^^^^^^^^^^^^^^^^^^^^
Some plots need to determine whether a protein can be compared between two different groups.
Four different results are possible when comparing a protein between group A and group B.

* Unique in A: above threshold in A and completely absent in B
* Unique in B: above threshold in B and completely absent in A
* Can be compared: above threshold in A and B
* Otherwise: not considered

The threshold is determined dynamically based on the number of samples in a group. The next two plots show the
thresholding for different numbers of samples.

.. ipython:: python

    import matplotlib.pyplot as plt
    from mspypeline.helpers import get_number_of_non_na_values, get_non_na_percentage
    x_values = [x for x in range(1, 31)]
    y_values = [get_number_of_non_na_values(x) for x in x_values]
    fig, ax = plt.subplots(1,1, figsize=(7,7))
    ax.plot(x_values, y_values, marker=".");
    ax.set_xticks(x_values);
    ax.set_yticks(y_values);
    ax.set_xlabel("Number of samples");
    @savefig plot_non_na_samples.png width=6in
    ax.set_ylabel("Required number of non na (missing) values");

The plot shows that the minimum number of non missing values is 3, then steadily increases.

An example: Group A has 7 samples, Group B has 8 Samples.

* Unique in A: Group A has equals or more than 5 non missing values and Group B has only missing values
* Unique in B: Group B has equals or more than 6 non missing values and Group A has only missing values
* Can be compared: Group A has equals or more than 5 non missing values and
  Group B has equals or more than 6 non missing values
* Not considered: In all other cases

This criterion is quite harsh, but the results will be dependable.

The next plot shows the required percentage of non zero values.

.. ipython:: python

    y_values = [get_non_na_percentage(x) for x in x_values]
    fig, ax = plt.subplots(1,1, figsize=(7,7))
    ax.plot(x_values, y_values, marker=".");
    ax.set_ylim(0, 1);
    ax.set_xlabel("Number of samples");
    ax.set_xticks(x_values);
    @savefig plot_non_zero_values.png width=6in
    ax.set_ylabel("Required percentage of non zero values");


File readers
~~~~~~~~~~~~
For the different file formats a file reader can be designed to translate the file into the correct internal format.
Currently only a reader for the MaxQuant format is provided.

Settings as YAML files
~~~~~~~~~~~~~~~~~~~~~~
The specified settings are saved as yaml files. Plots can be easily reproduced by reusing the same settings.
