import pandas as pd
import numpy as np
import os
import functools
import matplotlib.pyplot as plt
from scipy import stats
from itertools import combinations
from collections import defaultdict as ddict
import logging
import warnings
from typing import Dict, Type, Iterable, Optional, Union, Any, List, Tuple
from sklearn.decomposition import PCA
from copy import deepcopy

from mspypeline.core import MSPInitializer
from mspypeline.file_reader import BaseReader
from mspypeline.plotting_backend import matplotlib_plots
from mspypeline.modules import default_normalizers, Normalization, DataTree
from mspypeline.helpers import get_number_of_non_na_values, get_intersection_and_unique, \
    get_logger, dict_depth, add_end_docstrings, make_contrasts


def validate_input(f):
    """
    Ensures that the function input, dfs_to_use and levels are passed as lists.

    Parameters
    ----------
    f
        method that will be decorated
    """
    @functools.wraps(f)
    def wrapper(self, *args, **kwargs):
        assert len(args) <= 2
        dfs_to_use = kwargs.pop("dfs_to_use", None)
        levels = kwargs.pop("levels", None)

        # this is a convention, otherwise it wont work
        if dfs_to_use is None:
            dfs_to_use = args[0]
            if levels is None:
                levels = args[1]
        else:
            if levels is None:
                levels = args[0]

        if isinstance(dfs_to_use, str):
            dfs_to_use = dfs_to_use,
        if not isinstance(levels, Iterable):
            levels = levels,
        return f(self, dfs_to_use=dfs_to_use, levels=levels, **kwargs)
    return wrapper


plot_para_return_docstring = """

Parameters
----------
dfs_to_use
    which dataframes/intensities should be plotted
levels
    at which level of the data tree should the data be compared
kwargs
    accepts kwargs

Returns
-------
List
    A list of all created plots.

See Also
--------
{}

"""


class BasePlotter:
    """
    | Base plotter to create plots.
    | The two main methods of the Base plotter comprise *"get_"* functions to calculate and provide the data for the
      *"plot_"* functions. The latter incorporates the *"get_"* functions as well as functions from the matplotlib
      backend to combine data calculation, plotting and saving of the results in one method.
    """

    possible_plots = [
        "plot_detection_counts", "plot_detected_proteins_per_replicate", "plot_intensity_histograms",
        "plot_relative_std", "plot_rank", "plot_pathway_analysis",
        "plot_scatter_replicates", "plot_experiment_comparison", "plot_go_analysis", "plot_venn_results",
        "plot_venn_groups", "plot_r_volcano", "plot_pca_overview",
        "plot_normalization_overview_all_normalizers", "plot_heatmap_overview_all_normalizers"
    ]

    def __init__(
            self,
            start_dir: str,
            reader_data: Optional[Dict[str, Dict[str, pd.DataFrame]]] = None,
            intensity_df_name: str = "",
            interesting_proteins: Optional[Dict[str, pd.Series]] = None,
            go_analysis_gene_names: Optional[Dict[str, pd.Series]] = None,
            configs: Optional[dict] = None,
            required_reader: Optional[str] = None,
            intensity_entries: Tuple[str, str, str] = (),
            loglevel: int = logging.DEBUG
    ):
        """
        Parameters
        ----------
        start_dir
            location to save results
        reader_data
            mapping to provide input data
        intensity_df_name
            name/key to input data
        interesting_proteins
            mapping with pathway proteins to analyze
        go_analysis_gene_names
            mapping with go terms to analyze
        configs
            mapping of configuration
        required_reader
            name of the file reader
        intensity_entries
            tuple of (key in all_tree_dict, prefix in data, name in plot). See :meth:`add_intensity_column`.
        loglevel
            level of the logger
        """
        self.logger = get_logger(self.__class__.__name__, loglevel=loglevel)
        # general information
        # TODO make start dir optional
        self.start_dir = start_dir
        # read in optional arguments
        self.configs = {} if configs is None else deepcopy(configs)
        self.interesting_proteins = {} if interesting_proteins is None else interesting_proteins
        self.go_analysis_gene_names = {} if go_analysis_gene_names is None else go_analysis_gene_names
        self.normalizers = deepcopy(default_normalizers)
        self.selected_normalizer_name = self.configs.get("selected_normalizer", "None")
        self.selected_normalizer = self.normalizers.get(self.selected_normalizer_name, None)
        self.intensity_df = None
        if required_reader is not None:
            try:
                self.required_reader_data = reader_data[required_reader]
                self.intensity_df = self.required_reader_data[intensity_df_name]
                self.configs.update(self.configs.pop(required_reader, {}))
            except KeyError:
                self.logger.exception("Reader data does not provide information from %s reader", required_reader)
                raise

        # setup everything for all_intensity dict
        self.int_mapping = {}
        self.intensity_label_names = {}
        self.all_intensities_dict: Dict[str, pd.DataFrame] = {}
        self.all_tree_dict: Dict[str, DataTree] = {}
        self.analysis_design = self.configs.get("analysis_design", {})
        if not self.analysis_design:
            self.logger.warning("No analysis design was provided. Most plotting functions will not work")

        for option_name, name_in_file, name_in_plot in intensity_entries:
            self.add_intensity_column(option_name, name_in_file, name_in_plot)
            if self.selected_normalizer is not None:
                self.add_normalized_option(option_name, self.selected_normalizer, self.selected_normalizer_name)
                self.add_normalized_option(option_name, self.selected_normalizer, "normalized")

        # set all result dirs
        # path for venn diagrams
        self.file_dir_venn = os.path.join(self.start_dir, "venn")
        # path for descriptive plots
        self.file_dir_descriptive = os.path.join(self.start_dir, "outliers_detection_and_comparison")
        # path for normalization plots (heatmap and normalization overview)
        self.file_dir_normalization = os.path.join(self.start_dir, "normalization")
        # path for pathway analysis
        self.file_dir_pathway = os.path.join(self.start_dir, "pathway_analysis")
        # path for go analysis
        self.file_dir_go_analysis = os.path.join(self.start_dir, "go_analysis")
        # path for volcano plots
        self.file_dir_volcano = os.path.join(self.start_dir, "volcano")

    @classmethod
    def from_MSPInitializer(cls, mspinit_instance: MSPInitializer, **kwargs) -> "BasePlotter":
        """
        | Creates a BasePlotter from a :class:`~MSPInitializer`.

        Parameters
        ----------
        mspinit_instance
            instance of a :class:`~MSPInitializer` used to get correct inputs for the plotter.
        kwargs
            all kwargs, which are passed to the :func:`BasePlotter.__init__` can be overwritten by passing as kwargs.

        Returns
        -------
        BasePlotter
            functional plotter

        """
        default_kwargs = dict(
            start_dir=mspinit_instance.start_dir,
            reader_data=mspinit_instance.reader_data,
            intensity_df_name="",
            interesting_proteins=mspinit_instance.interesting_proteins,
            go_analysis_gene_names=mspinit_instance.go_analysis_gene_names,
            configs=mspinit_instance.configs,
            required_reader=None,
            intensity_entries=(),
            loglevel=mspinit_instance.logger.getEffectiveLevel()
        )
        default_kwargs.update(**kwargs)
        return cls(**default_kwargs)

    @classmethod
    def from_file_reader(cls, reader_instance: BaseReader, **kwargs):
        """
        | Creates a BasePlotter from a :class:`~BaseReader` (BasePlotter or MaxQuantPlotter).

        Parameters
        ----------
        reader_instance
            instance of a :class:`~BaseReader` used to get correct inputs for the plotter.
        kwargs
            all kwargs, which are passed to the :func:`BasePlotter.__init__` can be overwritten by passing as kwargs.

        Returns
        -------
        BasePlotter
            functional plotter

        """
        default_kwargs = dict(
            start_dir=reader_instance.start_dir,
            reader_data={reader_instance.name: reader_instance.full_data},
            intensity_df_name="",
            interesting_proteins=None,
            go_analysis_gene_names=None,
            configs=reader_instance.reader_config,
            required_reader=reader_instance.name,
            intensity_entries=(),
            loglevel=reader_instance.logger.getEffectiveLevel()
        )
        default_kwargs.update(**kwargs)
        return cls(**default_kwargs)

    def create_results(self):
        """
        Creates all plots that where chosen/set to True in the settings "create plot" (see :ref:`default-yaml`).
        """
        global_settings = self.configs.get("global_settings", {})
        self.logger.debug(f"got global settings: %s", global_settings)
        for plot_name in self.possible_plots:
            plot_settings_name = plot_name + "_settings"
            plot_settings = self.configs.get(plot_settings_name, {})
            plot_settings.update({k: v for k, v in global_settings.items() if k not in plot_settings})
            if plot_settings.pop("create_plot", False):
                self.logger.debug(f"creating plot {plot_name}")
                if self.selected_normalizer is not None:
                    dfs_to_use = plot_settings.get("dfs_to_use", [])
                    dfs_to_use = [x.replace("_", "_normalized_") for x in dfs_to_use]
                    plot_settings.update({"dfs_to_use": dfs_to_use})
                getattr(self, plot_name)(**plot_settings)
        self.logger.info("Done creating plots")

    def add_intensity_column(self, option_name: str, name_in_file: str, name_in_plot: str,
                             scale: str = "normal", df: Optional[pd.DataFrame] = None):
        """
        | Adds two options to all_intensities_dict and all_tree_dict, called option_name and option_name_log2.

        Parameters
        ----------
        option_name
            the name that the added data has internally, can be referred to via the df_to_use option e.g. *lfq* or
            *ibaq*
        name_in_file
            prefix of the columns e.g. *Intensity* or *LFQ intensity*
        name_in_plot
            shown name in the plots e.g. *LFQ Intensity* or "iBAQ values"
        scale
            is the data in "normal" or in "log2" scale
        df
            can be passed to use instead of :attr:`BasePlotter.intensity_df`

        """
        if df is None:
            if self.intensity_df is None:
                self.logger.warning("No intensity df provided")
                return
            df = self.intensity_df
        if not any((col.startswith(name_in_file) for col in df)):
            self.logger.warning("%s columns could not be found in data", name_in_file)
            return
        self.logger.debug("Adding option %s and %s_log2", option_name, option_name)
        if scale == "log2":
            option_name.replace("_log2", "")
        self.int_mapping.update({option_name: name_in_file, f"{option_name}_log2": name_in_file})
        self.intensity_label_names.update({option_name: name_in_plot,
                                           f"{option_name}_log2": rf"$Log_2$ {name_in_plot}"})

        # extract all raw intensities from the dataframe
        # replace all 0 with nan and remove the prefix from the columns
        intensities = df.loc[:, [c for c in df.columns if c.startswith(self.int_mapping[option_name])]
            ].replace({0: np.nan}).rename(lambda x: x.replace(self.int_mapping[option_name], ""), axis=1)
        if scale == "log2":
            intensities = np.exp2(intensities)
        # ensure data will not have faulty values after log2 transformation
        assert np.isinf(intensities).sum().sum() == 0
        assert (intensities < 1).sum().sum() == 0
        # filter all rows where all intensities are nan
        mask = (~intensities.isna()).sum(axis=1) != 0
        intensities = intensities[mask]

        # add log2 intensities
        intensities_log2 = np.log2(intensities)

        self.all_intensities_dict.update({
            option_name: intensities, f"{option_name}_log2": intensities_log2,
        })

        self.all_tree_dict.update({
            f"{option_name}": DataTree.from_analysis_design(
                self.analysis_design, intensities, self.configs.get("has_techrep", False)
            ),
            f"{option_name}_log2": DataTree.from_analysis_design(
                self.analysis_design, intensities_log2, self.configs.get("has_techrep", False)
            )
        })

    def add_normalized_option(self, df_to_use: str, normalizer: Union[Type[Normalization.BaseNormalizer], Any],
                              norm_option_name: str):
        """
        | Adds a new option/key of available data sets in all_intensities_dict and all_tree_dict by taking the data set
          all_intensities_dict[df_to_use], performing the normalization on the data and then adding the new option with
          :meth:`add_intensity_column`.

        Parameters
        ----------
        df_to_use
            data set that should be normalized
        normalizer
            normalizer either derived from :class:`~mspypeline.module.Normalization.BaseNormalizer` or a class with a
            :func:`fit_transform`
        norm_option_name
            suffix of the new option name
        """
        if df_to_use not in self.all_intensities_dict:
            self.logger.warning("normalization option %s could not be added", df_to_use)
            return
        df_to_use_no_log2 = df_to_use.replace("_log2", "")
        new_option_name = f"{df_to_use_no_log2}_{norm_option_name}"
        if new_option_name in self.all_tree_dict:
            self.logger.info("%s already exists as option", new_option_name)
            return
        import inspect
        if inspect.isclass(normalizer):
            normalizer = normalizer()
        if isinstance(normalizer, Normalization.BaseNormalizer):
            input_scale = "log2" if "log2" in df_to_use else "normal"
            setattr(normalizer, "input_scale", input_scale)
            setattr(normalizer, "output_scale", "normal")
            setattr(normalizer, "col_name_prefix", norm_option_name)
        assert hasattr(normalizer, "fit_transform"), "normalizer must have fit_transform method"
        data = self.all_intensities_dict[df_to_use].copy()
        data = normalizer.fit_transform(data)
        self.add_intensity_column(new_option_name, norm_option_name + " ",
                                  f"{norm_option_name.replace('_', ' ')} "
                                  f"{self.intensity_label_names[df_to_use_no_log2]}",
                                  scale="normal", df=data)

    def create_report(self):
        raise NotImplementedError

    def get_venn_group_data(self, df_to_use: str, level: int, non_na_function=get_number_of_non_na_values):
        """
        | Calculates which proteins can be compared between groups or are unique for a group of the selected level (see
          :ref:`thresholding`) and then counts these proteins per group.

        Parameters
        ----------
        df_to_use
            which dataframe/intensity should be analysed
        level
            at which level of the data tree should the data be compared
        non_na_function
            threshold function to determine if proteins can be compared, default: :func:`get_number_of_non_na_values`
        Returns
        -------
        Dict
            Dictionary containing the proteins that can be compared per group
        """
        n_children = {
            key: self.all_tree_dict[df_to_use][key].get_total_number_children()
            for key in self.all_tree_dict[df_to_use].level_keys_full_name[level]
        }
        counts_per_group = self.all_tree_dict[df_to_use].groupby(level, method=lambda x: (x > 0).sum())
        named_sets = {}
        for group, n_child in n_children.items():
            minimum = non_na_function(n_child)
            named_sets[group] = set(counts_per_group.index[counts_per_group[group] >= minimum])
        return named_sets

    def get_venn_data_per_key(self, df_to_use: str, key: str):
        """
        | Counts the protein intensity values greater than 0 (number of detected proteins) for each replicate of a group
          from the selected level.

        Parameters
        ----------
        df_to_use
            which dataframe/intensity should be analysed
        key
            which data node/group of samples should be compared

        Returns
        -------
        Dict
            Dictionary containing the proteins detected per sample
        """
        df = self.all_tree_dict[df_to_use].aggregate(key, method=None) > 0
        per_group_dict = {column: set(df.index[df[column]]) for column in df}
        return per_group_dict

    @add_end_docstrings(plot_para_return_docstring.format(
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_venn`, "
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_bar_venn`"
    ))
    @validate_input
    def plot_venn_groups(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
        """
        | Venn diagrams conduce the graphical illustration of set theory. In the ``mspypeline`` protein counts (greater
          than zero) constitute the sets and set relationships indicate the number of proteins that are shared between
          two or more sets. Thereby the similarity of detected proteins of a set can be assessed.
        | The method creates both a venn diagram and a bar-venn diagram comparing the similarity of the groups on
          the selected level (based on protein counts).
        | The ordinary venn diagram is quite intuitive, but it supports a maximum of three comparisons in the
          ``mspypeline``.
        | The bar-venn diagram holds the advantage of allowing an unlimited number of comparison sets. These figures
          consists of two combined graphs, an upper bar diagram, tha indicates the number of unique or shared proteins
          of a set or overlapping sets. The lower graph indicates which set or sets are being compared, respectively,
          which protein count (upper graph) belongs to which comparison (lower graph).

        | For overview of plots see :ref:`analysis options <detection-plots>`
        | For exemplary plot see :ref:`gallery <venn-group>`

        .. note::
            * A venn diagram can compare a maximum of 3 samples.
            * A bar-venn diagram can compare more than 3 samples.
            * If the selected level has more than 3 groups, only the bar-venn diagram will be created.
            * If the selected level has more than 6 groups no diagram will be created

        .. note::
            To determine which proteins can be compared between the groups and which are unique for one group an
            internal :ref:`threshold function <thresholding>` is applied.
        """
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use],
                                   ex=f"group_level_{level}", split_files=True,
                                   df_to_use=df_to_use, level=level, save_path=self.file_dir_venn)
                plot_kwargs.update(**kwargs)
                # create venn diagrams comparing all replicates within an experiment
                named_sets = self.get_venn_group_data(df_to_use, level)
                # save the resulting venn diagram
                plot = matplotlib_plots.save_venn(named_sets=named_sets, **plot_kwargs)
                plots.append(plot)
                # create a mixture of bar and venn diagram
                plot = matplotlib_plots.save_bar_venn(named_sets=named_sets, **plot_kwargs)
                plots.append(plot)
        return plots

    @add_end_docstrings(plot_para_return_docstring.format(
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_venn`, "
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_bar_venn`"
    ))
    @validate_input
    def plot_venn_results(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
        """
        | Venn diagrams conduce the graphical illustration of set theory. In the ``mspypeline`` protein counts (greater
          than zero) constitue the sets and set relationships indicate the number of proteins that are shared between
          two or more sets. Thereby the similarity of detected proteins of a set can be assessed.
        | The method creates both a venn diagram and a bar-venn diagram comparing the similarity of the replicates
          of each group from the selected level (based on protein counts).
        | The ordinary venn diagram is quite intuitive, but it supports a maximum of three comparisons in the
          ``mspypeline``.
        | The bar-venn diagram holds the advantage of allowing an unlimited number of comparison sets. These figures
          consists of two combined graphs, an upper bar diagram, tha indicates the number of unique or shared proteins
          of a set or overlapping sets. The lower graph indicates which set or sets are being compared, respectively,
          which protein count (upper graph) belongs to which comparison (lower graph).

        | For overview of plots see :ref:`analysis options <detection-plots>`
        | For exemplary plot see :ref:`gallery <venn-rep>`

        .. note::
            * A venn diagram can compare a maximum of 3 samples.
            * A bar-venn diagram can compare more than 3 samples.
            * If a group of the selected level has more than 3 replicates, only the bar-venn diagram will be created.
            * If the selected level has more than 6 groups no diagram will be created
        """
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                for key in self.all_tree_dict[df_to_use].level_keys_full_name[level]:
                    plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use], ex=key,
                                       df_to_use=df_to_use, level=level, save_path=self.file_dir_venn)
                    plot_kwargs.update(**kwargs)
                    named_sets = self.get_venn_data_per_key(df_to_use, key)
                    # save the resulting venn diagram
                    plot = matplotlib_plots.save_venn(named_sets, **plot_kwargs)
                    plots.append(plot)
                    # create a mixture of bar and venn diagram
                    plot = matplotlib_plots.save_bar_venn(named_sets, **plot_kwargs)
                    plots.append(plot)
        return plots

    def get_detection_counts_data(self, df_to_use: str, level: int, **kwargs) -> Dict[str, pd.DataFrame]:
        """
        | Counts the number of intensity values greater than 0 per protein (number of samples that the protein was
          detected in) per group of the selected level.

        Parameters
        ----------
        df_to_use
            which dataframes/intensities should be analysed
        level
            at which level of the data tree should the data be compared
        kwargs
            accepts kwargs
        Returns
        -------
        Dict
            Dictionary with key *"counts"* to a DataFrame containing the counts of proteins detected in a sample

        """
        level_values = self.all_tree_dict[df_to_use].level_keys_full_name[level]
        level_counts = []
        for full_name in level_values:
            intensities = self.all_tree_dict[df_to_use][full_name].aggregate(None)
            # from 0 to number of replicates, how often was each protein detected
            counts = (intensities > 0).sum(axis=1)
            counts = counts.value_counts().drop([0], errors="ignore").rename(full_name)
            level_counts.append(counts)
        level_counts = pd.concat(level_counts, axis=1).astype("Int64").sort_index()
        return {"counts": level_counts}

    @add_end_docstrings(plot_para_return_docstring.format(
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_detection_counts_results`"
    ))
    @validate_input
    def plot_detection_counts(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
        """
        | Bar diagram showing how often proteins were detected in a number of replicates for each group.

        | For overview of plots see :ref:`analysis options <detection-plots>`
        | For exemplary plot see :ref:`gallery <detection-counts>`

        """
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                data = self.get_detection_counts_data(df_to_use=df_to_use, level=level, **kwargs)
                if data:
                    plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use],
                                       df_to_use=df_to_use, level=level, save_path=self.file_dir_descriptive)
                    plot_kwargs.update(**kwargs)
                    plot = matplotlib_plots.save_detection_counts_results(**data, **plot_kwargs)
                    plots.append(plot)
        return plots

    def get_detected_proteins_per_replicate_data(self, df_to_use: str, level: int, **kwargs
                                                 ) -> Dict[str, Dict[str, pd.Series]]:
        """
        | Counts the number of protein intensity values greater than 0 (number of detected proteins) per sample of a
          group from the selected level.

        Parameters
        ----------
        df_to_use
            which dataframes/intensities should be analysed
        level
            at which level of the data tree should the data be compared
        kwargs
            accepts kwargs
        Returns
        -------
        Dict
            Dictionary with key *"all_height"* to a mapping of protein counts as Series per group
        """
        # determine number of rows and columns in the plot based on the number of experiments
        level_values = self.all_tree_dict[df_to_use].level_keys_full_name[level]
        all_heights = {}
        for experiment in level_values:
            intensities = self.all_tree_dict[df_to_use][experiment].aggregate(None)
            counts = (intensities > 0).sum(axis=1)
            counts = counts[counts > 0]
            heights = [len(counts)]
            for col in intensities:
                h = len(intensities[col][intensities[col] > 0])
                heights.append(h)
            all_heights[experiment] = pd.Series(heights, index=["Total"] + intensities.columns.to_list(),
                                                name=experiment)
        return {"all_heights": all_heights}

    @add_end_docstrings(plot_para_return_docstring.format(
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_detected_proteins_per_replicate_results`"
    ))
    @validate_input
    def plot_detected_proteins_per_replicate(self, dfs_to_use: Union[str, Iterable[str]],
                                             levels: Union[int, Iterable[int]], **kwargs):
        """
        | Bar diagram showing the number of detected proteins per sample (protein intensities greater than zero) as well
          as the total number of detected proteins for each group of a selected level.
        | The average number of detected proteins per group is indicated as gray
          dashed line.

        | For overview of plots see :ref:`analysis options <detection-plots>`
        | For exemplary plot see :ref:`gallery <detected-proteins>`

        """
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                data = self.get_detected_proteins_per_replicate_data(df_to_use=df_to_use, level=level, **kwargs)
                if data:
                    plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use],
                                       df_to_use=df_to_use, level=level, save_path=self.file_dir_descriptive)
                    plot_kwargs.update(**kwargs)
                    plot = matplotlib_plots.save_detected_proteins_per_replicate_results(**data, **plot_kwargs)
                    plots.append(plot)
        return plots

    def get_intensity_histograms_data(self, df_to_use: str, level: int, **kwargs):
        """
        | Get protein intensity values for each sample per group of the selected level.

        Parameters
        ----------
        df_to_use
            which dataframes/intensities should be analysed
        level
            at which level of the data tree should the data be compared
        kwargs
            accepts kwargs
        Returns
        -------
        Dict
            Dictionary with key *"hist_data"* to a DataFrame containing the protein intensity values per group
        """
        return {"hist_data": self.all_tree_dict[df_to_use].groupby(level, method=None)}

    @add_end_docstrings(plot_para_return_docstring.format(
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_intensity_histogram_results`"
    ))
    @validate_input
    def plot_intensity_histograms(self, dfs_to_use: Union[str, Iterable[str]],
                                  levels: Union[int, Iterable[int]], **kwargs):
        """
        | For each group of the selected level a histogram is created that counts the occurrence of the binned intensity
          values of each sample.
        | If *"show_mean"* is set to True in the :ref:`configs <default-yaml>` the mean intensity of the plotted samples
          of a group will be shown as gray dashed line.

        | For overview of plots see :ref:`analysis options <detection-plots>`
        | For exemplary plot see :ref:`gallery <int-hist>`
        """
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                data = self.get_intensity_histograms_data(df_to_use=df_to_use, level=level, **kwargs)
                if data:
                    plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use],
                                       df_to_use=df_to_use, level=level, save_path=self.file_dir_descriptive)
                    plot_kwargs.update(**kwargs)
                    plot = matplotlib_plots.save_intensity_histogram_results(**data, **plot_kwargs)
                    plots.append(plot)
        return plots

    def get_scatter_replicates_data(self, df_to_use: str, full_name: str) -> Dict[str, pd.DataFrame]:
        """
        | Get protein intensity values for each sample of a selected group.

        Parameters
        ----------
        df_to_use
            which dataframes/intensities should be analysed
        full_name
            which data node/group of samples should be compared
        Returns
        -------
        Dict
            Dictionary with key *"scatter_data"* to a DataFrame containing the protein intensity values per replicate
        """
        data = self.all_tree_dict[df_to_use][full_name].aggregate(None)
        if data.empty:
            return {}
        return {"scatter_data": data}

    @add_end_docstrings(plot_para_return_docstring.format(
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_scatter_replicates_results`"
    ))
    @validate_input
    def plot_scatter_replicates(self, dfs_to_use: Union[str, Iterable[str]],
                                levels: Union[int, Iterable[int]], **kwargs):
        """
        | For all replicates per group of the selected level, pairwise comparisons of the protein intensities are
          plotted and their correlation, calculated with the the Pearson’s correlation coefficient r2, is indicated.
        | Unique proteins per replicate are shown at the bottom and right side of the graph (substitution of na values
          by min value of data set).
        | For a group with more than 2 replicates, each pairwise comparison of the replicates is calculated and plotted
          together in one graph. For every group of the selected level one plot will be created.

        | For overview of plots see :ref:`analysis options <detection-plots>`
        | For exemplary plot see :ref:`gallery <scatter-rep>`
        """
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                for full_name in self.all_tree_dict[df_to_use].level_keys_full_name[level]:
                    data = self.get_scatter_replicates_data(df_to_use=df_to_use, full_name=full_name)
                    if data:
                        plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use], full_name=full_name,
                                           df_to_use=df_to_use, level=level, save_path=self.file_dir_descriptive)
                        plot_kwargs.update(**kwargs)
                        plot = matplotlib_plots.save_scatter_replicates_results(**data, **plot_kwargs)
                        plots.append(plot)
        return plots

    def get_rank_data(self, df_to_use: str, full_name: str, **kwargs) -> Dict[str, pd.Series]:
        """
        | Get protein intensity values of the selected group and rank the proteins by their intensity value.

        Parameters
        ----------
        df_to_use
            which dataframes/intensities should be analysed
        full_name
            which data node/group of samples should be compared
        Returns
        -------
        Dict
            Dictionary with key *"rank_data"* to Series containing the protein intensities of the group ranked by
            intensity value
        """
        return {"rank_data": self.all_tree_dict[df_to_use][full_name].aggregate().sort_values(ascending=False)}

    @add_end_docstrings(plot_para_return_docstring.format(
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_rank_results`"
    ))
    @validate_input
    def plot_rank(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
        """
        | In the rank plot all proteins are sorted by intensity value and plotted against their rank. For every group
          of the selected level one plot is created, averaging the protein intensities of the replicates of a group.
        | The highest intensity accounts for rank 0, the lowest intensity for the number of proteins - 1 whereby
          proteins with missing values are neglected. The median intensity of all proteins is given in the legend.
        | :ref:`Pathway analysis protein lists <pathway-proteins>` can be applied to the rank plot to provide
          information about the median intensity or rank of pathways of interest. If a protein is part of a selected
          pathway it will be presented in color and the median rank of all proteins of a given pathway is indicated.
          Multiple pathways can be selected and will be represented in the same graph as distinct groups.

        | For overview of plots see :ref:`analysis options <detection-plots>`
        | For exemplary plot see :ref:`gallery <rank>`
        """
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                for level_key in self.all_tree_dict[df_to_use].level_keys_full_name[level]:
                    data = self.get_rank_data(df_to_use=df_to_use, full_name=level_key, **kwargs)
                    if data:
                        plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use], full_name=level_key,
                                           df_to_use=df_to_use, level=level, save_path=self.file_dir_descriptive,
                                           interesting_proteins=self.interesting_proteins)
                        plot_kwargs.update(**kwargs)
                        plot = matplotlib_plots.save_rank_results(**data, **plot_kwargs)
                        plots.append(plot)
        return plots

    def get_relative_std_data(self, df_to_use: str, full_name: str, **kwargs) -> Dict[str, pd.DataFrame]:
        """
        | Calculate which proteins of a group can be used for the analysis (see :ref:`thresholding`) and filters proteins
          below the threshold out.

        Parameters
        ----------
        df_to_use
            which dataframes/intensities should be analysed
        full_name
            which data node/group of samples should be compared
        Returns
        -------
        Dict
            Dictionary with key *"intensities"* to a DataFrame containing the protein intensities of the group

        """
        intensities = self.all_tree_dict[df_to_use][full_name].aggregate(None)
        non_na = get_number_of_non_na_values(intensities.shape[1])
        mask = (intensities > 0).sum(axis=1) >= non_na
        intensities = intensities[mask]
        if intensities.empty:
            self.logger.warning("data for %s is empty", full_name)
            return {}
        return {"intensities": intensities}

    @add_end_docstrings(plot_para_return_docstring.format(
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_relative_std_results`"
    ))
    @validate_input
    def plot_relative_std(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
        """
        | Illustrates the relative standard deviation of the proteins between samples of a group which can help to
          understand how much fluctuation of the measured intensities is present between the replicates. Low deviation
          indicates that measured intensities are stable over multiple samples.
        | For each group of the selected level one plot will be created.

        | For overview of plots see :ref:`analysis options <detection-plots>`
        | For exemplary plot see :ref:`gallery <rel-std>`

        .. note::
            To determine which proteins can be compared between the two samples an internal :ref:`threshold function
            <thresholding>` is applied.
        """

        plots = []
        # TODO check with log2 thresholds
        for level in levels:
            for df_to_use in dfs_to_use:
                for full_name in self.all_tree_dict[df_to_use].level_keys_full_name[level]:
                    data = self.get_relative_std_data(df_to_use=df_to_use, full_name=full_name, **kwargs)
                    if data:
                        plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use],
                                           experiment_name=full_name, df_to_use=df_to_use, level=level,
                                           save_path=self.file_dir_descriptive)
                        plot_kwargs.update(**kwargs)
                        plot = matplotlib_plots.save_relative_std_results(**data, **plot_kwargs)
                        plots.append(plot)
        return plots

    def get_pathway_analysis_data(self, df_to_use: str, level: int, pathway: str, equal_var=True, **kwargs):
        """
        | Filters out all proteins of the given pathways for all samples per group of the selected level, then
          calculates the pairwise significances between the groups with an independent t-test (see
          :meth:`plot_pathway_analysis`) for all those proteins that can be compared (see :ref:`thresholding`).

        Parameters
        ----------
        df_to_use
            which dataframes/intensities should be analysed
        level
            at which level of the data tree should the data be compared
        pathway
            which pathway should be analysed
        equal_var
            should equal variance be assumed
        kwargs
            accepts kwargs
        Returns
        -------
        Dict
            Dictionary with keys *"protein_intensities"* to a DataFrame containing the protein intensities of detected
            proteins from all given pathways per group and *"significances"* to a DataFrame containing the calculated
            significances between groups for each protein of all given pathways
        """
        level_keys = self.all_tree_dict[df_to_use].level_keys_full_name[level]
        found_proteins = set(self.interesting_proteins[pathway])
        found_proteins &= set(self.all_intensities_dict[df_to_use].index)
        found_proteins = list(found_proteins)
        if len(found_proteins) < 1:
            self.logger.warning("Skipping pathway %s in pathway analysis because no proteins were found", pathway)
            return {}
        protein_intensities = self.all_tree_dict[df_to_use].groupby(level, method=None, index=found_proteins).\
            sort_index(0).sort_index(1, ascending=False)
        significances = []
        for protein in protein_intensities.index:
            per_protein_significant = []
            for e1, e2 in combinations(level_keys, 2):
                v1 = protein_intensities.loc[protein].loc[e1]
                v2 = protein_intensities.loc[protein].loc[e2]
                # filter entries with too many nans based on function
                non_na_group_1 = get_number_of_non_na_values(v1.shape[0])
                non_na_group_2 = get_number_of_non_na_values(v2.shape[0])
                mask_1 = (v1 > 0).sum(axis=0) >= non_na_group_1
                mask_2 = (v2 > 0).sum(axis=0) >= non_na_group_2
                mask = np.logical_and(mask_1, mask_2)
                if not mask:
                    per_protein_significant.append(np.nan)
                else:
                    test = stats.ttest_ind(v1, v2, equal_var=equal_var, nan_policy="omit")
                    per_protein_significant.append(test[1])
            significances.append(per_protein_significant)
        significances = pd.DataFrame(significances, index=protein_intensities.index,
                                     columns=pd.MultiIndex.from_tuples([(e1, e2) for e1, e2 in
                                                                        combinations(level_keys, 2)]))
        return {"protein_intensities": protein_intensities, "significances": significances}

    @add_end_docstrings(plot_para_return_docstring.format(
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_pathway_analysis_results`"
    ))
    @validate_input
    def plot_pathway_analysis(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
        """
        | In the pathway analysis, for each protein of a desired :ref:`pathway <pathway-proteins>` a subplot is created
          displaying the intensities of the protein for all groups of the selected level.
        | Additionally, significances are calculated for each pairwise comparison between groups with an independent
          `t-test <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html>`__. For every
          selected pathway, two figures are created, one displaying the significances and the other not displaying them.
        | For a group of multiple samples, the protein intensity is plotted for each sample (single scatter dot) which
          are jointly presented in uniform coloring.

        | For overview of plots see :ref:`analysis options <statistic-plots>`
        | For exemplary plot see :ref:`gallery <pathway-analysis>`

        .. note::
            To determine which proteins can be compared between two groups an internal :ref:`threshold function
            <thresholding>` is applied.
        """
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                for pathway in list(self.interesting_proteins.keys()):
                    data = self.get_pathway_analysis_data(level=level, df_to_use=df_to_use, pathway=pathway, **kwargs)
                    if data:
                        plot_kwargs = dict(pathway=pathway, save_path=self.file_dir_pathway, df_to_use=df_to_use,
                                           level=level, intensity_label=self.intensity_label_names[df_to_use])
                        plot_kwargs.update(**kwargs)
                        plot = matplotlib_plots.save_pathway_analysis_results(**data, **plot_kwargs)
                        plots.append(plot)
        return plots

    def get_pathway_timeline_data(self):
        pass

    def plot_pathway_timecourse(self, df_to_use: str = "raw", show_suptitle: bool = False,
                                levels: Iterable = (2,), **kwargs):
        """
        not yet implemented
        """
        """
        group_colors = {
            "SD": "#808080",
            "4W": "#0b8040",
            "6W": "#ec2024",
            "8W": "#4378bb"
        }
        for level in levels:
            level_keys = self.all_tree_dict[df_to_use].level_keys_full_name[level]
            groups = {k: "SD" if "SD" in k else k.split("_")[1] for k in level_keys}
            if len(set(groups.values()) - set(group_colors.keys())):
                self.logger.warning("Skipping pathway timeline plot because of incorrect match between found groups and 
                target groups")
                continue
            x_values = {}
            # for key in level_keys:
            #     sample_names = self.all_tree_dict[df_to_use][key].aggregate(None).columns
            #     x_values.update({sample: sum([int(s.replace("W", "")) for s in sample.split("_") if s.endswith("W")]) 
            for sample in sample_names})
            for key in level_keys:
                x_values.update({key: sum([int(s.replace("W", "")) for s in key.split("_") if s.endswith("W")])})
            max_time = max(x_values.values())
            for pathway in self.interesting_proteins:
                plt.close("all")
                found_proteins = set(self.interesting_proteins[pathway])
                found_proteins &= set(self.all_intensities_dict[df_to_use].index)
                found_proteins = sorted(list(found_proteins))
                if len(found_proteins) < 1:
                    self.logger.warning("Skipping pathway %s in pathway timeline because no proteins were found", 
                    pathway)
                    continue
        """

    def get_experiment_comparison_data(self, df_to_use: str, full_name1: str, full_name2: str):
        """
        | Gets protein intensities for all samples of a given group, then calculates the proteins that can be compared
          between groups and those that are unique for each group (see :ref:`thresholding`) and takes the mean intensity
          of these proteins.

        Parameters
        ----------
        df_to_use
            which dataframes/intensities should be analysed
        full_name1
            name of the first data node/group that should be compared to 'full_name2'
        full_name2
            name of the second data node/group that should be compared to 'full_name1'

        Returns
        -------
        Dict
            Dictionary with keys *"protein_intensities_sample1"* and *"protein_intensities_sample2"* to Series
            containing the mean protein intensities of sample 1 and sample 2  and *"exclusive_sample1"* and
            *"exclusive_sample2"* to Series containing the mean intensities of unique proteins for sample 1 and
            sample 2.
        """
        protein_intensities_sample1 = self.all_tree_dict[df_to_use][full_name1].aggregate(None)
        protein_intensities_sample2 = self.all_tree_dict[df_to_use][full_name2].aggregate(None)
        mask, exclusive_1, exclusive_2 = get_intersection_and_unique(protein_intensities_sample1,
                                                                     protein_intensities_sample2)
        # flatten all replicates
        exclusive_sample1 = protein_intensities_sample1[exclusive_1].mean(axis=1)
        exclusive_sample2 = protein_intensities_sample2[exclusive_2].mean(axis=1)
        protein_intensities_sample1 = protein_intensities_sample1[mask].mean(axis=1)
        protein_intensities_sample2 = protein_intensities_sample2[mask].mean(axis=1)
        if protein_intensities_sample1.empty and protein_intensities_sample2.empty:
            self.logger.warning("protein samples of %s and %s are both empty", full_name1, full_name2)
            return {}
        return {
            "protein_intensities_sample1": protein_intensities_sample1,
            "protein_intensities_sample2": protein_intensities_sample2,
            "exclusive_sample1": exclusive_sample1, "exclusive_sample2": exclusive_sample2
        }

    @add_end_docstrings(plot_para_return_docstring.format(
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_experiment_comparison_results`"
    ))
    @validate_input
    def plot_experiment_comparison(self, dfs_to_use: Union[str, Iterable[str]],
                                   levels: Union[int, Iterable[int]], **kwargs):
        """
        | For all groups of the selected level, pairwise comparisons of the protein intensities are
          plotted and their correlation, calculated with the the Pearson’s correlation coefficient r2, is indicated.
        | Unique proteins per group are shown at the bottom and right side of the graph (substitution of na values
          by min value of data set).
        | For every pairwise comparison of the groups from the selected level, one plot will be created.

        | For overview of plots see :ref:`analysis options <detection-plots>`
        | For exemplary plot see :ref:`gallery <scatter-group>`

        .. note::
            To determine which proteins can be compared between the two groups and which are unique for one group an
            internal :ref:`threshold function <thresholding>` is applied.
        """
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                for ex1, ex2 in combinations(self.all_tree_dict[df_to_use].level_keys_full_name[level], 2):
                    data = self.get_experiment_comparison_data(df_to_use=df_to_use, full_name1=ex1, full_name2=ex2)
                    if data:
                        plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use], sample1=ex1,
                                           sample2=ex2, df_to_use=df_to_use, level=level,
                                           save_path=self.file_dir_descriptive)
                        plot_kwargs.update(**kwargs)
                        plot = matplotlib_plots.save_experiment_comparison_results(**data, **plot_kwargs)
                        plots.append(plot)
        return plots

    def get_go_analysis_data(self, df_to_use: str, level: int):
        """
        | Calculates an enrichment analysis for all samples per group of the selected level and for each given GO list
          (see :meth:`plot_go_analysis`). Significances are calculated with a fisher exact test.

        Parameters
        ----------
        df_to_use
            which dataframes/intensities should be analysed
        level
            at which level of the data tree should the data be compared
        Returns
        -------
        Dict
            Dictionary with keys *"heights"* to a ddict containing the counts of proteins per sample of each given
            GO list, *"test_results"* to a ddict containing the corresponding Fisher's exact test results and
            *"go_length"* to a list containing the total number of proteins of each chosen GO list

        """
        if not self.go_analysis_gene_names:
            return {}
        background = set(self.all_intensities_dict[df_to_use].index)
        heights = ddict(list)
        test_results = ddict(list)
        go_length = []
        for compartiment, all_pathway_genes in self.go_analysis_gene_names.items():
            all_pathway_genes = set(all_pathway_genes)
            # get all pathway genes that were detected throughout all experiments
            pathway_genes = all_pathway_genes & background
            # all genes that are not in the pathway
            not_pathway_genes = background - pathway_genes
            # sanity check
            assert pathway_genes | not_pathway_genes == background
            heights["Total"].append(len(pathway_genes))
            go_length.append(len(all_pathway_genes))
            for experiment in self.all_tree_dict[df_to_use].level_keys_full_name[level]:
                # create df with intensity means for specific experiment over all replicates
                mean_intensity = self.all_tree_dict[df_to_use][experiment].aggregate()
                # get all proteins with mean intensity > 0
                experiment_genes = set(mean_intensity[mean_intensity > 0].index)
                # all other experiments are not detected
                not_experiment_genes = background - experiment_genes
                # sanity check
                assert experiment_genes | not_experiment_genes == background
                # append height for the plot
                heights[experiment].append(len(experiment_genes & pathway_genes))
                table = pd.DataFrame([
                    [len(experiment_genes & pathway_genes), len(experiment_genes & not_pathway_genes)],
                    [len(not_experiment_genes & pathway_genes), len(not_experiment_genes & not_pathway_genes)]
                ], columns=["in pathway", "not in pathway"], index=["in experiment", "not in experiment"])
                oddsratio, pvalue = stats.fisher_exact(table, alternative='greater')
                # chi2, p, dof, ex = stats.chi2_contingency(table, correction=True)
                # self.logger.debug(f"{chi2}, {dof}, {ex}")
                # self.logger.debug(f"{compartiment} {experiment}: table: {table} fisher: {pvalue:.4f}, chi2: {p:.4f}")
                test_results[experiment].append(pvalue)
        return {"heights": heights, "test_results": test_results, "go_length": go_length}

    @add_end_docstrings(plot_para_return_docstring.format(
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_go_analysis_results`"
    ))
    @validate_input
    def plot_go_analysis(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
        """
        | In the GO analysis, an enrichment analysis is performed for each selected
          :ref:`GO Term file <go-term-proteins>` (based on protein counts). The number of detected proteins from a GO
          term found in each group of the analysis design is illustrated as the length of the corresponding bar. P
          values shown at the end of a bar indicate the calculated significance. Samples referred to as "Total"
          represent the complete data set and numbers at the top of the graph accord to the count of detected proteins
          in all samples over the total number of proteins in the GO term.
        | For p-value calculateion, first, for each GO term, a list *"pathway_genes"* is created by taking the
          intersection of the proteins from the GO list and the total detected proteins.
        | Secondly, a list of *"non_pathway_genes"* is created which comprises total detected proteins but proteins in
          *"pathway_genes"*.
        | Third, a list of *"experiment_genes"* and *"non_experiment_genes"* is created in a similar fashion where an
          experiment references to a sample/group of samples of the data set.
        | Lastly, a one-tailed `fisher exact test
          <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html>`__ is calculated to
          retrieve statistical significances based on the following contingency table:

        +------------------------+--------------------------------------+------------------------------------------+
        |                        | in pathway                           | not in pathway                           |
        +========================+======================================+==========================================+
        | **in experiment**      | experiment_genes & pathway_genes     | experiment_genes & not_pathway_genes     |
        +------------------------+--------------------------------------+------------------------------------------+
        | **not in experiment**  | not_experiment_genes & pathway_genes | not_experiment_genes & not_pathway_genes |
        +------------------------+--------------------------------------+------------------------------------------+

        | The resulting p-value is thus, also dependent on the overall protein count of the sample/group of samples.

        | For overview of plots see :ref:`analysis options <statistic-plots>`
        | For exemplary plot see :ref:`gallery <go-analysis>`
        """
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                data = self.get_go_analysis_data(df_to_use=df_to_use, level=level)
                if data:
                    plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use],
                                       go_analysis_gene_names=self.go_analysis_gene_names,
                                       df_to_use=df_to_use, level=level, save_path=self.file_dir_go_analysis)
                    plot_kwargs.update(**kwargs)
                    plot = matplotlib_plots.save_go_analysis_results(**data, **plot_kwargs)
                    plots.append(plot)
        return plots

    def get_r_volcano_data(self, g1: str, g2: str, df_to_use: str):
        """
        | Gets the protein intensities for all samples of the two given groups, then calculates the proteins that can be
          compared between groups and those unique for each group (see :ref:`thresholding`).
        | Hands over the protein intensities to be compared to the R package ``limma`` that outputs the logFC, p-value,
          adjusted p value and other data which is calculated based on a `moderated t-statistic
          <https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf>`__.
        | Results are converted back to python format afterwards.

        .. note::
            This function uses the R package limma which will be automatically downloaded the first time this analysis
            is performed

        Parameters
        ----------
        g1
            first sample that should be analysed (downregulated)
        g2
            second sample that should be analysed (upregulated)
        df_to_use
            which dataframes/intensities should be analysed

        Returns
        -------
        Dict
            Dictionary with keys *"volcano_data"* to a DataFrame containing processed output of the ``limma.eBayes``
            analysis, *"unique_g1"* and *"unique_g2"* to Series containing the unique protein intensities per group

        """
        # install r packages
        from mspypeline.helpers.Utils import install_r_dependencies
        r_package_names = ("BiocManager", )
        r_bioconducter_package_names = ("limma", )
        install_r_dependencies(r_package_names, r_bioconducter_package_names)

        # import r interface package
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.packages import importr
        from rpy2.robjects.conversion import localconverter
        # allow conversion of pd objects to r
        pandas2ri.activate()

        # import r packages
        limma = importr("limma")
        # get data of the samples
        v1 = self.all_tree_dict[df_to_use][g1].aggregate(None)
        v2 = self.all_tree_dict[df_to_use][g2].aggregate(None)
        if v2.shape[1] < 2 or v2.shape[1] < 2:
            self.logger.warning("Skipping Volcano plot for comparison: %s, %s because the groups contain only "
                                "%s and %s experiments", g1, g2, v1.shape[1], v2.shape[1])
            return {}
        mask, exclusive_1, exclusive_2 = get_intersection_and_unique(v1, v2)

        df = pd.concat([v1[mask], v2[mask]], axis=1)
        design = pd.DataFrame([[0] * v1.shape[1] + [1] * v2.shape[1],
                               [1] * v1.shape[1] + [0] * v2.shape[1]], index=[g2, g1]).T

        if df.empty:
            return {}
        # transform to r objects
        with localconverter(ro.default_converter + pandas2ri.converter):
            r_df = ro.conversion.py2rpy(df)
            r_design = ro.conversion.py2rpy(design)
        # run the r code
        fit = limma.lmFit(r_df, r_design)
        c_matrix = make_contrasts(g1, g2)
        contrast_fit = limma.contrasts_fit(fit, c_matrix)
        fit_bayes = limma.eBayes(contrast_fit, trend=True)
        res = limma.topTable(fit_bayes, adjust="BH", number=df.shape[0])
        # transform back to python
        with localconverter(ro.default_converter + pandas2ri.converter):
            # positive is upregulated in v2 / g2
            ress = ro.conversion.rpy2py(res)
        # possible names are keys of this dict
        plot_data = ress.loc[:, ["logFC", "AveExpr", "P.Value", "adj.P.Val"]]
        plot_data = plot_data.rename({"P.Value": "pval", "adj.P.Val": "adjpval"}, axis=1)
        # calculate mean intensity for unique genes
        unique_g1 = v1[exclusive_1].mean(axis=1).rename(f"{df_to_use} mean intensity")
        unique_g2 = v2[exclusive_2].mean(axis=1).rename(f"{df_to_use} mean intensity")

        return {"volcano_data": plot_data, "unique_g1": unique_g1, "unique_g2": unique_g2}

    @validate_input
    def plot_r_volcano(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]],
                       sample1: str = None, sample2: str = None, **kwargs):
        """
        | A volcano plot illustrates the statistical inferences from a pairwise comparison of the two groups.
        | The plot shows the log2 fold change between two different conditions against the -log10(p-value)
          (based on protein intensities). The p-value is determined using the R limma package (`moderated t-statistic
          <https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf>`__).
        | Dashed lines indicate the fold change cutoff (default = log2(2) and p-value cutoff (default = p < 0.05) by
          which proteins are considered significant (blue and red) or non significant (gray). Measured intensities of
          unique proteins are indicated at the sides of the volcano plot for each groups (light blue and orange).
        | Volcano plots also permit the annotation of mapped proteins. This can be achieved by labeling a number of
          the most significant proteins for each group or by selecting a
          :ref:`pathway analysis protein list <pathway-proteins>`.
        | For every pairwise comparison of the groups of the selected level two volcano plots are created, where one
          plot has a set of proteins annotated and the other does not.

        | For overview of plots see :ref:`analysis options <statistic-plots>`
        | For exemplary plot see :ref:`gallery <volcano>`

        .. note::
           * should be used with log2 intensities
           * minimum of 3 samples per group required

        .. note::
            To determine which proteins can be compared between the two groups and which are unique for one group an
            internal :ref:`threshold function <thresholding>` is applied.

        Parameters
        ----------
        dfs_to_use
            which dataframes/intensities should be plotted
        levels
            at which level of the data tree should the data be compared
        sample1
            first sample that should be compared (downregulated)
        sample2
            second sample that should be compared (upregulated)
        kwargs
            accepts kwargs
        Returns
        -------
        List
            A list of all created plots.
        """
        plots = []
        plot_once = False
        for level in levels:
            for df_to_use in dfs_to_use:
                if sample1 and sample2:
                    level_keys = [sample1, sample2]
                    plot_once = True
                else:
                    level_keys = self.all_tree_dict[df_to_use].level_keys_full_name[level]
                for g1, g2 in combinations(level_keys, 2):
                    data = self.get_r_volcano_data(g1, g2, df_to_use)
                    if data:
                        plot_kwargs = dict(g1=g1, g2=g2, save_path=self.file_dir_volcano, df_to_use=df_to_use,
                                           intensity_label=self.intensity_label_names[df_to_use],
                                           interesting_proteins=self.interesting_proteins, split_files=True)
                        plot_kwargs.update(**kwargs)
                        plot = matplotlib_plots.save_volcano_results(**data, **plot_kwargs)
                        plots.append(plot)
            if plot_once:
                break
        return plots

    def get_pca_data(self, df_to_use: str, level: int, n_components: int = 2, fill_value: float = 0,
                     no_missing_values: bool = True, fill_na_before_norm: bool = False, **kwargs):
        """
        | Gets protein intensities for all samples per group processes data according to given arguments and then
          performs a dimensionality reduction (PCA) using ``sklearn.decomposition.PCA``.

        Parameters
        ----------
        df_to_use
            which dataframes/intensities should be analysed
        level
            at which level of the data tree should the data be compared
        n_components
            how many principal components should be calculated
        fill_value
            if data should be interpolated, which fill value should be used
        no_missing_values
            should missing values be neglected
        fill_na_before_norm
            if data should be interpolated, should this be done before normalisation
        kwargs
            accepts kwargs
        Returns
        -------
        Dict
            Dictionary with keys *"pca_data"* to a DataFrame containing the output of a PCA using `
            ``sklearn.decomposition`` and *"pca_fit"* to a PCA object that was fitted to normalized input data
        """
        data_input = self.all_tree_dict[df_to_use].groupby(level, method=None)
        if no_missing_values:
            data_input = data_input.dropna(axis=0)
        else:
            if fill_na_before_norm:
                data_input.fillna(fill_value, inplace=True)
        data_norm = data_input.subtract(data_input.mean(axis=1), axis=0).divide(data_input.std(axis=1), axis=0)
        if not fill_na_before_norm:
            data_norm.fillna(fill_value, inplace=True)
        data_transform = data_norm.T

        pca: PCA = PCA(n_components=n_components).fit(data_transform)
        df = pd.DataFrame(pca.transform(data_transform).T, columns=data_input.columns,
                          index=[f"PC_{i}" for i in range(1, n_components + 1)])

        return {"pca_data": df, "pca_fit": pca}

    @add_end_docstrings(plot_para_return_docstring.format(
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_pca_results`"
    ))
    @validate_input
    def plot_pca_overview(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
        """
        | With the option to perform PCA, data can be studied for its variance and in doing so, parameters can be
          determined that have most strongly affected the variability between samples. The created PCA compares all
          components against each other, the default is set to 2 components where only PC 1 and PC 2 are compared.
        | PCA results do not change in dependence on the chosen level, however, determining the level on which the
          data should be compared influences the coloring of the scatter elements. Each group of the selected level is
          colored differently.
        | Multiple different analysis options can be chosen to generate a PCA (see: :ref:`multiple option config
          <default-yaml>`).

        | For overview of plots see :ref:`analysis options <detection-plots>`
        | For exemplary plot see :ref:`gallery <pca>`
       """
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                data = self.get_pca_data(level=level, df_to_use=df_to_use, **kwargs)
                if data:
                    plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use],
                                       df_to_use=df_to_use, level=level, save_path=self.file_dir_descriptive)
                    plot_kwargs.update(**kwargs)
                    plot = matplotlib_plots.save_pca_results(**data, **plot_kwargs)
                    plots.append(plot)
        return plots

    def get_boxplot_data(self, df_to_use: str, level: int, **kwargs) -> dict:
        """
        | Get protein intensities for all samples per group of the selected level and then sorts samples by their median
          intensity.

        Parameters
        ----------
        df_to_use
            which dataframes/intensities should be analysed
        level
            at which level of the data tree should the data be compared
        kwargs
            accepts kwargs
        Returns
        -------
        Dict
            Dictionary with key *"protein_intensities"* to a DataFrame containing the protein intensities per group
            sorted by median intensity
        """
        df = self.all_tree_dict[df_to_use].groupby(level)
        # sort columns by median intensity
        df = df[df.median().sort_values(ascending=True).index]
        return {"protein_intensities": df}

    @add_end_docstrings(plot_para_return_docstring.format(
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_boxplot_results`"
    ))
    @validate_input
    def plot_boxplot(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
        """
        | A standard boxplot displaying the five quantile distribution per group of the selected level and ranking the
          groups by median intensity from the bottom of the graph to the top.
        | The boxplot is part of the :ref:`Normalization overview <norm-overview>`.

        | For overview of plots see :ref:`analysis options <add-python-plots>`
        | For exemplary plot see :ref:`gallery <boxplot>`

        """
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                data = self.get_boxplot_data(df_to_use=df_to_use, level=level, **kwargs)
                if data:
                    plot_kwargs = dict(level=level, df_to_use=df_to_use, save_path=self.file_dir_descriptive,
                                       intensity_label=self.intensity_label_names[df_to_use])
                    plot_kwargs.update(**kwargs)
                    plot = matplotlib_plots.save_boxplot_results(**data, **plot_kwargs)
                    plots.append(plot)
        return plots

    def get_n_protein_vs_quantile_data(self, df_to_use: str, level: int, quantile_range: Optional[np.array] = None,
                                       **kwargs):
        """
        | Gets protein intensities for all samples per group, counts the number of intensity values greater than 0
          (total number of detected proteins) and the quantiles per sample.

        Parameters
        ----------
        df_to_use
            which dataframes/intensities should be analysed
        level
            at which level of the data tree should the data be compared
        quantile_range
            which quantile range should be used for analysis
        kwargs
            accepts kwargs
        Returns
        -------
        Dict
            Dictionary with keys *"quantiles"* to a DataFrame of calculated quantiles per sample and *"n_proteins"* to a
            Series of total number of identified proteins per sample
        """
        if quantile_range is None:
            quantile_range = np.arange(0.05, 1, 0.05)
        df = self.all_tree_dict[df_to_use].groupby(level)
        n_proteins = (df > 0).sum(axis=0)
        quantiles = df.quantile(quantile_range)
        return {"quantiles": quantiles, "n_proteins": n_proteins}

    @add_end_docstrings(plot_para_return_docstring.format(
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_n_proteins_vs_quantile_results`"
    ))
    @validate_input
    def plot_n_proteins_vs_quantile(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]],
                                    **kwargs):
        """
        | Plots the quantile protein intensities against the number of identified proteins per sample.
        | Samples are indicated as a horizontal line of scatter dots where the color anf x position of a dot indicate
          the intensity value of the respective quantile. The y position of the dots of a sample point to the total
          number of detected proteins in that sample.
        | Solid, rather vertical lines indicate a linear fit of each quantile for all the samples.
        | This plot is part of the :ref:`Normalization overview <norm-overview>`.

        | For overview of plots see :ref:`analysis options <add-python-plots>`
        | For exemplary plot see :ref:`gallery <proteins-vs-quantiles>`

        """
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                data = self.get_n_protein_vs_quantile_data(df_to_use=df_to_use, level=level, **kwargs)
                if data:
                    plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use],
                                       df_to_use=df_to_use, level=level, save_path=self.file_dir_descriptive)
                    plot_kwargs.update(**kwargs)
                    plot = matplotlib_plots.save_n_proteins_vs_quantile_results(**data, **plot_kwargs)
                    plots.append(plot)
        return plots

    def get_kde_data(self, df_to_use: str, level: int, **kwargs) -> Dict[str, pd.DataFrame]:
        """
        | Gets the protein intensities for all samples per group of the selected level.

        Parameters
        ----------
        df_to_use
            which dataframes/intensities should be analysed
        level
            at which level of the data tree should the data be compared
        kwargs
            accepts kwargs
        Returns
        -------
        Dict
            Dictionary with key *"intensities"* to a DataFrame containing the protein intensities per group
        """
        intensities = self.all_tree_dict[df_to_use].groupby(level)
        return {"intensities": intensities}

    @add_end_docstrings(plot_para_return_docstring.format(
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_kde_results`"
    ))
    @validate_input
    def plot_kde(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
        """
        | In the kernel density estimate (KDE) plot, one density graph per sample is plotted indicating the intensity on
          the x axis and the density on the y axis. These plots should be presented on a log2 scale.
        | The KDE is well suited to study the influence of different :ref:`normalization methods <hyperparameter>` and
          :ref:`protein intensities <hyperparameter>` on the data which is why it is part if the
          :ref:`Normalization overview <norm-overview>`.

        | For overview of plots see :ref:`analysis options <add-python-plots>`
        | For exemplary plot see :ref:`gallery <kde>`
        """
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                data = self.get_kde_data(df_to_use=df_to_use, level=level, **kwargs)
                if data:
                    plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use],
                                       df_to_use=df_to_use, level=level, save_path=self.file_dir_descriptive)
                    plot_kwargs.update(**kwargs)
                    plot = matplotlib_plots.save_kde_results(**data, **plot_kwargs)
                    plots.append(plot)
        return plots

    @add_end_docstrings(plot_para_return_docstring.format(
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_normalization_overview_results`"
    ))
    @validate_input
    def plot_normalization_overview(
            self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs
    ) -> List[Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes, plt.Axes, plt.Axes]]]:
        """
        | The Normalization overview offers the opportunity to examine different aspects of the data in three distinct
          plots. For each :ref:`normalization method <hyperparameter>` provided an additional page will be attached to
          the resulting pdf file starting with the raw or not normalized data. That way it is possible to get a better
          understanding of the effects of the normalization methods on the data, to inspect the different approaches and
          to find the best suitable normalization for the data.
        | The normalization overview combines the plots :meth:`~mspypeline.BasePlotter.plot_kde`,
          :meth:`~mspypeline.BasePlotter.plot_n_proteins_vs_quantile` and :meth:`~mspypeline.BasePlotter.plot_boxplot`.

        | For overview of plots see :ref:`analysis options <norm-plots>`
        | For exemplary plot see :ref:`gallery <norm-overview>`
        """
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                n_prot_data = self.get_n_protein_vs_quantile_data(df_to_use=df_to_use, level=level, **kwargs)
                kde_data = self.get_kde_data(df_to_use=df_to_use, level=level, **kwargs)
                boxplot_data = self.get_boxplot_data(df_to_use=df_to_use, level=level, **kwargs)

                if n_prot_data and kde_data and boxplot_data:
                    plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use],
                                       df_to_use=df_to_use, level=level, save_path=self.file_dir_descriptive)
                    plot_kwargs.update(**kwargs)
                    plot = matplotlib_plots.save_normalization_overview_results(
                        **n_prot_data, **kde_data, **boxplot_data, **plot_kwargs
                    )
                    plots.append(plot)
        return plots

    def get_intensity_heatmap_data(self, df_to_use: str, level: int, sort_index: bool = False,
                                   sort_index_by_missing: bool = True, sort_columns_by_missing: bool = True, **kwargs):
        """
        | Get the protein intensities for all samples per group of the selected level and sorts samples and proteins
          according to settings.

        Parameters
        ----------
        df_to_use
            which dataframes/intensities should be analysed
        level
            at which level of the data tree should the data be compared
        sort_index
            should proteins be sorted alphanumerically
        sort_index_by_missing
            should proteins be sorted by number of missing values across samples
        sort_columns_by_missing
            should samples be sorted by number of missing values
        kwargs
            accepts kwargs
        Returns
        -------
        Dict
            Dictionary with key *"intensities"* to a DataFrame containing protein intensities of samples
        """
        intensities = self.all_tree_dict[df_to_use].groupby(level)
        if sort_index_by_missing:
            index = intensities.isna().sum(axis=1).sort_values().index
        elif sort_index:
            index = intensities.index.sort_values()
        else:
            index = intensities.index
        if sort_columns_by_missing:
            columns = intensities.isna().sum(axis=0).sort_values().index
        else:
            columns = intensities.columns

        return {"intensities": intensities.loc[index, columns]}

    @add_end_docstrings(plot_para_return_docstring.format(
        ":func:`~mspypeline.plotting_backend.matplotlib_plots.save_intensities_heatmap_result`"
    ))
    @validate_input
    def plot_intensity_heatmap(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]],
                               **kwargs):
        """
        | The intensity heatmap demonstrates protein intensities, where samples are given in rows on the y axis and
          proteins on the x axis. Missing values are colored in gray.
        | The heatmap can be used to spot patterns in the different :ref:`normalization methods <hyperparameter>` and to
          understand how different :ref:`intensity types <hyperparameter>` affect the data.
        | The :ref:`Heatmap overview <heatmap-overview>` is created from a series of intensity heatmap plots.

        | For overview of plots see :ref:`analysis options <add-python-plots>`
        | For exemplary plot see :ref:`gallery <int-heatmap>`
        """
        plots = []
        for df_to_use in dfs_to_use:
            for level in levels:
                data = self.get_intensity_heatmap_data(df_to_use=df_to_use, level=level, **kwargs)
                if data:
                    plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use],
                                       df_to_use=df_to_use, level=level, save_path=self.file_dir_normalization)
                    plot_kwargs.update(**kwargs)
                    plot = matplotlib_plots.save_intensities_heatmap_result(**data, **plot_kwargs)
                    plots.append(plot)
        return plots

    def plot_all_normalizer_overview(self, dfs_to_use, levels, plot_function, file_name, normalizers=None, **kwargs):
        """
        | Helper method to create a multipaged file containing one plot per normalization option.
        
        | For overview of plots see :ref:`analysis options <norm-plots>`
        | For exemplary plot see :ref:`gallery <norm-plots-gallery>`

        Parameters
        ----------
        dfs_to_use
            which dataframes/intensities should be plotted
        levels
            at which level of the data tree should the data be compared
        plot_function
            which plot should be created
        file_name
            name of the file that will be crated and saved
        normalizers
            normalizers either derived from :class:`~mspypeline.module.Normalization.BaseNormalizer` or a class with a
            :func:`fit_transform`
        kwargs
            accepts kwargs
        Returns
        -------
        list
            A list of all created plots
        """
        max_depth = dict_depth(self.analysis_design)
        if self.configs.get("has_techrep", False):
            max_depth -= 1
        plots = []
        normalizers = normalizers if normalizers is not None else {}
        normalizers.update(self.normalizers)
        plot_kwargs = dict()
        plot_kwargs.update(**kwargs)
        save_path = plot_kwargs.pop("save_path", self.file_dir_normalization)
        for df_to_use in dfs_to_use:
            for normaliser_name, normalizer in normalizers.items():
                self.add_normalized_option(df_to_use, normalizer, normaliser_name)
            dfs = [x for x in self.all_tree_dict if x.startswith(df_to_use.replace("_log2", ""))]
            if "log2" in df_to_use:
                dfs = [x for x in dfs if x.endswith("log2")]
            plot_kwargs["save_path"] = None  # make sure the plots dont get saved
            df_plots = plot_function(dfs, max_depth - 1, **plot_kwargs)
            plots += df_plots
            if save_path is not None:
                plot_kwargs.pop("save_path")  # use the other save path instead here
                save_path, result_name = matplotlib_plots.get_path_and_name_from_kwargs(
                    file_name, **plot_kwargs, df_to_use=df_to_use, save_path=save_path)
                matplotlib_plots.collect_plots_to_pdf(os.path.join(save_path, result_name), *df_plots)
        return plots

    @validate_input
    def plot_normalization_overview_all_normalizers(self, dfs_to_use, levels, **kwargs):
        """
        | Will create the :meth:`plot_normalization_overview` for all normalization methods.

        | For overview of plots see :ref:`analysis options <norm-plots>`
        | For exemplary plot see :ref:`gallery <norm-overview>`

        Parameters
        ----------
        dfs_to_use
            which dataframes/intensities should be plotted
        levels
            at which level of the data tree should the data be compared
        kwargs
            accepts kwargs
        """
        plot_kwargs = dict()
        plot_kwargs.update(kwargs)
        return self.plot_all_normalizer_overview(
            dfs_to_use=dfs_to_use, levels=levels, plot_function=self.plot_normalization_overview,
            file_name="normalization_overview_all_normalizers", **plot_kwargs
        )

    @validate_input
    def plot_heatmap_overview_all_normalizers(self, dfs_to_use, levels, **kwargs):
        """
        | Will create the :meth:`plot_intensity_heatmap` for all normalization methods.
        | The intensity heatmap demonstrates protein intensities, where samples are given in rows on the y axis and
          proteins on the x axis. Missing values are colored in gray.
        | The heatmap can be used to spot patterns in the different :ref:`normalization methods <hyperparameter>` and to
          understand how different :ref:`intensity types <hyperparameter>` affect the data.

        | For overview of plots see :ref:`analysis options <norm-plots>`
        | For exemplary plot see :ref:`gallery <heatmap-overview>`

        Parameters
        ----------
        dfs_to_use
            which dataframes/intensities should be plotted
        levels
            at which level of the data tree should the data be compared
        kwargs
            accepts kwargs
        """
        plot_kwargs = dict(sort_index=False, sort_index_by_missing=True, sort_columns_by_missing=False)
        plot_kwargs.update(kwargs)
        return self.plot_all_normalizer_overview(
            dfs_to_use=dfs_to_use, levels=levels, plot_function=self.plot_intensity_heatmap,
            file_name="heatmap_overview_all_normalizers", **plot_kwargs
        )
