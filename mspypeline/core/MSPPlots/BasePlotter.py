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
from typing import Dict, Type, Iterable, Optional, Union, Any
from sklearn.decomposition import PCA
from copy import deepcopy

from mspypeline.core import MSPInitializer
from mspypeline.file_reader import BaseReader
from mspypeline.plotting_backend import matplotlib_plots
from mspypeline.modules import default_normalizers, Normalization, DataTree
from mspypeline.helpers import get_number_rows_cols_for_fig, get_number_of_non_na_values, \
    get_intersection_and_unique, get_logger, dict_depth

# TODO VALIDATE descriptive plots not changing between log2 and non log2

FIG_FORMAT = ".pdf"


def validate_input(f):
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


class BasePlotter:
    possible_plots = [
        "plot_detection_counts", "plot_number_of_detected_proteins", "plot_intensity_histograms",
        "plot_relative_std", "plot_rank", "plot_pathway_analysis", "plot_pathway_timeline",
        "plot_scatter_replicates", "plot_experiment_comparison", "plot_go_analysis", "plot_venn_results",
        "plot_venn_groups", "plot_r_volcano", "plot_pca_overview", "plot_boxplot",
        "plot_normalization_overview_all_normalizers", "plot_heatmap_overview_all_normalizers"
    ]

    def __init__(
            self,
            start_dir: str,
            reader_data: dict,
            intensity_df_name: str = "",
            interesting_proteins: Optional[Dict[str, pd.Series]] = None,
            go_analysis_gene_names: Optional[Dict[str, pd.Series]] = None,
            configs: Optional[dict] = None,
            required_reader: Optional[str] = None,
            intensity_entries=(),
            loglevel=logging.DEBUG
    ):
        self.logger = get_logger(self.__class__.__name__, loglevel=loglevel)
        # general information
        # TODO make start dir optional
        self.start_dir = start_dir
        # read in optional arguments
        self.configs = {} if configs is None else configs
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
        # create file structure and folders
        # path for venn diagrams
        self.file_dir_venn = os.path.join(self.start_dir, "venn")
        # path for descriptive plots
        self.file_dir_descriptive = os.path.join(self.start_dir, "descriptive")
        # path for pathway analysis
        self.file_dir_pathway = os.path.join(self.start_dir, "pathway_analysis")
        # path for go analysis
        self.file_dir_go_analysis = os.path.join(self.start_dir, "go_analysis")
        # path for volcano plots
        self.file_dir_volcano = os.path.join(self.start_dir, "volcano")

    @classmethod
    def from_MSPInitializer(cls, mspinti_instance: MSPInitializer, **kwargs):
        default_kwargs = dict(
            start_dir=mspinti_instance.start_dir,
            reader_data=mspinti_instance.reader_data,
            intensity_df_name="",
            interesting_proteins=mspinti_instance.interesting_proteins,
            go_analysis_gene_names=mspinti_instance.go_analysis_gene_names,
            configs=mspinti_instance.configs,
            required_reader=None,
            intensity_entries=(),
            loglevel=mspinti_instance.logger.getEffectiveLevel()
        )
        default_kwargs.update(**kwargs)
        return cls(**default_kwargs)

    @classmethod
    def from_file_reader(cls, reader_instance: BaseReader, **kwargs):
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
        global_settings = self.configs.get("global_settings", {})
        self.logger.debug(f"got global settings: %s", global_settings)
        for plot_name in self.possible_plots:
            plot_settings_name = plot_name + "_settings"
            plot_settings = self.configs.get(plot_settings_name, {})
            plot_settings.update({k: v for k, v in global_settings.items() if k not in plot_settings})
            if plot_settings.pop("create_plot", False):
                self.logger.debug(f"creating plot {plot_name}")
                getattr(self, plot_name)(**plot_settings)
        self.logger.info("Done creating plots")

    def add_intensity_column(self, option_name: str, name_in_file: str, name_in_plot: str,
                             scale: str = "normal", df: Optional[pd.DataFrame] = None):
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
        self.intensity_label_names.update({option_name: name_in_plot, f"{option_name}_log2": rf"$Log_2$ {name_in_plot}"})

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
                self.analysis_design, intensities, self.configs.get("has_replicates", False)
            ),
            f"{option_name}_log2": DataTree.from_analysis_design(
                self.analysis_design, intensities_log2, self.configs.get("has_replicates", False)
            )
        })

    def add_normalized_option(self, df_to_use: str, normalizer: Union[Type[Normalization.BaseNormalizer], Any], norm_option_name: str):
        df_to_use_no_log2 = df_to_use.replace("_log2", "")
        new_option_name = f"{df_to_use_no_log2}_{norm_option_name}"
        if new_option_name in self.all_tree_dict:
            self.logger.info("%s already exists as option")
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
                                  f"{norm_option_name.replace('_', ' ')} {self.intensity_label_names[df_to_use_no_log2]}",
                                  scale="normal", df=data)

    def create_report(self):
        raise NotImplementedError

    def get_venn_group_data(self, df_to_use: str, level: int, non_na_function=get_number_of_non_na_values):
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

    def get_venn_data_per_key(self, df_to_use: str = "raw", key: str = None):
        df = self.all_tree_dict[df_to_use].aggregate(key, method=None) > 0
        per_group_dict = {column: set(df.index[df[column]]) for column in df}
        return per_group_dict

    @validate_input
    def plot_venn_groups(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use],
                                   ex=f"group_level_{level}",
                                   df_to_use=df_to_use, level=level, save_path=self.file_dir_venn)
                plot_kwargs.update(**kwargs)
                # create venn diagrams comparing all replicates within an experiment
                named_sets = self.get_venn_group_data(df_to_use, level)
                # save the resulting venn diagram
                plot = matplotlib_plots.save_venn(named_sets, **plot_kwargs)
                plots.append(plot)
                # create a mixture of bar and venn diagram
                plot = matplotlib_plots.save_bar_venn(named_sets, **plot_kwargs)
                plots.append(plot)
        return plots

    @validate_input
    def plot_venn_results(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
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
        Counts the number values greater than 0 per protein

        Parameters
        ----------
        df_to_use
            DataFrame to use
        level
            level to use
        kwargs
            accepts kwargs

        Returns
        -------
        Dictionary with DataFrame containing the counts

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

    @validate_input
    def plot_detection_counts(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
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

    def get_number_of_detected_proteins_data(self, df_to_use: str, level: int, **kwargs) -> Dict[str, Dict[str, pd.Series]]:
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
            all_heights[experiment] = pd.Series(heights, index=["Total"] + intensities.columns.to_list(), name=experiment)
        return {"all_heights": all_heights}

    @validate_input
    def plot_number_of_detected_proteins(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                data = self.get_number_of_detected_proteins_data(df_to_use=df_to_use, level=level, **kwargs)
                if data:
                    plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use],
                                       df_to_use=df_to_use, level=level, save_path=self.file_dir_descriptive)
                    plot_kwargs.update(**kwargs)
                    plot = matplotlib_plots.save_number_of_detected_proteins_results(**data, **plot_kwargs)
                    plots.append(plot)
        return plots

    def get_intensity_histograms_data(self, df_to_use: str, level: int, **kwargs):
        return {"hist_data": self.all_tree_dict[df_to_use].groupby(level, method=None)}

    @validate_input
    def plot_intensity_histograms(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
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
        data = self.all_tree_dict[df_to_use][full_name].aggregate(None)
        if data.empty:
            return {}
        return {"scatter_data": data}

    @validate_input
    def plot_scatter_replicates(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
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
        return {"rank_data": self.all_tree_dict[df_to_use][full_name].aggregate().sort_values(ascending=False)}

    @validate_input
    def plot_rank(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
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

        Parameters
        ----------
        df_to_use
            which dataframe/intensity should be used
        full_name
            name of the experiment which should be accepted
        kwargs
            accepts kwargs

        Returns
        -------
        Dictionary with the intensity of the experiment and the relative standard deviations

        """
        intensities = self.all_tree_dict[df_to_use][full_name].aggregate(None)
        non_na = get_number_of_non_na_values(intensities.shape[1])
        mask = (intensities > 0).sum(axis=1) >= non_na
        intensities = intensities[mask]
        if intensities.empty:
            self.logger.warning("data for %s is empty", full_name)
            return {}
        return {"intensities": intensities}

    @validate_input
    def plot_relative_std(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
        plots = []
        # TODO check with log2 thresholds
        for level in levels:
            for df_to_use in dfs_to_use:
                for full_name in self.all_tree_dict[df_to_use].level_keys_full_name[level]:
                    data = self.get_relative_std_data(df_to_use=df_to_use, full_name=full_name, **kwargs)
                    if data:
                        plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use], experiment_name=full_name,
                                           df_to_use=df_to_use, level=level, save_path=self.file_dir_descriptive)
                        plot_kwargs.update(**kwargs)
                        plot = matplotlib_plots.save_relative_std_results(**data, **plot_kwargs)
                        plots.append(plot)
        return plots

    def get_pathway_analysis_data(self, df_to_use: str, level: int, pathway: str, equal_var=True, **kwargs):
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
                                     columns=pd.MultiIndex.from_tuples([(e1, e2) for e1, e2 in combinations(level_keys, 2)]))
        return {"protein_intensities": protein_intensities, "significances": significances}

    @validate_input
    def plot_pathway_analysis(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
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

    def plot_pathway_timeline(self, df_to_use: str = "raw", show_suptitle: bool = False, levels: Iterable = (2,), **kwargs):
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
                self.logger.warning("Skipping pathway timeline plot because of incorrect match between found groups and target groups")
                continue
            x_values = {}
            # for key in level_keys:
            #     sample_names = self.all_tree_dict[df_to_use][key].aggregate(None).columns
            #     x_values.update({sample: sum([int(s.replace("W", "")) for s in sample.split("_") if s.endswith("W")]) for sample in sample_names})
            for key in level_keys:
                x_values.update({key: sum([int(s.replace("W", "")) for s in key.split("_") if s.endswith("W")])})
            max_time = max(x_values.values())
            for pathway in self.interesting_proteins:
                plt.close("all")
                found_proteins = set(self.interesting_proteins[pathway])
                found_proteins &= set(self.all_intensities_dict[df_to_use].index)
                found_proteins = sorted(list(found_proteins))
                if len(found_proteins) < 1:
                    self.logger.warning("Skipping pathway %s in pathway timeline because no proteins were found", pathway)
                    continue
                n_rows, n_cols = get_number_rows_cols_for_fig(found_proteins)
                fig, axarr = plt.subplots(n_rows, n_cols, figsize=(n_cols * int(max_time / 5), 4 * n_rows))
                if show_suptitle:
                    fig.suptitle(pathway)
                try:
                    axiterator = axarr.flat
                except AttributeError:
                    axiterator = [axarr]
                protein_minimum = self.all_intensities_dict[df_to_use].max().max()
                protein_maximum = self.all_intensities_dict[df_to_use].min().min()
                for protein, ax in zip(found_proteins, axiterator):
                    ax.set_title(protein)
                    ax.set_xlabel(f"Age [weeks]")
                    ax.set_ylabel(f"{self.intensity_label_names[df_to_use]}")
                    for idx, experiment in enumerate(level_keys):
                        protein_intensities = self.all_tree_dict[df_to_use][experiment].aggregate(None, index=protein)
                        mask = protein_intensities > 0
                        protein_minimum = min(protein_minimum, protein_intensities[mask].min())
                        protein_maximum = max(protein_maximum, protein_intensities[mask].max())
                        ax.scatter([x_values[experiment]] * sum(mask), protein_intensities[mask],
                                   label=f"{groups[experiment]}", color=group_colors[groups[experiment]])
                # adjust labels based on overall min and max of the pathway
                try:
                    axiterator = axarr.flat
                except AttributeError:
                    axiterator = [axarr]
                for protein, ax in zip(found_proteins, axiterator):
                    ax.set_ylim(bottom=protein_minimum * 0.99, top=protein_maximum * 1.01)
                    ax.set_xlim(left=0, right=max_time + 1)
                handles, labels = axiterator[0].get_legend_handles_labels()
                fig.legend(handles, labels, bbox_to_anchor=(1.04, 0.5), loc="center left")
                fig.tight_layout(rect=[0, 0.03, 1, 0.95])
                res_path = os.path.join(self.file_dir_pathway, f"pathway_timeline_{pathway}" + FIG_FORMAT)
                fig.savefig(res_path, dpi=200, bbox_inches="tight")

    def get_experiment_comparison_data(self, df_to_use: str, full_name1: str, full_name2: str):
        protein_intensities_sample1 = self.all_tree_dict[df_to_use][full_name1].aggregate(None)
        protein_intensities_sample2 = self.all_tree_dict[df_to_use][full_name2].aggregate(None)
        mask, exclusive_1, exclusive_2 = get_intersection_and_unique(protein_intensities_sample1, protein_intensities_sample2)
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

    @validate_input
    def plot_experiment_comparison(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
        # TODO correlation of log2 and not log 2 data is different
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                for ex1, ex2 in combinations(self.all_tree_dict[df_to_use].level_keys_full_name[level], 2):
                    data = self.get_experiment_comparison_data(df_to_use=df_to_use, full_name1=ex1, full_name2=ex2)
                    if data:
                        plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use], sample1=ex1, sample2=ex2,
                                           df_to_use=df_to_use, level=level, save_path=self.file_dir_descriptive)
                        plot_kwargs.update(**kwargs)
                        plot = matplotlib_plots.save_experiment_comparison_results(**data, **plot_kwargs)
                        plots.append(plot)
        return plots

    def get_go_analysis_data(self, df_to_use: str, level: int):
        if not self.go_analysis_gene_names:
            return {}
        background = set(self.all_intensities_dict[df_to_use].index)
        heights = ddict(list)
        test_results = ddict(list)
        for compartiment, all_pathway_genes in self.go_analysis_gene_names.items():
            all_pathway_genes = set(all_pathway_genes)
            # get all pathway genes that were detected throughout all experiments
            pathway_genes = all_pathway_genes & background
            # all genes that are not in the pathway
            not_pathway_genes = background - pathway_genes
            # sanity check
            assert pathway_genes | not_pathway_genes == background
            heights["background"].append(len(pathway_genes))
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
        return {"heights": heights, "test_results": test_results}

    @validate_input
    def plot_go_analysis(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                data = self.get_go_analysis_data(df_to_use=df_to_use, level=level)
                if data:
                    plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use],
                                       go_analysis_gene_names=self.go_analysis_gene_names,
                                       df_to_use=df_to_use, level=level, save_path=self.file_dir_descriptive)
                    plot_kwargs.update(**kwargs)
                    plot = matplotlib_plots.save_go_analysis_results(**data, **plot_kwargs)
                    plots.append(plot)
        return plots

    def get_r_volcano_data(self, g1: str, g2: str, df_to_use: str, level: int):
        # import r interface package
        from rpy2.robjects.packages import importr
        from rpy2.robjects import pandas2ri
        # allow conversion of pd objects to r
        pandas2ri.activate()

        # install r packages
        from mspypeline.helpers.Utils import install_r_dependencies
        r_package_names = ("BiocManager", )
        r_bioconducter_package_names = ("limma", )
        install_r_dependencies(r_package_names, r_bioconducter_package_names)

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
        r_df = pandas2ri.py2ri(df)
        r_design = pandas2ri.py2ri(design)
        r_rownames = pandas2ri.py2ri(df.index)
        # add the index to the dataframe because it will be sorted
        r_df.rownames = r_rownames
        # run the r code
        fit = limma.lmFit(r_df, r_design)
        # TODO replace this function with one that handles special characters
        c_matrix = limma.makeContrasts(f"{g2}-{g1}", levels=r_design)
        contrast_fit = limma.contrasts_fit(fit, c_matrix)
        fit_bayes = limma.eBayes(contrast_fit)
        res = limma.topTable(fit_bayes, adjust="BH", number=df.shape[0])
        # transform back to python
        with warnings.catch_warnings():
            # catch a warning from ri2py where DataFrame.from_items is being used
            warnings.simplefilter("ignore", FutureWarning)
            # positive is upregulated in v2 / g2
            ress = pandas2ri.ri2py(res)
        # extract index
        ress.index = pd.Index([x for x in res.rownames], name="Gene_names")
        # possible names are keys of this dict
        plot_data = ress.loc[:, ["logFC", "AveExpr", "P.Value", "adj.P.Val"]]
        plot_data = plot_data.rename({"P.Value": "pval", "adj.P.Val": "adjpval"}, axis=1)
        # calculate mean intensity for unique genes
        unique_g1 = v1[exclusive_1].mean(axis=1).rename(f"{df_to_use} mean intensity")
        unique_g2 = v2[exclusive_2].mean(axis=1).rename(f"{df_to_use} mean intensity")

        return {"volcano_data": plot_data, "unique_g1": unique_g1, "unique_g2": unique_g2}

    @validate_input
    def plot_r_volcano(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
        # TODO both adj and un adj should be available
        plots = []
        for level in levels:
            for df_to_use in dfs_to_use:
                level_keys = self.all_tree_dict[df_to_use].level_keys_full_name[level]
                for g1, g2 in combinations(level_keys, 2):
                    data = self.get_r_volcano_data(g1, g2, df_to_use, level)
                    if data:
                        plot_kwargs = dict(g1=g1, g2=g2, save_path=self.file_dir_volcano, df_to_use=df_to_use, level=level,
                                           intensity_label=self.intensity_label_names[df_to_use])
                        plot_kwargs.update(**kwargs)
                        plot = matplotlib_plots.save_volcano_results(**data, **plot_kwargs)
                        plots.append(plot)
        return plots

    def get_pca_data(self, df_to_use: str, level: int, n_components: int = 4, fill_value: float = 0, fill_na_before_norm: bool = False, **kwargs):
        data_input = self.all_tree_dict[df_to_use].groupby(level, method=None)
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

    @validate_input
    def plot_pca_overview(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
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
        Generates data for the boxplot, where columns are sorted by median intensities

        Parameters
        ----------
        df_to_use
            which dataframe should be used
        level
            on which level should the data be grouped
        kwargs
            accepts kwargs

        Returns
        -------
        A dictionary with the the protein intensities
        """
        df = self.all_tree_dict[df_to_use].groupby(level)
        # sort columns by median intensity
        df = df[df.median().sort_values(ascending=True).index]
        return {"protein_intensities": df}

    @validate_input
    def plot_boxplot(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
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

    def get_n_protein_vs_quantile_data(self, df_to_use: str, level: int, quantile_range: Optional[np.array] = None, **kwargs):
        if quantile_range is None:
            quantile_range = np.arange(0.05, 1, 0.05)
        df = self.all_tree_dict[df_to_use].groupby(level)
        n_proteins = (df > 0).sum(axis=0)
        quantiles = df.quantile(quantile_range)
        return {"quantiles": quantiles, "n_proteins": n_proteins}

    @validate_input
    def plot_n_proteins_vs_quantile(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
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
        intensities = self.all_tree_dict[df_to_use].groupby(level)
        return {"intensities": intensities}

    @validate_input
    def plot_kde(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
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

    @validate_input
    def plot_normalization_overview(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
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

    @validate_input
    def plot_intensity_heatmap(self, dfs_to_use: Union[str, Iterable[str]], levels: Union[int, Iterable[int]], **kwargs):
        plots = []
        for df_to_use in dfs_to_use:
            for level in levels:
                data = self.get_intensity_heatmap_data(df_to_use=df_to_use, level=level, **kwargs)
                if data:
                    plot_kwargs = dict(intensity_label=self.intensity_label_names[df_to_use],
                                       df_to_use=df_to_use, level=level, save_path=self.file_dir_descriptive)
                    plot_kwargs.update(**kwargs)
                    plot = matplotlib_plots.save_intensities_heatmap_result(**data, **plot_kwargs)
                    plots.append(plot)
        return plots

    def plot_all_normalizer_overview(self, dfs_to_use, levels, plot_function, file_name, **kwargs):
        max_depth = dict_depth(self.analysis_design)
        if self.configs.get("has_replicates", False):
            max_depth -= 1
        plots = []
        normalizers = kwargs.get("normalizers", {})
        normalizers.update(self.normalizers)
        plot_kwargs = dict()
        plot_kwargs.update(**kwargs)
        plot_kwargs.update({"save_path": None})
        for df_to_use in dfs_to_use:
            for normaliser_name, normalizer in normalizers.items():
                self.add_normalized_option(df_to_use, normalizer, normaliser_name)
            dfs = [x for x in self.all_tree_dict if x.startswith(df_to_use.replace("_log2", ""))]
            if "log2" in df_to_use:
                dfs = [x for x in dfs if x.endswith("log2")]
            plots += plot_function(dfs, max_depth - 1, **plot_kwargs)
        matplotlib_plots.collect_plots_to_pdf(os.path.join(self.file_dir_descriptive, file_name), *plots)
        return plots

    @validate_input
    def plot_normalization_overview_all_normalizers(self, dfs_to_use, levels, **kwargs):
        plot_kwargs = dict()
        plot_kwargs.update(kwargs)
        return self.plot_all_normalizer_overview(
            dfs_to_use=dfs_to_use, levels=levels, plot_function=self.plot_normalization_overview,
            file_name="normalization_overview_all_normalizers.pdf", **plot_kwargs
        )

    @validate_input
    def plot_heatmap_overview_all_normalizers(self, dfs_to_use, levels, **kwargs):
        plot_kwargs = dict(sort_index=True, sort_index_by_missing=False, sort_columns_by_missing=False)
        plot_kwargs.update(kwargs)
        return self.plot_all_normalizer_overview(
            dfs_to_use=dfs_to_use, levels=levels, plot_function=self.plot_intensity_heatmap,
            file_name="heatmap_overview_all_normalizers.pdf", **plot_kwargs
        )
