import pandas as pd
import numpy as np
import os
import functools
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib_venn import venn3, venn2
import matplotlib.pyplot as plt
from scipy import stats
from itertools import combinations
from collections import defaultdict as ddict
from _collections_abc import Iterable
import logging
import warnings
from typing import Dict
from sklearn.decomposition import PCA

from mspypeline import MSPInitializer, matplotlib_plots
from mspypeline.helpers import get_number_rows_cols_for_fig, venn_names, get_number_of_non_na_values, plot_annotate_line,\
    get_intersection_and_unique, DataTree, get_logger

# TODO VALIDATE descriptive plots not changing between log2 and non log2

# plt.style.use('ggplot')

# TODO move these to the yml file
FIG_FORMAT = ".pdf"
# Venn diagram settings
# TODO figsize
VENN_TITLE_FONT_SIZE = 20
VENN_SET_LABEL_FONT_SIZE = 16
VENN_SUBSET_LABEL_FONT_SIZE = 14
# descriptive plots settings


class MSPPlots:
    possible_plots = [
        "plot_detection_counts", "plot_number_of_detected_proteins", "plot_intensity_histograms",
        "plot_relative_std", "plot_rank", "plot_pathway_analysis", "plot_pathway_timeline",
        "plot_scatter_replicates", "plot_experiment_comparison", "plot_go_analysis", "plot_venn_results",
        "plot_venn_groups", "plot_r_volcano", "plot_pca_overview", "plot_boxplot"
    ]

    def exception_handler(f):
        @functools.wraps(f)
        def wrapper(self, *args, **kwargs):
            try:
                ret = f(self, *args, **kwargs)
                return ret
            except PermissionError:
                self.logger.warning("Permission error in function %s. Did you forget to close the file?",
                                    str(f).split(" ")[1])
                return
        return wrapper

    def __init__(
            self,
            start_dir: str,
            reader_data: dict,
            intensity_df_name: str = "",
            interesting_proteins: dict = None,
            go_analysis_gene_names: dict = None,
            configs: dict = None,
            required_reader=None,
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

        if required_reader is not None:
            try:
                self.required_reader_data = reader_data[required_reader]
            except KeyError:
                self.logger.exception("Reader data does not provide information from %s reader", required_reader)
                raise

        self.intensity_df = None
        for reader in reader_data:
            for key in reader_data[reader]:
                if key == intensity_df_name:
                    self.intensity_df = reader_data[reader][key]
                    break

        # setup everything for all_intensity dict
        self.int_mapping = {}
        self.intensity_label_names = {}
        self.all_intensities_dict: Dict[str, pd.DataFrame] = {}
        for option_name, name_in_file, name_in_plot in intensity_entries:
            self.add_intensity_column(option_name, name_in_file, name_in_plot)

        # create the data trees
        tech_rep = self.configs.get("has_replicates", False)
        analysis_design = self.configs.get("analysis_design", {})
        if analysis_design:
            self.all_tree_dict: Dict[str, DataTree] = {
                intensity: DataTree.from_analysis_design(analysis_design, data, tech_rep)
                for intensity, data in self.all_intensities_dict.items()
            }
        else:
            self.all_tree_dict = {
                intensity: {}
                for intensity, data in self.all_intensities_dict.items()
            }
            self.logger.warning("No analysis design was provided. Most plotting functions will not work")

        # set all result dirs
        # create file structure and folders
        # TODO: for now just create all of them
        # TODO: dont create them on init
        # path for venn diagrams
        self.file_dir_venn = os.path.join(self.start_dir, "venn")
        os.makedirs(self.file_dir_venn, exist_ok=True)
        # path for descriptive plots
        self.file_dir_descriptive = os.path.join(self.start_dir, "descriptive")
        os.makedirs(self.file_dir_descriptive, exist_ok=True)
        # path for pathway analysis
        self.file_dir_pathway = os.path.join(self.start_dir, "pathway_analysis")
        os.makedirs(self.file_dir_pathway, exist_ok=True)
        # path for go analysis
        self.file_dir_go_analysis = os.path.join(self.start_dir, "go_analysis")
        os.makedirs(self.file_dir_go_analysis, exist_ok=True)
        # path for volcano plots
        self.file_dir_volcano = os.path.join(self.start_dir, "volcano")
        os.makedirs(self.file_dir_volcano, exist_ok=True)

    @classmethod
    def from_MSPInitializer(cls, mspinti_instance: MSPInitializer, intensity_entries=()):
        return cls(
            start_dir=mspinti_instance.start_dir,
            reader_data=mspinti_instance.reader_data,
            intensity_df_name="",
            interesting_proteins=mspinti_instance.interesting_proteins,
            go_analysis_gene_names=mspinti_instance.go_analysis_gene_names,
            configs=mspinti_instance.configs,
            required_reader=None,
            intensity_entries=intensity_entries,
            loglevel=mspinti_instance.logger.getEffectiveLevel()
            )

    def create_results(self):
        for plot_name in self.possible_plots:
            plot_settings_name = plot_name + "_settings"
            plot_settings = self.configs.get(plot_settings_name, {})
            if plot_settings.pop("create_plot", False):
                self.logger.debug(f"creating plot {plot_name}")
                getattr(self, plot_name)(**plot_settings)
        self.logger.info("Done creating plots")

    def add_intensity_column(self, option_name, name_in_file, name_in_plot):
        if not any((col.startswith(name_in_file) for col in self.intensity_df)):
            self.logger.warning("%s columns could not be found in data", name_in_file)
            return
        self.int_mapping.update({option_name: name_in_file, f"{option_name}_log2": name_in_file})
        self.intensity_label_names.update({option_name: name_in_plot, f"{option_name}_log2": rf"$Log_2$ {name_in_plot}"})

        # extract all raw intensities from the dataframe
        # replace all 0 with nan and remove the prefix from the columns
        intensities = self.intensity_df.loc[:,
                        [c for c in self.intensity_df.columns if c.startswith(self.int_mapping[option_name])]
        ].replace({0: np.nan}).rename(lambda x: x.replace(self.int_mapping[option_name], ""), axis=1)
        # filter all rows where all intensities are nan
        mask = (~intensities.isna()).sum(axis=1) != 0
        intensities = intensities[mask]

        # add log2 intensities
        intensities_log2 = np.log2(intensities)

        self.all_intensities_dict.update({
            option_name: intensities, f"{option_name}_log2": intensities_log2,
        })

    def create_report(self):
        raise NotImplementedError

    @exception_handler
    def save_venn(self, ex: str, named_sets: Dict[str, set], show_suptitle: bool = True):
        # TODO legend with colors and numbers
        creates_figure = True
        plt.close("all")
        plt.figure(figsize=(14, 7))
        title = ex
        if show_suptitle:
            plt.title(title, fontsize=VENN_TITLE_FONT_SIZE)

        # create venn diagram based on size of set
        sets = named_sets.values()
        set_names = named_sets.keys()
        if len(sets) < 2:
            creates_figure = False
            self.logger.warning(f"Could not create venn diagram for {ex} because it has less than 2 replicates")
        elif len(sets) == 2:
            venn = venn2(subsets=sets, set_labels=set_names)
        elif len(sets) == 3:
            venn = venn3(subsets=sets, set_labels=set_names)
        else:
            self.logger.warning(f"Could not create venn diagram for {ex}"
                                f" because it has more than 3 replicates ({len(sets)})")
            creates_figure = False

        # if a figure was created, do some further configuration and save it
        if creates_figure:
            res_path = os.path.join(
                self.file_dir_venn, f"venn_replicate_{title.replace(' ', '_')}" + FIG_FORMAT
            )

            for text in venn.set_labels:
                try:
                    text.set_fontsize(VENN_SET_LABEL_FONT_SIZE)
                except AttributeError:
                    pass
            for text in venn.subset_labels:
                try:
                    text.set_fontsize(VENN_SUBSET_LABEL_FONT_SIZE)
                except AttributeError:
                    pass
            plt.savefig(res_path, dpi=200, bbox_inches="tight")
        plt.close("all")

    @exception_handler
    def save_bar_venn(self, ex: str, named_sets: Dict[str, set], show_suptitle: bool = True):
        plt.close("all")
        if len(named_sets) > 6:
            self.logger.warning("Skipping bar-venn for %s because it has more than 6 experiments", ex)
            return

        # create a mapping from name to a y coordinate
        y_mappings = {name: i for i, name in enumerate(named_sets)}
        # get all the heights and other info required for the plot
        heights = []
        x = []
        ys = []
        for i, (intersected, unioned, result) in enumerate(venn_names(named_sets)):
            heights.append(len(result))
            x.append(i)
            ys.append([y_mappings[x] for x in intersected])

        # initial figure setup
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(1 * len(heights), 7))
        name = ex
        if show_suptitle:
            fig.suptitle(name, fontsize=20)
        # create the bar plot
        ax1.bar(x, heights, color="skyblue")
        # add text to the bar plot
        for x_level, height in zip(x, heights):
            ax1.text(x_level, max(heights) / 2, height, verticalalignment='center', horizontalalignment='center')
        ax1.set_ylabel("Number of proteins")

        # create the line plots
        for x_level, y in zip(x, ys):
            # we just want to draw a straight line every time so we repeat x as often as needed
            ax2.plot([x_level] * len(y), y, linestyle="-", color="gray", marker=".")
        # replace the yticks with the names of the samples
        ax2.set_yticks([i for i in range(len(y_mappings))])
        ax2.set_yticklabels(y_mappings)
        ax2.set_ylabel("Sample name")
        ax2.set_xlabel("Number of comparison")

        # save the result with the according name
        res_path = os.path.join(
            self.file_dir_venn, f"venn_bar_{name.replace(' ', '_')}" + FIG_FORMAT
        )
        fig.savefig(res_path, dpi=200, bbox_inches="tight")
        plt.close("all")

    @exception_handler
    def save_venn_names(self, named_sets: dict):
        if len(named_sets) > 6:
            self.logger.warning("Skipping save_venn_names because more than 6 experiments were passed at once")
            return

        for intersected, unioned, result in venn_names(named_sets):
            # create name based on the intersections and unions that were done
            intersected_name = "&".join(sorted(intersected))
            unioned_name = "-" + "-".join(sorted(unioned)) if unioned else ""
            res_path = os.path.join(
                self.file_dir_venn,
                f"set_{intersected_name}{unioned_name}.txt"
            )
            # write all names line by line into the file
            with open(res_path, "w") as out:
                for re in result:
                    out.write(re + "\n")

    def get_venn_group_data(self, df_to_use: str = "raw", level: int = 0, non_na_function=get_number_of_non_na_values):
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

    @exception_handler
    def plot_venn_groups(self, df_to_use: str = "raw", levels: tuple = (0,), show_suptitle=True):
        for level in levels:
            # create venn diagrams comparing all replicates within an experiment
            named_sets = self.get_venn_group_data(df_to_use, level)
            ex = f"group_level_{level}"
            # save the resulting venn diagram
            self.save_venn(ex, named_sets, show_suptitle=show_suptitle)
            # save the sets to txt files
            self.save_venn_names(named_sets)
            # create a mixture of bar and venn diagram
            self.save_bar_venn(ex, named_sets, show_suptitle=show_suptitle)

    @exception_handler
    def plot_venn_results(self, df_to_use: str = "raw", levels: tuple = (0,), show_suptitle=True):
        for level in levels:
            level_keys = self.all_tree_dict[df_to_use].level_keys_full_name[level]
            for key in level_keys:
                named_sets = self.get_venn_data_per_key(df_to_use, key)
                # save the resulting venn diagram
                self.save_venn(key, named_sets, show_suptitle=show_suptitle)
                # save the sets to txt files
                self.save_venn_names(named_sets)
                # create a mixture of bar and venn diagram
                self.save_bar_venn(key, named_sets, show_suptitle=show_suptitle)

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

    @exception_handler
    def plot_detection_counts(self, df_to_use: str = "raw", levels: Iterable = (0,), **kwargs):
        for level in levels:
            data = self.get_detection_counts_data(df_to_use=df_to_use, level=level, **kwargs)
            if data:
                matplotlib_plots.save_detection_counts_results(
                    **data, intensity_label=self.intensity_label_names[df_to_use], df_to_use=df_to_use,
                    level=level, save_path=self.file_dir_descriptive, **kwargs
                )

    @exception_handler
    def plot_number_of_detected_proteins(self, df_to_use: str = "raw", show_suptitle: bool = True, levels: Iterable = (0,)):
        for level in levels:
            plt.close("all")
            # determine number of rows and columns in the plot based on the number of experiments
            level_values = self.all_tree_dict[df_to_use].level_keys_full_name[level]
            # determine number of rows and columns in the plot based on the number of experiments
            n_rows_experiment, n_cols_experiment = get_number_rows_cols_for_fig(level_values)
            fig, axarr = plt.subplots(n_rows_experiment, n_cols_experiment,
                                      figsize=(5 * n_cols_experiment, 3 * n_rows_experiment))
            if show_suptitle:
                fig.suptitle(f"Number of detected proteins from {self.intensity_label_names[df_to_use]}")

            for experiment, ax in zip(level_values, axarr.flat):
                intensities = self.all_tree_dict[df_to_use][experiment].aggregate(None)

                # how many proteins were detected per replicate and in total
                counts = (intensities > 0).sum(axis=1)
                counts = counts[counts > 0]
                heights = [len(counts)]
                # labels start at 0 so we prepend one empty string
                labels = ["Total"]
                for col in intensities:
                    h = len(intensities[col][intensities[col] > 0])
                    heights.append(h)
                    labels.append(col)
                mean_height = np.mean(heights[1:])
                # self.logger.debug(labels)
                y_pos = [x for x in range(len(heights))]
                ax.barh(y_pos, heights, color="skyblue")

                for y, value in zip(y_pos, heights):
                    ax.text(heights[0] / 2, y, value,
                            verticalalignment='center', horizontalalignment='center')
                ax.set_title(experiment)
                ax.axvline(mean_height, linestyle="--", color="black", alpha=0.6)
                ax.set_yticks([i for i in range(len(labels))])
                ax.set_yticklabels(labels)
                ax.set_xlabel("Counts")
                #ax1.tick_params(rotation=70)
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            res_path = os.path.join(self.file_dir_descriptive, "detection_per_replicate" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")

    @exception_handler
    def plot_intensity_histograms(self, df_to_use: str = "raw", show_suptitle: bool = False, levels=()):
        # TODO levels should overlay histograms
        plt.close("all")
        columns = self.all_tree_dict[df_to_use].level_keys_full_name[max(self.all_tree_dict[df_to_use].level_keys_full_name.keys())]
        n_rows_replicates, n_cols_replicates = get_number_rows_cols_for_fig(columns)
        # make a intensity histogram for every replicate
        fig, axarr = plt.subplots(n_rows_replicates, n_cols_replicates,
                                  figsize=(5 * n_cols_replicates, 5 * n_rows_replicates))
        if show_suptitle:
            fig.suptitle(f"Replicate {self.intensity_label_names[df_to_use]} histograms")

        for col, ax in zip(columns, axarr.flat):
            intensities = self.all_tree_dict[df_to_use][col].aggregate()

            if "log2" in df_to_use:
                bins = np.linspace(intensities.min(), intensities.max(), 25)
            else:
                bins = np.logspace(np.log2(intensities.min()), np.log2(intensities.max()), 25, base=2)

            ax.set_title(f"{col}".replace(self.int_mapping[df_to_use], ""))
            ax.hist(intensities, bins=bins)
            if "log2" not in df_to_use:
                ax.set_xscale("log", basex=2)
            ax.set_xlabel(f"{self.intensity_label_names[df_to_use]}")
            ax.set_ylabel("Counts")

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        res_path = os.path.join(self.file_dir_descriptive, "intensity_histograms" + FIG_FORMAT)
        fig.savefig(res_path, dpi=200, bbox_inches="tight")

    @exception_handler
    def plot_scatter_replicates(self, df_to_use: str = "raw"):
        for experiment in self.replicates:
            plt.close("all")
            # scatter plot of all replicates vs all replicates
            fig, ax = plt.subplots(1, 1, figsize=(7, 7))
            for rep1, rep2 in combinations(self.replicates[experiment], 2):
                x1 = self.all_intensities_dict[df_to_use][self.int_mapping[df_to_use] + rep1]
                x2 = self.all_intensities_dict[df_to_use][self.int_mapping[df_to_use] + rep2]
                corr_mask = np.logical_and(~x1.isna(), ~x2.isna())
                plot_mask = np.logical_or(~x1.isna(), ~x2.isna())
                exp = r"$r^{2}$"
                ax.scatter(x1.fillna(x2.min() * 0.95)[plot_mask], x2.fillna(x2.min() * 0.95)[plot_mask], label=f"{rep1} vs {rep2}, "
                           fr"{exp}: {stats.pearsonr(x1[corr_mask], x2[corr_mask])[0] ** 2:.4f}",
                           alpha=0.5, marker=".")
                ax.set_xlabel(f"{self.intensity_label_names[df_to_use]}")
                ax.set_ylabel(f"{self.intensity_label_names[df_to_use]}")

            fig.legend(frameon=False)
            if "log2" not in df_to_use:
                ax.set_xscale("log")
                ax.set_yscale("log")

            res_path = os.path.join(self.file_dir_descriptive,
                                    f"{self.replicate_representation[experiment].replace(' ', '_')}_scatter" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")

    @exception_handler
    def plot_rank(self, df_to_use: str = "raw", levels: Iterable = (0,)):
        if self.interesting_proteins.values():
            all_pathway_proteins = set.union(*(set(x) for x in self.interesting_proteins.values()))
        else:
            all_pathway_proteins = set()
        for level in levels:
            for level_key in self.all_tree_dict[df_to_use].level_keys_full_name[level]:
                plt.close("all")
                # get the experiment intensities calculate mean intensity for the experiment and sort from highest to lowest
                m_intensity = self.all_tree_dict[df_to_use][level_key].aggregate().sort_values(ascending=False)
                # TODO apply filter for rare proteins before here?
                # protein ranks vs intensity
                # create dict to map each protein its respective rank and mean intensity
                dic = {idx: (i, value) for i, (idx, value) in enumerate(m_intensity.items())}

                found_proteins = set(m_intensity.index)
                # get all proteins that are not part of any pathway
                non_pathway_proteins = found_proteins - all_pathway_proteins
                # get all proteins that are part of any pathway
                pathway_proteins = found_proteins & all_pathway_proteins
                rank_identified_proteins = [dic[protein][0] for protein in pathway_proteins]
                # plot the non pathway proteins
                x = [dic[protein][0] for protein in non_pathway_proteins]
                y = [dic[protein][1] for protein in non_pathway_proteins]

                fig, ax = plt.subplots(1, 1, figsize=(14, 7))
                ax.scatter(x, y, c=f"darkgray", s=10, alpha=0.3, marker=".", label="no pathway")
                # plot all proteins of a specific pathway
                for i, (pathway, proteins) in enumerate(self.interesting_proteins.items()):
                    proteins = set(proteins) & found_proteins
                    x = [dic[protein][0] for protein in proteins]
                    y = [dic[protein][1] for protein in proteins]
                    ax.scatter(x, y, c=f"C{i}", s=80, alpha=0.6, marker=".", label=pathway)
                #
                # only if more than 0 proteins are identified
                if rank_identified_proteins:
                    median_pathway_rank = int(np.median(rank_identified_proteins))
                    median_intensity = m_intensity.iloc[median_pathway_rank]
                    xmin, xmax = ax.get_xbound()
                    xm = (median_pathway_rank + abs(xmin)) / (abs(xmax) + abs(xmin))
                    ymin, ymax = ax.get_ybound()
                    ym = ((median_intensity) - (ymin) ) / ((ymax) - (ymin))
                    # plot the median rank and intensity at that rank
                    ax.axvline(median_pathway_rank, ymax=ym, linestyle="--", color="black", alpha=0.6)
                    ax.axhline(median_intensity, xmax=xm, linestyle="--", color="black", alpha=0.6)
                    ax.text(xmin * 0.9, median_intensity * 0.9,
                            f"median rank: {median_pathway_rank} ({median_pathway_rank/len(m_intensity) * 100 :.1f}%) "
                            f"with intensity: {median_intensity:.2E}",  # TODO case for log and non log
                            verticalalignment="top", horizontalalignment="left")
                    # ax.annotate(f"median rank: {median_pathway_rank} with intensity: {median_intensity}",
                    #             xy=(median_pathway_rank, median_intensity), xycoords='data',
                    #             xytext=(0.1, 0.1), textcoords='axes fraction',
                    #             arrowprops=dict(facecolor='black', shrink=0.05),
                    #             horizontalalignment='right', verticalalignment='top',
                    #             )

                ax.set_xlabel("Protein rank")
                ax.set_ylabel(f"{level_key} mean")
                if "log2" not in df_to_use:
                    ax.set_yscale("log")
                fig.legend(bbox_to_anchor=(1.02, 0.5), loc="center left")

                res_path = os.path.join(self.file_dir_descriptive,
                                        f"{level_key}_rank" + FIG_FORMAT)
                fig.savefig(res_path, dpi=200, bbox_inches="tight")

    def get_relative_std_data(self, df_to_use: str, full_name: str, **kwargs):
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
        return {"intensities": intensities}

    @exception_handler
    def plot_relative_std(self, df_to_use: str = "raw", levels: Iterable = (0,), **kwargs):
        # TODO check with log2 thresholds
        for level in levels:
            for full_name in self.all_tree_dict[df_to_use].level_keys_full_name[level]:
                data = self.get_relative_std_data(df_to_use=df_to_use, full_name=full_name, **kwargs)
                if data:
                    matplotlib_plots.save_relative_std_results(
                        **data, df_to_use=df_to_use, name=full_name, save_path=self.file_dir_descriptive,
                        intensity_label=self.intensity_label_names[df_to_use], **kwargs
                    )

    def get_pathway_analysis_data(self, df_to_use, level, pathway, equal_var=True, **kwargs):
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

    @exception_handler
    def plot_pathway_analysis(self, df_to_use, levels, **kwargs):
        for level in levels:
            for pathway in list(self.interesting_proteins.keys()):
                data = self.get_pathway_analysis_data(level=level, df_to_use=df_to_use, pathway=pathway, **kwargs)
                if data:
                    matplotlib_plots.save_pathway_analysis_results(**data, pathway=pathway, level=level, intensity_label=self.intensity_label_names[df_to_use], save_path=self.file_dir_pathway, **kwargs)

    @exception_handler
    def plot_pathway_timeline(self, df_to_use: str = "raw", show_suptitle: bool = False, levels: Iterable = (2,)):
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

    @exception_handler
    def plot_experiment_comparison(self, df_to_use: str = "raw", levels: Iterable = (0,)):
        # TODO correlation of log2 and not log 2 data is different
        for level in levels:
            for ex1, ex2 in combinations(self.all_tree_dict[df_to_use].level_keys_full_name[level], 2):
                plt.close("all")
                protein_intensities_ex1 = self.all_tree_dict[df_to_use][ex1].aggregate(None)
                protein_intensities_ex2 = self.all_tree_dict[df_to_use][ex2].aggregate(None)
                mask, exclusive_1, exclusive_2 = get_intersection_and_unique(protein_intensities_ex1, protein_intensities_ex2)
                # flatten all replicates
                exclusive_ex1 = protein_intensities_ex1[exclusive_1].mean(axis=1)
                exclusive_ex2 = protein_intensities_ex2[exclusive_2].mean(axis=1)
                protein_intensities_ex1 = protein_intensities_ex1[mask].mean(axis=1)
                protein_intensities_ex2 = protein_intensities_ex2[mask].mean(axis=1)
                # calculate r
                try:
                    r = stats.pearsonr(protein_intensities_ex1, protein_intensities_ex2)
                except ValueError:
                    self.logger.warning("Could not calculate pearson r for %s vs %s", ex1, ex2)
                    r = (np.nan, )

                fig, ax = plt.subplots(1, 1, figsize=(7, 7))
                exp = r"$r^{2}$"
                ax.scatter(protein_intensities_ex1, protein_intensities_ex2, s=8, alpha=0.6, marker=".",
                           label=f"{ex1} vs {ex2}, {exp}: {r[0] ** 2:.4f}")
                ax.scatter(exclusive_ex1, [np.min(protein_intensities_ex2) * 0.95] * exclusive_ex1.shape[0],
                           s=8, alpha=0.6, marker=".", label=f"exclusive for {ex1}")
                ax.scatter([np.min(protein_intensities_ex1) * 0.95] * exclusive_ex2.shape[0], exclusive_ex2,
                           s=8, alpha=0.6, marker=".", label=f"exclusive for {ex2}")

                ax.set_xlabel(ex1)
                ax.set_ylabel(ex2)
                fig.legend(frameon=False)
                if "log2" not in df_to_use:
                    ax.set_xscale("log")
                    ax.set_yscale("log")

                xmin, xmax = ax.get_xbound()
                ymin, ymax = ax.get_ybound()
                ax.set_xlim(min(xmin, ymin), max(xmax, ymax))
                ax.set_ylim(min(xmin, ymin), max(xmax, ymax))

                res_path = os.path.join(self.file_dir_descriptive,
                                        f"{ex1}_vs_{ex2}" + FIG_FORMAT)
                fig.savefig(res_path, dpi=200, bbox_inches="tight")
                plt.close("all")

    @exception_handler
    def plot_go_analysis(self, df_to_use: str = "raw", levels: Iterable = (0,)):
        plt.close("all")

        # all proteins that were detected in any replicate
        background = set(self.all_intensities_dict[df_to_use].index)
        heights = ddict(list)
        test_results = ddict(list)

        for level in levels:
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

            # TODO only show pvalue when significant
            # TODO also create table
            # TODO move labels to bars, remove legend
            fig, ax = plt.subplots(1, 1, figsize=(7, int(len(heights) * len(self.go_analysis_gene_names) / 3)))

            bar_width = 0.25
            for i, experiment in enumerate(heights):
                y_pos = np.array([(i + (len(heights) + 1) * x) * bar_width for x in range(len(self.go_analysis_gene_names))])
                ax.barh(y_pos, heights[experiment], height=bar_width, edgecolor='white', label=experiment)
                for x, y, text in zip(heights[experiment], y_pos, test_results[experiment]):
                    if text > 0.05:
                        continue
                    text = f"{text:.4f}" if text > 0.0005 else "< 0.0005"
                    ax.text(x, y, f" p: {text}", verticalalignment="center", fontsize="8")

            ax.set_ylabel('compartiment')
            ax.set_xlabel('number of proteins')
            # add yticks on the middle of the group bars
            ax.set_yticks([(len(heights) / 2 + (len(heights) + 1) * x) * bar_width
                           for x in range(len(self.go_analysis_gene_names))])
            # replace the y ticks with the compartiment names
            ax.set_yticklabels([x for x in self.go_analysis_gene_names])
            plt.legend()
            res_path = os.path.join(self.file_dir_go_analysis, f"go_analysis_level_{level}" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")

    def get_r_volcano_data(self, g1: str, g2: str, df_to_use: str):
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

    @exception_handler
    def plot_r_volcano(self, df_to_use: str, levels, **kwargs):  # TODO Iterable[int]
        # TODO both adj and un adj should be available
        for level in levels:
            level_keys = self.all_tree_dict[df_to_use].level_keys_full_name[level]
            for g1, g2 in combinations(level_keys, 2):
                data = self.get_r_volcano_data(g1, g2, df_to_use)
                if data:
                    matplotlib_plots.save_volcano_results(
                        **data, g1=g1, g2=g2, save_path=self.file_dir_volcano,
                        intensity_label=self.intensity_label_names[df_to_use], **kwargs
                    )

    def get_pca_data(self, df_to_use: str = "raw_log2", level: int = 0, n_components: int = 4, fill_value: float = 0, fill_na_before_norm: bool = False, **kwargs):
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

    @exception_handler
    def plot_pca_overview(self, df_to_use: str, levels, **kwargs):  # TODO Iterable[int]
        for level in levels:
            data = self.get_pca_data(level=level, df_to_use=df_to_use, **kwargs)
            if data:
                matplotlib_plots.save_pca_results(**data, save_path=self.file_dir_descriptive, **kwargs)

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

    def plot_boxplot(self, df_to_use, levels, **kwargs):
        for level in levels:
            data = self.get_boxplot_data(df_to_use=df_to_use, level=level, **kwargs)
            if data:
                matplotlib_plots.save_boxplot_results(
                    **data, level=level, intensity_label=self.intensity_label_names[df_to_use],
                    save_path=self.file_dir_descriptive, **kwargs
                )


class MaxQuantPlotter(MSPPlots):
    def __init__(
            self,
            start_dir: str,
            reader_data: dict,
            intensity_df_name: str = "proteinGroups.txt",
            interesting_proteins: dict = None,
            go_analysis_gene_names: dict = None,
            configs: dict = None,
            required_reader="mqreader",
            intensity_entries=(("raw", "Intensity ", "Intensity"), ("lfq", "LFQ intensity ", "LFQ intensity"), ("ibaq", "iBAQ ", "iBAQ intensity")),
            loglevel=logging.DEBUG
    ):
        super().__init__(
            start_dir,
            reader_data,
            intensity_df_name,
            interesting_proteins,
            go_analysis_gene_names,
            configs,
            required_reader,
            intensity_entries,
            loglevel
        )

    @classmethod
    def from_MSPInitializer(cls, mspinti_instance: MSPInitializer, intensity_entries=(("raw", "Intensity ", "Intensity"), ("lfq", "LFQ intensity ", "LFQ intensity"), ("ibaq", "iBAQ ", "iBAQ intensity"))):
        return cls(
            start_dir=mspinti_instance.start_dir,
            reader_data=mspinti_instance.reader_data,
            intensity_df_name="proteinGroups.txt",
            interesting_proteins=mspinti_instance.interesting_proteins,
            go_analysis_gene_names=mspinti_instance.go_analysis_gene_names,
            configs=mspinti_instance.configs,
            required_reader="mqreader",
            intensity_entries=intensity_entries,
            loglevel=mspinti_instance.logger.getEffectiveLevel()
            )

    def create_report(self):
        def bar_from_counts(ax, counts, title=None, relative=False, yscale=None):
            if relative:
                counts = counts / counts.sum()
            if title is not None:
                ax.set_title(title)
            ax.bar([x for x in range(len(counts))], counts.values)
            ax.set_xticks([i for i in range(len(counts))])
            ax.set_xticklabels(counts.index)
            if yscale is not None:
                if isinstance(yscale, str):
                    ax.set_yscale(yscale)
                elif isinstance(yscale, dict):
                    ax.set_yscale(**yscale)

        with PdfPages(os.path.join(self.start_dir, "MaxQuantReport2.pdf")) as pdf:
            first_page = plt.figure(figsize=(14, 7))
            text_conf = dict(transform=first_page.transFigure, size=24, ha="center")
            first_page.text(0.5, 0.92, "MaxQuant report", **text_conf)
            text_conf.update({"size": 20})
            first_page.text(0.5, 0.85, "parameter.txt info", **text_conf)
            text_conf.pop("size")
            first_page.text(0.5, 0.8, f"Version: {self.required_reader_data['parameters.txt']['Version']}, "
                            f"run at: {self.required_reader_data['parameters.txt']['Date of writing']}", **text_conf)
            first_page.text(0.5, 0.75, f"Fasta File: {os.path.split(self.required_reader_data['parameters.txt']['Fasta file'])[1]}, "
                                       f"Match between runs: {self.required_reader_data['parameters.txt']['Match between runs']}", **text_conf)
            first_page.text(0.5, 0.7, "Min. to Max. peptide length for unspecific search: "
                                      f"{self.required_reader_data['parameters.txt']['Min. peptide length for unspecific search']} to {self.required_reader_data['parameters.txt']['Max. peptide length for unspecific search']}", **text_conf)

            text_conf.update({"size": 20})
            first_page.text(0.5, 0.65, "summary.txt info", **text_conf)
            text_conf.pop("size")
            first_page.text(0.5, 0.6, f"Used Enzyme: {self.required_reader_data['summary.txt'].loc[1, 'Enzyme']}", **text_conf)
            first_page.text(0.5, 0.55, f"Variable modifications: {self.required_reader_data['summary.txt'].loc[1, 'Variable modifications']}", **text_conf)
            first_page.text(0.5, 0.5, f"Mass Standard Deviation: mean {self.required_reader_data['summary.txt'].loc[:, 'Mass Standard Deviation [ppm]'].mean()}, max {self.required_reader_data['summary.txt'].loc[:, 'Mass Standard Deviation [ppm]'].max()}", **text_conf)

            # TODO LFQ, Identified proteins, and peptides
            pdf.savefig()
            plt.close(first_page)

            fig, axarr = plt.subplots(3, 1, figsize=(14, 7))
            bar_from_counts(axarr[0], self.required_reader_data["peptides.txt"]["Missed cleavages"].value_counts(), title="Missed Cleavages", relative=True)
            bar_from_counts(axarr[1], self.required_reader_data["peptides.txt"]["Amino acid before"].value_counts(), title="Amino acid before", yscale="log")
            bar_from_counts(axarr[2], self.required_reader_data["peptides.txt"]["Last amino acid"].value_counts(), title="Last amino acid", yscale="log")

            fig.tight_layout(rect=[0, 0.03, 1, 0.95])

            pdf.savefig()
            plt.close(fig)

            fig, axarr = plt.subplots(3, 1, figsize=(14, 7))
            bar_from_counts(axarr[0], self.required_reader_data["peptides.txt"]["Charges"].value_counts(), title="Peptide Charges")

            fig.tight_layout(rect=[0, 0.03, 1, 0.95])

            pdf.savefig()
            plt.close(fig)
