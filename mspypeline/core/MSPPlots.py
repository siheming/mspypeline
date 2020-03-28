import pandas as pd
import numpy as np
import os
import functools
from matplotlib_venn import venn3, venn2
import matplotlib.pyplot as plt
from scipy import stats
from itertools import combinations
from collections import defaultdict as ddict
from collections.abc import Iterable
import logging
import warnings
from typing import Dict
from adjustText import adjust_text
import matplotlib.gridspec as gridspec

from mspypeline import MSPInitializer
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
        "plot_venn_groups", "plot_r_volcano"
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
        self, start_dir, configs, reader_data,
        interesting_proteins, go_analysis_gene_names,
        additional_entries=(),
        loglevel=logging.DEBUG
    ):
        self.logger = get_logger(self.__class__.__name__, loglevel=loglevel)
        # general information
        self.start_dir = start_dir
        self.configs = configs
        self.int_mapping = {}
        self.intensity_label_names = {}
        self.all_intensities_dict: Dict[str, pd.DataFrame] = {}
        # TODO atm only one reader is being used. All required information is read here atm
        df_protein_names = reader_data["mqreader"]["proteinGroups.txt"]
        tech_rep = self.configs["has_replicates"]
        self.all_replicates = self.configs["all_replicates"]
        self.analysis_design = self.configs["analysis_design"]
        # data frames
        self.df_protein_names = df_protein_names.set_index(df_protein_names["Gene name fasta"], drop=False)
        # self.df_peptide_names = df_peptide_names
        # dicts
        self.interesting_proteins = interesting_proteins
        self.go_analysis_gene_names = go_analysis_gene_names

        self.add_intensity_column("raw", "Intensity ", "Intensity")
        self.add_intensity_column("lfq", "LFQ intensity ", "LFQ intensity")
        self.add_intensity_column("ibaq", "iBAQ ", "iBAQ intensity")
        for option_name, name_in_file, name_in_plot in additional_entries:
            self.add_intensity_column(option_name, name_in_file, name_in_plot)

        # create the data trees
        self.all_tree_dict: Dict[str, DataTree] = {
            intensity: DataTree.from_analysis_design(self.analysis_design, data, tech_rep)
            for intensity, data in self.all_intensities_dict.items()
        }

        # set all result dirs
        # create file structure and folders
        # TODO: for now just create all of them
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
    def from_MSPInitializer(cls, mspinti_instance: MSPInitializer, additional_entries=()):
        return cls(
            start_dir = mspinti_instance.start_dir,
            configs = mspinti_instance.configs,
            reader_data = mspinti_instance.reader_data,
            interesting_proteins = mspinti_instance.interesting_proteins,
            go_analysis_gene_names = mspinti_instance.go_analysis_gene_names,
            additional_entries = additional_entries,
            loglevel = mspinti_instance.logger.getEffectiveLevel()
            )

    def create_results(self):
        for plot_name in MSPPlots.possible_plots:
            plot_settings_name = plot_name + "_settings"
            plot_settings = self.configs.get(plot_settings_name, {})
            if plot_settings.pop("create_plot", False):
                self.logger.debug(f"creating plot {plot_name}")
                getattr(self, plot_name)(**plot_settings)
        self.logger.info("Done creating plots")

    def add_intensity_column(self, option_name, name_in_file, name_in_plot):
        if not any((col.startswith(name_in_file) for col in self.df_protein_names)):
            self.logger.warning("%s columns could not be found in data", name_in_file)
            return
        self.int_mapping.update({option_name: name_in_file, f"{option_name}_log2": name_in_file})
        self.intensity_label_names.update({option_name: name_in_plot, f"{option_name}_log2": rf"$Log_2$ {name_in_plot}"})

        # extract all raw intensities from the dataframe
        # replace all 0 with nan and remove the prefix from the columns
        intensities = self.df_protein_names.loc[:,
            [f"{self.int_mapping[option_name]}{rep}" for rep in self.all_replicates]
        ].replace({0: np.nan}).rename(lambda x: x.replace(self.int_mapping[option_name], ""), axis=1)
        # filter all rows where all intensities are nan
        mask = (~intensities.isna()).sum(axis=1) != 0
        intensities = intensities[mask]

        # add log2 intensities
        intensities_log2 = np.log2(intensities)

        self.all_intensities_dict.update({
            option_name: intensities, f"{option_name}_log2": intensities_log2,
        })

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

    @exception_handler
    def plot_detection_counts(self, df_to_use: str = "raw", show_suptitle: bool = True, levels: Iterable = (0,)):
        for level in levels:
            plt.close("all")
            # determine number of rows and columns in the plot based on the number of experiments
            level_values = self.all_tree_dict[df_to_use].level_keys_full_name[level]
            n_rows_experiment, n_cols_experiment = get_number_rows_cols_for_fig(level_values)
            fig, axarr = plt.subplots(n_rows_experiment, n_cols_experiment, squeeze=True,
                                      figsize=(5 * n_cols_experiment, 3 * n_rows_experiment))
            if show_suptitle:
                fig.suptitle(f"Detection counts from {self.intensity_label_names[df_to_use]} intensities")
            for experiment, ax in zip(level_values, axarr.flat):
                intensities = self.all_tree_dict[df_to_use][experiment].aggregate(None)
                # from 0 to number of replicates, how often was each protein detected
                counts = (intensities > 0).sum(axis=1)
                counts = counts[counts > 0]
                heights = [len(counts[counts == value]) for value in sorted(counts.unique())]
                # y_pos = [value for value in sorted(counts.unique())]  # old version
                y_pos = [f"detected in {value} replicates" for value in sorted(counts.unique())]
                max_val = max(heights)

                ax.set_title(f"{experiment},\ntotal detected: {len(counts)}")
                ax.barh(y_pos, heights, color="skyblue")
                for y, value in zip(y_pos, heights):
                    ax.text(max_val / 2, y, value,
                            verticalalignment='center', horizontalalignment='center')
                ax.set_yticks(y_pos)
                ax.set_xlabel("Counts")

            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            res_path = os.path.join(self.file_dir_descriptive, f"detected_counts" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")

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
        n_rows_replicates, n_cols_replicates = get_number_rows_cols_for_fig(self.all_replicates)
        # make a intensity histogram for every replicate
        fig, axarr = plt.subplots(n_rows_replicates, n_cols_replicates,
                                  figsize=(5 * n_cols_replicates, 5 * n_rows_replicates))
        if show_suptitle:
            fig.suptitle(f"Replicate {self.intensity_label_names[df_to_use]} histograms")

        for col, ax in zip(self.all_replicates, axarr.flat):
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

    @exception_handler
    def plot_relative_std(self, df_to_use: str = "raw", levels: Iterable = (0,)):
        plt.close("all")
        # TODO check with log2 thresholds
        for level in levels:
            for full_name in self.all_tree_dict[df_to_use].level_keys_full_name[level]:
                intensities = self.all_tree_dict[df_to_use][full_name].aggregate(None)
                non_na = get_number_of_non_na_values(intensities.shape[1])
                mask = (intensities > 0).sum(axis=1) >= non_na
                intensities = intensities[mask]
                relative_std_percent = intensities.std(axis=1) / intensities.mean(axis=1) * 100

                bins = np.array([10, 20, 30])
                if "log2" in df_to_use:
                    bins = np.log2(bins)
                inds = np.digitize(relative_std_percent, bins).astype(int)

                cm = {0: "navy", 1: "royalblue", 2: "skyblue", 3: "darkgray"}
                colors = pd.Series([cm.get(x, "black") for x in inds], index=relative_std_percent.index)
                color_counts = {color: (colors == color).sum() for color in colors.unique()}

                # intensity vs relative standard deviation
                fig, ax = plt.subplots(1, 1, figsize=(14, 7))
                ax.scatter(intensities.mean(axis=1), relative_std_percent, c=colors, marker="o", s=(2* 72./fig.dpi)**2, alpha=0.8)
                ax.set_xlabel(f"Mean {self.intensity_label_names[df_to_use]}")
                ax.set_ylabel("Relative Standard deviation [%]")
                if "log2" not in df_to_use:
                    ax.set_xscale('log')
                xmin, xmax = ax.get_xbound()
                cumulative_count = 0
                for i, bin_ in enumerate(bins):
                    cumulative_count += color_counts.get(cm[i], 0)
                    ax.axhline(bin_, color=cm[i])
                    ax.text(xmin, bin_, cumulative_count)

                res_path = os.path.join(self.file_dir_descriptive,
                                        f"{full_name}_rel_std" + FIG_FORMAT)
                fig.savefig(res_path, dpi=200, bbox_inches="tight")
                plt.close(fig)

    @exception_handler
    def plot_pathway_analysis(self, df_to_use: str = "raw", show_suptitle: bool = False, levels: Iterable = (0,)):
        plot_ns = False
        threshhold = 0.05
        for level in levels:
            level_keys = self.all_tree_dict[df_to_use].level_keys_full_name[level]
            ex_list = sorted(level_keys, reverse=True)
            for pathway in self.interesting_proteins:
                plt.close("all")
                found_proteins = set(self.interesting_proteins[pathway])
                found_proteins &= set(self.all_intensities_dict[df_to_use].index)
                found_proteins = sorted(list(found_proteins))
                if len(found_proteins) < 1:
                    self.logger.warning("Skipping pathway %s in pathway analysis because no proteins were found", pathway)
                    continue
                n_rows, n_cols = get_number_rows_cols_for_fig(found_proteins)
                fig, axarr = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, int(n_rows * len(level_keys) / 1.5)))
                if show_suptitle:
                    fig.suptitle(pathway)
                all_heights = {}
                try:
                    axiterator = axarr.flat
                except AttributeError:
                    axiterator = [axarr]
                for protein, ax in zip(found_proteins, axiterator):
                    ax.set_title(protein)
                    heights = []
                    for idx, level_key in enumerate(sorted(level_keys, reverse=True)):
                        protein_intensities = self.all_tree_dict[df_to_use][level_key].aggregate(None, index=protein)
                        ax.scatter(protein_intensities, [idx] * len(protein_intensities), label=f"{level_key}")
                        heights.append(np.max(protein_intensities))
                    all_heights[protein] = heights
                    ax.set_ylim((-1, len(level_keys)))

                    ax.set_yticks([i for i in range(len(level_keys))])
                    # ax.set_yticklabels([self.replicate_representation[r] for r in sorted(level_keys, reverse=True)])
                    ax.set_yticklabels(sorted(level_keys, reverse=True))
                    ax.set_xlabel(self.intensity_label_names[df_to_use])
                fig.tight_layout(rect=[0, 0.03, 1, 0.95])
                res_path = os.path.join(self.file_dir_pathway, f"pathway_analysis_{pathway}_no_labels" + FIG_FORMAT)
                fig.savefig(res_path, dpi=200, bbox_inches="tight")
                try:
                    axiterator = axarr.flat
                except AttributeError:
                    axiterator = [axarr]
                significances = []
                for protein, ax in zip(found_proteins, axiterator):
                    per_protein_significant = []
                    for e1, e2 in combinations(level_keys, 2):
                        v1 = self.all_tree_dict[df_to_use][e1].aggregate(None, index=protein)
                        v2 = self.all_tree_dict[df_to_use][e2].aggregate(None, index=protein)
                        # filter entries with too many nans based on function
                        non_na_group_1 = get_number_of_non_na_values(v1.shape[0])
                        non_na_group_2 = get_number_of_non_na_values(v2.shape[0])
                        mask_1 = (v1 > 0).sum(axis=0) >= non_na_group_1
                        mask_2 = (v2 > 0).sum(axis=0) >= non_na_group_2
                        mask = np.logical_and(mask_1, mask_2)
                        if not mask:
                            continue
                        test = stats.ttest_ind(v1[~v1.isna()], v2[~v2.isna()], equal_var=self.configs.get("equal_var", False))
                        if plot_ns or test[1] <= threshhold:
                            significances.append((protein, e1, e2, test[1]))
                            per_protein_significant.append((e1, e2, test[1]))
                    per_protein_df = pd.DataFrame(per_protein_significant, columns=["experiment1", "experiment2", "pvalue"])
                    per_protein_df = per_protein_df.sort_values("pvalue")#.head(20)
                    # adjust axis height
                    xmin, xmax = ax.get_xbound()
                    ax.set_xlim(right=xmax * (1 + per_protein_df.shape[0] * 0.015))
                    for i, (e1, e2, pval) in enumerate(zip(per_protein_df["experiment1"], per_protein_df["experiment2"], per_protein_df["pvalue"])):
                        plot_annotate_line(ax, ex_list.index(e1), ex_list.index(e2), xmax * (1 + i * 0.015) - 0.005, pval)
                df = pd.DataFrame(significances, columns=["protein", "experiment1", "experiment2", "pvalue"])
                df.to_csv(os.path.join(self.file_dir_pathway, f"pathway_analysis_{pathway}_table.csv"), index=False)

                fig.tight_layout(rect=[0, 0.03, 1, 0.95])
                res_path = os.path.join(self.file_dir_pathway, f"pathway_analysis_{pathway}" + FIG_FORMAT)
                fig.savefig(res_path, dpi=200, bbox_inches="tight")

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

    def create_report(self):
        pass

    @exception_handler
    def plot_go_analysis(self, df_to_use: str = "raw", levels: Iterable = (0,)):
        self.logger.info("Creating go analysis plots")
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

    @exception_handler
    def plot_r_volcano(self, df_to_use: str = "lfq_log2", show_suptitle: bool = True, levels: Iterable = (1,), p_value="pval", fchange_threshold=2, scatter_size=10):
        # TODO both adj and un adj should be available
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

        for level in levels:
            level_keys = self.all_tree_dict[df_to_use].level_keys_full_name[level]
            for g1, g2 in combinations(level_keys, 2):
                plt.close("all")
                # get data of the samples
                v1 = self.all_tree_dict[df_to_use][g1].aggregate(None)
                v2 = self.all_tree_dict[df_to_use][g2].aggregate(None)
                if v2.shape[1] < 2 or v2.shape[1] < 2:
                    self.logger.warning("Skipping Volcano plot for comparison: %s, %s because the groups contain only "
                                        "%s and %s experiments", g1, g2, v1.shape[1], v2.shape[1])
                    continue
                mask, exclusive_1, exclusive_2 = get_intersection_and_unique(v1, v2)

                df = pd.concat([v1[mask], v2[mask]], axis=1)
                design = pd.DataFrame([[0] * v1.shape[1] + [1] * v2.shape[1],
                                       [1] * v1.shape[1] + [0] * v2.shape[1]], index=[g2, g1]).T

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
                sig_to_col = {"pval": "P.Value", "adjpval": "adj.P.Val"}
                col = sig_to_col[p_value]
                col_mapping = {"adj.P.Val": "adjusted p value", "P.Value": "unadjusted p value"}
                plot_data = ress.loc[:, ["logFC", "AveExpr", "P.Value", "adj.P.Val"]]
                # save all values
                plot_data.to_csv(os.path.join(self.file_dir_volcano, f"volcano_plot_data_{g1}_vs_{g2}_full_{col_mapping[col].replace(' ', '_')}.csv"))
                # save significant values
                plot_data[plot_data[col] < 0.05].to_csv(
                    os.path.join(self.file_dir_volcano, f"volcano_plot_data_{g1}_vs_{g2}_significant_{col_mapping[col].replace(' ', '_')}.csv"))
                # calculate mean intensity for unique genes for plotting
                unique_g1 = v1[exclusive_1].mean(axis=1).rename(f"{df_to_use} mean intensity")
                unique_g2 = v2[exclusive_2].mean(axis=1).rename(f"{df_to_use} mean intensity")
                # save unique values
                unique_g1.to_csv(os.path.join(self.file_dir_volcano, f"volcano_plot_data_{g1}_vs_{g2}_unique_{g1}.csv"), header=True)
                unique_g2.to_csv(os.path.join(self.file_dir_volcano, f"volcano_plot_data_{g1}_vs_{g2}_unique_{g2}.csv"), header=True)

                def get_volcano_significances(fchange, pval, fchange_threshold):
                    if pval > 0.05 or abs(fchange) < np.log2(fchange_threshold):
                        return "ns"
                    elif fchange >= 0:
                        return "up"
                    elif fchange < 0:
                        return "down"
                    else:
                        raise ValueError("heisenbug")

                significance_to_color = {"ns": "gray", "up": "red", "down": "blue"}
                significance_to_label = {"ns": "non-significant", "up": f"upregulated in {g2}", "down": f"upregulated in {g1}"}

                # plot
                fig = plt.figure(figsize=(7, 7))

                gs = gridspec.GridSpec(1, 3, width_ratios=[1, 8, 1])
                ax_unique_down = plt.subplot(gs[0])
                ax = plt.subplot(gs[1])
                ax_unique_up = plt.subplot(gs[2])

                # hide the spines between ax and ax2
                ax_unique_down.spines['right'].set_visible(False)
                ax_unique_up.spines['left'].set_visible(False)
                ax_unique_down.yaxis.tick_left()
                ax_unique_up.yaxis.tick_right()
                ax_unique_up.yaxis.set_label_position("right")
                # hide the xticks
                ax_unique_down.tick_params(which='both', bottom=False, labelbottom=False)
                ax_unique_up.tick_params(which='both', bottom=False, labelbottom=False)

                # non sign gray, left side significant blue right side red
                significances_plot = [get_volcano_significances(log_fold_change, p_val, fchange_threshold)
                                      for log_fold_change, p_val in zip(plot_data["logFC"], plot_data[col])]
                for regulation in significance_to_color:
                    mask = [x == regulation for x in significances_plot]
                    ax.scatter(plot_data["logFC"][mask], -np.log10(plot_data[col])[mask], s=scatter_size,
                               color=significance_to_color[regulation], label=f"{sum(mask)} {significance_to_label[regulation]}")
                ymin, ymax = ax.get_ybound()
                xmin, xmax = ax.get_xbound()
                m = max(abs(xmin), xmax)
                ax.set_xlim(left=-1 * m, right=m)
                xmin, xmax = ax.get_xbound()
                x_offset = (np.log2(fchange_threshold) / xmax) / 2
                # add line at significance threshold
                axline_kwargs = dict(linestyle="--", color="black", alpha=0.5, linewidth=1)
                if any(plot_data[col] < 0.05):
                    ax.axhline(-np.log10(0.05), **axline_kwargs, xmin=0, xmax=0.5 - x_offset)
                    ax.axhline(-np.log10(0.05), **axline_kwargs, xmin=0.5 + x_offset, xmax=1)

                # add lines for minimum fold change threshold
                y_percentage = (-np.log10(0.05) + abs(ymin)) / (ymax + abs(ymin))
                if fchange_threshold > 0:
                    ax.axvline(-np.log2(fchange_threshold), **axline_kwargs, ymin=y_percentage, ymax=1)
                    ax.axvline(np.log2(fchange_threshold), **axline_kwargs, ymin=y_percentage, ymax=1)
                # plot unique values with mean intensity at over maximum
                ax_unique_down.scatter([0] * len(unique_g1), unique_g1, s=scatter_size, color="dodgerblue", label=f"{len(unique_g1)} unique in {g1}")
                ax_unique_up.scatter([0] * len(unique_g2), unique_g2, s=scatter_size, color="coral", label=f"{len(unique_g2)} unique in {g2}")
                # adjust bounds
                ymin_down, ymax_down = ax_unique_down.get_ybound()
                ymin_up, ymax_up = ax_unique_up.get_ybound()
                ax_unique_down.set_ylim(bottom=min(ymin_down, ymin_up), top=max(ymax_down, ymax_up))
                ax_unique_up.set_ylim(bottom=min(ymin_down, ymin_up), top=max(ymax_down, ymax_up))

                # figure stuff
                if show_suptitle:
                    fig.suptitle(f"{g1} vs {g2}")
                ax.set_xlabel(r"$Log_2$ Fold Change")
                ax.set_ylabel(r"-$Log_{10}$" + f" {col_mapping[col]}")
                ax_unique_down.set_ylabel(self.intensity_label_names[df_to_use])
                ax_unique_up.set_ylabel(self.intensity_label_names[df_to_use])
                fig.legend(bbox_to_anchor=(1.02, 0.5), loc="center left", frameon=False)
                res_path = os.path.join(self.file_dir_volcano, f"volcano_{g1}_{g2}_no_annotation_{col_mapping[col].replace(' ', '_')}" + FIG_FORMAT)
                fig.tight_layout(rect=[0, 0.03, 1, 0.95])
                fig.savefig(res_path, dpi=200, bbox_inches="tight")
                significant_upregulated = plot_data[(plot_data["logFC"] > np.log2(fchange_threshold)) & (plot_data[col] < 0.05)].sort_values(by=[col], ascending=True).head(10)
                significant_downregulated = plot_data[(plot_data["logFC"] < -np.log2(fchange_threshold)) & (plot_data[col] < 0.05)].sort_values(by=[col], ascending=True).head(10)
                significant = pd.concat([significant_upregulated, significant_downregulated])
                texts = []
                for log_fold_change, p_val, gene_name in zip(significant["logFC"], significant[col], significant.index):
                    texts.append(ax.text(log_fold_change, -np.log10(p_val), gene_name, ha="center", va="center", fontsize=8))
                adjust_text(texts, arrowprops=dict(width=0.15, headwidth=0, color='gray', alpha=0.6), ax=ax)
                res_path = os.path.join(self.file_dir_volcano, f"volcano_{g1}_{g2}_annotation_{col_mapping[col].replace(' ', '_')}" + FIG_FORMAT)
                fig.savefig(res_path, dpi=200, bbox_inches="tight")
                # TODO scatter plot of significant genes
