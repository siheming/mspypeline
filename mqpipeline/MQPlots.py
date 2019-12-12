import pandas as pd
import numpy as np
import os
import functools
from mqpipeline.Logger import Logger
from matplotlib_venn import venn3, venn2
import matplotlib.pyplot as plt
from scipy import stats
from itertools import combinations
from collections import defaultdict as ddict
from mqpipeline.Utils import barplot_annotate_brackets, get_number_rows_cols_for_fig, venn_names, get_number_of_non_na_values, \
    string_similarity_ratio, get_overlap
import logging
import warnings
from adjustText import adjust_text

# TODO VALIDATE descriptive plots not changing between log2 and non log2

# plt.style.use('ggplot')

FIG_FORMAT = ".pdf"
# Venn diagram settings
# TODO figsize
VENN_TITLE_FONT_SIZE = 20
VENN_SET_LABEL_FONT_SIZE = 16
VENN_SUBSET_LABEL_FONT_SIZE = 14
# descriptive plots settings


class MQPlots(Logger):
    possible_plots = [
        "plot_detection_counts", "plot_number_of_detected_proteins", "plot_intensity_histograms",
        "plot_relative_std", "plot_rank", "plot_pathway_analysis", "plot_pathway_timeline", "plot_pathway_proportions",
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
                self.logger.warning("Permission error in function %s. Did you forget to close the file?", f.split(" ")[1])
                return
        return wrapper

    def __init__(
        self, start_dir, replicates, experiment_groups, configs,
        df_protein_names, df_peptide_names,
        interesting_proteins, go_analysis_gene_names,
        loglevel=logging.DEBUG
    ):
        super().__init__(self.__class__.__name__, loglevel=loglevel)
        # property variables
        self._replicate_representation = None
        self._min_number_replicates = None
        self._max_number_replicates = None
        self._replicates_representation = None
        # general information
        self.start_dir = start_dir
        self.replicates = replicates
        self.experiment_groups = experiment_groups
        self.configs = configs
        self.int_mapping = {
            "raw": "Intensity ", "raw_log2": "Intensity ",
            "lfq": "LFQ intensity ", "lfq_log2": "LFQ intensity ",
            "ibaq": "iBAQ ", "ibaq_log2": "iBAQ "
        }
        self.intensity_label_names = {
            "raw": "Intensity", "raw_log2": r"$Log_2$ Intensity",
            "lfq": "LFQ intensity", "lfq_log2": r"$Log_2$ LFQ intensity",
            "ibaq": "iBAQ intensity", "ibaq_log2": r"$Log_2$ iBAQ intensity",
        }
        # data frames
        self.df_protein_names = df_protein_names.set_index(df_protein_names["Gene name fasta"], drop=False)
        if any((col.startswith("LFQ intensity") for col in self.df_protein_names)):
            self.has_lfq = True
        else:
            self.has_lfq = False
            self.logger.warning("LFQ intensities were not found. Raw intensities are used for all plots")
        if any((col.startswith("iBAQ") for col in self.df_protein_names)):
            self.has_ibaq = True
        else:
            self.has_ibaq = False
            self.logger.warning("iBAQ intensities were not found. Raw intensities are used for all plots")
        self.df_peptide_names = df_peptide_names
        # dicts
        self.interesting_proteins = interesting_proteins
        self.go_analysis_gene_names = go_analysis_gene_names

        # extract all raw intensites from the dataframe
        self.all_intensities_raw = self.df_protein_names[
            [f"Intensity {rep}" for exp in self.replicates for rep in self.replicates[exp]]
        ]
        # filter all rows where all intensities are 0
        mask = (self.all_intensities_raw != 0).sum(axis=1) != 0
        self.all_intensities_raw = self.all_intensities_raw[mask]

        # write all intensities of every experiment to one dict entry
        self.intensities_per_experiment_raw = {
            exp: self.all_intensities_raw[
                [f"Intensity {rep}" for rep in self.replicates[exp]]
            ] for exp in self.replicates
        }

        # add log2 intensities
        self.all_intensities_raw_log2 = np.log2(self.all_intensities_raw.replace({0: np.nan}))
        self.intensities_per_experiment_raw_log2 = {
            experiment: np.log2(df.replace({0: np.nan}))
            for experiment, df in self.intensities_per_experiment_raw.items()
        }

        if self.has_lfq:
            # extract all lfq intensities from the dataframe
            self.all_intensities_lfq = self.df_protein_names[
                [f"LFQ intensity {rep}" for exp in self.replicates for rep in self.replicates[exp]]
            ]
            # filter all rows where all intensities are 0
            mask = (self.all_intensities_lfq != 0).sum(axis=1) != 0
            self.all_intensities_lfq = self.all_intensities_lfq[mask]

            # write all intensities of every experiment to one dict entry
            self.intensities_per_experiment_lfq = {
                exp: self.all_intensities_lfq[
                    [f"LFQ intensity {rep}" for rep in self.replicates[exp]]
                ] for exp in self.replicates
            }

            # add log2 values
            self.all_intensities_lfq_log2 = np.log2(self.all_intensities_lfq.replace({0: np.nan}))
            self.intensities_per_experiment_lfq_log2 = {
                experiment: np.log2(df.replace({0: np.nan}))
                for experiment, df in self.intensities_per_experiment_lfq.items()
            }

        else:
            self.all_intensities_lfq = self.all_intensities_raw
            self.intensities_per_experiment_lfq = self.intensities_per_experiment_raw
            self.all_intensities_lfq_log2 = self.all_intensities_raw_log2
            self.intensities_per_experiment_lfq_log2 = self.intensities_per_experiment_raw_log2

        if self.has_ibaq:
            # extract all lfq intensities from the dataframe
            self.all_intensities_ibaq = self.df_protein_names[
                [f"iBAQ {rep}" for exp in self.replicates for rep in self.replicates[exp]]
            ]
            # filter all rows where all intensities are 0
            mask = (self.all_intensities_ibaq != 0).sum(axis=1) != 0
            self.all_intensities_ibaq = self.all_intensities_ibaq[mask]

            # write all intensites of every experiment to one dict entry
            self.intensities_per_experiment_ibaq = {
                exp: self.all_intensities_ibaq[
                    [f"iBAQ {rep}" for rep in self.replicates[exp]]
                ] for exp in self.replicates
            }

            # add log2 values
            self.all_intensities_ibaq_log2 = np.log2(self.all_intensities_ibaq.replace({0: np.nan}))
            self.intensities_per_experiment_ibaq_log2 = {
                experiment: np.log2(df.replace({0: np.nan}))
                for experiment, df in self.intensities_per_experiment_ibaq.items()
            }

        else:
            self.all_intensities_ibaq = self.all_intensities_raw
            self.intensities_per_experiment_ibaq = self.intensities_per_experiment_raw
            self.all_intensities_ibaq_log2 = self.all_intensities_raw_log2
            self.intensities_per_experiment_ibaq_log2 = self.intensities_per_experiment_raw_log2

        # create a dict matching the raw and lfq dfs by string
        self.all_intensities_dict = {
            "lfq": self.all_intensities_lfq, "lfq_log2": self.all_intensities_lfq_log2,
            "raw": self.all_intensities_raw, "raw_log2": self.all_intensities_raw_log2,
            "ibaq": self.all_intensities_ibaq, "ibaq_log2": self.all_intensities_ibaq_log2
        }
        self.intensities_per_experiment_dict = {
            "lfq": self.intensities_per_experiment_lfq, "lfq_log2": self.intensities_per_experiment_lfq_log2,
            "raw": self.intensities_per_experiment_raw, "raw_log2": self.intensities_per_experiment_raw_log2,
            "ibaq": self.intensities_per_experiment_ibaq, "ibaq_log2": self.intensities_per_experiment_ibaq_log2
        }

        # create all sets that are required for plotting
        # this could be turned into a property
        venn_int = self.configs.get("plot_venn_results_intensity", "raw")
        self.protein_ids = ddict(dict)
        self.whole_experiment_protein_ids = {}
        for experiment in self.replicates:
            exp_prot_ids = set()
            for rep in self.replicates[experiment]:
                series = self.intensities_per_experiment_dict[venn_int][experiment].loc[:, self.int_mapping[venn_int] + rep]
                idx = series.loc[series > 0].index
                rep_set = set(self.df_protein_names.loc[idx, "Protein name"])
                self.protein_ids[experiment][rep] = rep_set
                exp_prot_ids |= rep_set
            self.whole_experiment_protein_ids[experiment] = exp_prot_ids

        # set all result dirs
        # TODO: for now just create all of them
        # path for venn diagrams
        self.file_dir_venn = os.path.join(self.start_dir, "venn")
        # create file structure and folder
        os.makedirs(self.file_dir_venn, exist_ok=True)
        # path for descriptive plots
        self.file_dir_descriptive = os.path.join(self.start_dir, "descriptive")
        # create file structure and folder
        os.makedirs(self.file_dir_descriptive, exist_ok=True)
        self.file_dir_go_analysis = os.path.join(self.start_dir, "go_analysis")
        # create file structure and folder
        os.makedirs(self.file_dir_go_analysis, exist_ok=True)
        self.file_dir_volcano = os.path.join(self.start_dir, "volcano")
        # create file structure and folder
        os.makedirs(self.file_dir_volcano, exist_ok=True)

    @classmethod
    def from_MQInitializer(cls, mqinti_instance):
        return cls(
            start_dir = mqinti_instance.start_dir,
            replicates = mqinti_instance.replicates,
            experiment_groups=mqinti_instance.experiment_groups,
            configs = mqinti_instance.configs,
            df_protein_names = mqinti_instance.df_protein_names,
            df_peptide_names = mqinti_instance.df_peptide_names,
            interesting_proteins = mqinti_instance.interesting_proteins,
            go_analysis_gene_names = mqinti_instance.go_analysis_gene_names,
            loglevel = mqinti_instance.logger.getEffectiveLevel()
            )

    @property
    def replicate_representation(self):
        if self._replicate_representation is None:
            d = {}
            for experiment in self.replicates:
                if len(self.replicates[experiment]) > 1:
                    experiment_rep = experiment.rstrip("0123456789_")
                    experiment_rep = " ".join(experiment_rep.split("_"))
                    d[experiment] = experiment_rep
                else:
                    d[experiment] = experiment
            self._replicate_representation = d
        return self._replicate_representation

    @property
    def min_number_replicates(self):
        if self._min_number_replicates is None:
            self._min_number_replicates = min([len(self.replicates[experiment]) for experiment in self.replicates])
        return self._min_number_replicates

    @property
    def max_number_replicates(self):
        if self._max_number_replicates is None:
            self._max_number_replicates = max([len(self.replicates[experiment]) for experiment in self.replicates])
        return self._max_number_replicates

    def create_results(self):
        for plot_name in MQPlots.possible_plots:
            intensity_name = plot_name + "_intensity"
            if self.configs.get(plot_name, False):
                self.logger.debug(f"creating plot {plot_name}")
                getattr(self, plot_name)(self.configs.get(intensity_name, "raw"))
        self.logger.info("Done creating plots")

    @exception_handler
    def save_venn(self, ex: str, sets, set_names):
        creates_figure = True
        plt.close("all")
        plt.figure(figsize=(14, 7))
        title = self.replicate_representation.get(ex, ex)
        plt.title(title, fontsize=VENN_TITLE_FONT_SIZE)

        # create venn diagram based on size of set
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
    def save_bar_venn(self, ex: str, named_sets):
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
        name = self.replicate_representation.get(ex, ex)
        fig.suptitle(name, fontsize=20)
        # create the bar plot
        ax1.bar(x, heights)
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

    @exception_handler
    def plot_venn_results(self, *args):
        # TODO the *args should not be used?
        self.logger.info("Creating venn diagrams")

        # create venn diagrams comparing all replicates within an experiment
        for ex in self.protein_ids:
            self.logger.debug("Creating venn diagram for experiment %s", ex)
            set_names = self.protein_ids[ex].keys()
            sets = self.protein_ids[ex].values()
            # save the resulting venn diagram
            self.save_venn(ex, sets, set_names)
            # save the sets to txt files
            self.save_venn_names(self.protein_ids[ex])
            # create a mixture of bar and venn diagram
            self.save_bar_venn(ex, self.protein_ids[ex])
        # create venn diagrams comparing all experiments
        # first compare only proteins which are found between all replicates
        self.logger.debug("Creating venn diagram for intersection")
        experiment_intersection_sets = {
            exp: set.intersection(*(self.protein_ids[exp][rep] for rep in self.protein_ids[exp]))
            for exp in self.protein_ids
        }
        self.save_bar_venn("All experiments intersection", experiment_intersection_sets)
        # then compare all proteins that are found at all
        self.logger.debug("Creating venn diagram for union")
        experiment_union_sets = {
            exp: set.union(*(self.protein_ids[exp][rep] for rep in self.protein_ids[exp]))
            for exp in self.protein_ids
        }
        self.save_bar_venn("All experiments union", experiment_union_sets)

    @exception_handler
    def plot_venn_groups(self, *args):
        if self.experiment_groups:
            if len(self.experiment_groups) <= 3:
                self.logger.debug("Creating venn diagram for group comparison intersection")
                group_proteins = [
                    set.union(*[
                        set.intersection(*(self.protein_ids[experiment][rep] for rep in self.protein_ids[experiment]))
                        for experiment in self.experiment_groups[group]
                    ])
                    for group in self.experiment_groups
                ]
                self.save_venn("group_comparison_intersection", group_proteins, [self.experiment_groups[group][0] for group in self.experiment_groups])
                self.logger.debug("Creating venn diagram for group comparison union")
                group_proteins = [
                    set.union(*(self.protein_ids[experiment][rep]
                                for experiment in self.experiment_groups[group] for rep in self.protein_ids[experiment]))
                    for group in self.experiment_groups
                ]
                self.save_venn("group_comparison_union", group_proteins, [self.experiment_groups[group][0] for group in self.experiment_groups])
            else:
                self.logger.debug("More than 3 groups for venn group plot. Creating pairwise comparison instead")
                for g1, g2 in combinations(self.experiment_groups, 2):
                    self.logger.debug("Creating venn diagram for group comparison intersection")
                    group_proteins = [
                        set.union(*[
                            set.intersection(
                                *(self.protein_ids[experiment][rep] for rep in self.protein_ids[experiment]))
                            for experiment in self.experiment_groups[group]
                        ])
                        for group in [g1, g2]
                    ]
                    self.save_venn(f"{g1}_vs_{g2}_intersection", group_proteins, [g1, g2])
                    self.save_venn_names({g1 + "_intersection": group_proteins[0], g2 + "_intersection": group_proteins[1]})
                    self.logger.debug("Creating venn diagram for group comparison union")
                    group_proteins = [
                        set.union(*(self.protein_ids[experiment][rep]
                                    for experiment in self.experiment_groups[group]
                                    for rep in self.protein_ids[experiment]))
                        for group in [g1, g2]
                    ]
                    self.save_venn(f"{g1}_vs_{g2}_union", group_proteins, [g1, g2])
                    self.save_venn_names({g1 + "_union": group_proteins[0], g2 + "_union": group_proteins[1]})

        else:
            self.logger.warning("Skipping venn group plot because grouping failed")

    @exception_handler
    def plot_detection_counts(self, df_to_use: str = "raw"):
        plt.close("all")
        # determine number of rows and columns in the plot based on the number of experiments
        n_rows_experiment, n_cols_experiment = get_number_rows_cols_for_fig(self.replicates)
        fig, axarr = plt.subplots(n_rows_experiment, n_cols_experiment, sharex=True, squeeze=True,
                                  figsize=(4 * n_rows_experiment, 4 * n_cols_experiment))
        fig.suptitle(f"Detection counts from {self.intensity_label_names[df_to_use]} intensities")
        for experiment, ax in zip(self.replicates, axarr.flat):
            intensities = self.intensities_per_experiment_dict[df_to_use][experiment]
            # from 0 to number of replicates, how often was each protein detected
            counts = (intensities > 0).sum(axis=1)
            counts = counts[counts > 0]
            heights = [len(counts[counts == value]) for value in sorted(counts.unique())]
            y_pos = [value for value in sorted(counts.unique())]
            max_val = max(heights)

            ax.set_title(f"{experiment}, total detected: {len(counts)}")
            ax.barh(y_pos, heights)
            for y, value in zip(y_pos, heights):
                ax.text(max_val / 2, y, value,
                        verticalalignment='center', horizontalalignment='center')
            ax.set_yticks(y_pos)
            ax.set_xlabel("Counts")

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        res_path = os.path.join(self.file_dir_descriptive, f"detected_counts" + FIG_FORMAT)
        fig.savefig(res_path, dpi=200, bbox_inches="tight")

    @exception_handler
    def plot_number_of_detected_proteins(self, df_to_use: str = "raw"):
        plt.close("all")
        # determine number of rows and columns in the plot based on the number of experiments
        n_rows_experiment, n_cols_experiment = get_number_rows_cols_for_fig(self.replicates)
        fig, axarr = plt.subplots(n_rows_experiment, n_cols_experiment,
                                  figsize=(5 * n_cols_experiment, 3 * n_rows_experiment))
        fig.suptitle(f"Number of detected proteins from {self.intensity_label_names[df_to_use]}")

        for experiment, ax in zip(self.replicates, axarr.flat):
            intensities = self.intensities_per_experiment_dict[df_to_use][experiment]

            # how many proteins were detected per replicate and in total
            counts = (intensities > 0).sum(axis=1)
            counts = counts[counts > 0]
            heights = [len(counts), 0]
            # labels start at 0 so we prepend one empty string
            labels = ["", "Total", ""]
            for col in intensities:
                h = len(intensities[col][intensities[col] > 0])
                heights.append(h)
                labels.append(col.replace(self.int_mapping[df_to_use], ""))
            mean_height = np.mean(heights[2:])
            # self.logger.debug(labels)
            ax.barh([x for x in range(len(heights))], heights)
            ax.text(mean_height * 0.99, 1, f"{mean_height:.2f}", horizontalalignment='right')
            ax.set_title(f"{experiment}")
            ax.axvline(mean_height, linestyle="--", color="black", alpha=0.6)
            ax.set_yticklabels(labels)
            ax.set_xlabel("Counts")
            #ax1.tick_params(rotation=70)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        res_path = os.path.join(self.file_dir_descriptive, "detection_per_replicate" + FIG_FORMAT)
        fig.savefig(res_path, dpi=200, bbox_inches="tight")

    @exception_handler
    def plot_intensity_histograms(self, df_to_use: str = "raw"):
        plt.close("all")
        n_rows_replicates, n_cols_replicates = get_number_rows_cols_for_fig([x for l in self.replicates.values() for x in l])
        # make a intensity histogram for every replicate
        fig, axarr = plt.subplots(n_rows_replicates, n_cols_replicates,
                                  figsize=(5 * n_cols_replicates, 5 * n_rows_replicates))
        fig_index = 0
        fig.suptitle(f"Replicate {self.intensity_label_names[df_to_use]} histograms")
        for experiment in self.replicates:
            intensities = self.intensities_per_experiment_dict[df_to_use][experiment]
            for col in intensities:
                mask = intensities[col] > 0
                ax = axarr.flat[fig_index]
                ax.set_title(f"{col}".replace(self.int_mapping[df_to_use], ""))
                ax.hist(intensities[col][mask],
                        bins=np.logspace(np.log10(intensities[col][mask].min()),
                                         np.log10(intensities[col][mask].max()),
                                         25)
                        )
                ax.set_xscale("log")
                ax.set_xlabel(f"{self.intensity_label_names[df_to_use]}")
                ax.set_ylabel("Counts")
                fig_index += 1

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        res_path = os.path.join(self.file_dir_descriptive, "intensity_histograms" + FIG_FORMAT)
        fig.savefig(res_path, dpi=200, bbox_inches="tight")

    @exception_handler
    def plot_scatter_replicates(self, df_to_use: str = "raw"):
        # TODO where to put uniquely found proteins
        for experiment in self.replicates:
            plt.close("all")
            # scatter plot of all replicates vs all replicates
            fig, ax = plt.subplots(1, 1, figsize=(7, 7))
            for rep1, rep2 in combinations(self.replicates[experiment], 2):
                x1 = self.all_intensities_dict[df_to_use][self.int_mapping[df_to_use] + rep1]
                x2 = self.all_intensities_dict[df_to_use][self.int_mapping[df_to_use] + rep2]
                mask = np.logical_or(x1.fillna(0) != 0, x2.fillna(0) != 0)
                corr_mask = np.logical_and(x1.fillna(0) != 0, x2.fillna(0) != 0)
                exp = r"$r^{2}$"
                ax.scatter(x1[mask] + 1e2, x2[mask] + 1e2, label=f"{rep1} vs {rep2}, "
                           fr"{exp}: {stats.pearsonr(x1[corr_mask], x2[corr_mask])[0] ** 2:.4f}",
                           alpha=0.5, marker=".")
                if "log2" in df_to_use:
                    pass
                    #ax.set_xscale('log', basex=2)
                    #ax.set_yscale('log', basey=2)
                ax.set_xlabel(f"{self.intensity_label_names[df_to_use]} intensity")
                ax.set_ylabel(f"{self.intensity_label_names[df_to_use]} intensity")

            fig.legend()

            res_path = os.path.join(self.file_dir_descriptive,
                                    f"{self.replicate_representation[experiment].replace(' ', '_')}_scatter" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")

    @exception_handler
    def plot_rank(self, df_to_use: str = "raw"):
        all_pathway_proteins = set.union(*(set(x) for x in self.interesting_proteins.values()))
        for experiment in self.replicates:
            plt.close("all")
            intensities = self.intensities_per_experiment_dict[df_to_use][experiment]
            # protein ranks vs intensity
            # calculate mean intensity for the experiment and sort from highest to lowest
            m_intensity = intensities.mean(axis=1).sort_values(ascending=False)
            # filter all 0
            m_intensity = m_intensity[m_intensity > 0]
            # create dict to map each protein its respective rank and mean intensity
            dic = {idx: (i, value) for i, (idx, value) in enumerate(m_intensity.items())}

            found_proteins = set(m_intensity.index)
            # get all proteins that are not part of any pathway
            non_pathway_proteins = found_proteins - all_pathway_proteins
            # get all proteins that are part of any pathway
            pathway_proteins = found_proteins & all_pathway_proteins
            rank_identified_proteins = [dic[protein][0] for protein in pathway_proteins]
            # only if more than 0 proteins are identified
            if rank_identified_proteins:
                median_pathway_rank = int(np.median(rank_identified_proteins))
                median_intensity = m_intensity.iloc[median_pathway_rank]
            # plot the non pathway proteins
            x = [dic[protein][0] for protein in non_pathway_proteins]
            y = [dic[protein][1] for protein in non_pathway_proteins]

            fig, ax = plt.subplots(1, 1, figsize=(14, 7))
            if "log2" in df_to_use:
                pass
                #ax.set_yscale("log", basey=2)
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
                xmin, xmax = ax.get_xbound()
                xm = (median_pathway_rank + abs(xmin)) / (abs(xmax) + abs(xmin))
                ymin, ymax = ax.get_ybound()
                ym = ((median_intensity) - (ymin) ) / ((ymax) - (ymin))
                # plot the median rank and intensity at that rank
                ax.axvline(median_pathway_rank, ymax=ym, linestyle="--", color="black", alpha=0.6)
                ax.axhline(median_intensity, xmax=xm, linestyle="--", color="black", alpha=0.6)
                ax.text(xmin * 0.9, median_intensity * 0.9,
                        f"median rank: {median_pathway_rank} ({median_pathway_rank/len(m_intensity) * 100 :.1f}%) "
                        f"with intensity: {median_intensity:.2E}",
                        verticalalignment="top", horizontalalignment="left")
                # ax.annotate(f"median rank: {median_pathway_rank} with intensity: {median_intensity}",
                #             xy=(median_pathway_rank, median_intensity), xycoords='data',
                #             xytext=(0.1, 0.1), textcoords='axes fraction',
                #             arrowprops=dict(facecolor='black', shrink=0.05),
                #             horizontalalignment='right', verticalalignment='top',
                #             )

            ax.set_xlabel("Protein rank")
            ax.set_ylabel(f"{self.intensity_label_names[df_to_use]} mean")
            fig.legend()

            res_path = os.path.join(self.file_dir_descriptive,
                                    f"{self.replicate_representation[experiment].replace(' ', '_')}_rank" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")

    @exception_handler
    def plot_relative_std(self, df_to_use: str = "raw"):
        plt.close("all")
        for experiment in self.replicates:
            intensities = self.intensities_per_experiment_dict[df_to_use][experiment]
            y = intensities.copy()
            # below cutoff remove zeros, above threshold include 0
            cutoff = 1e6 if "log2" not in df_to_use else np.log2(1e6)
            index = y.replace(0, np.nan).mean(axis=1) < cutoff
            y[index] = y[index].replace(0, np.nan)
            relative_std_percent = y.std(axis=1) / y.mean(axis=1) * 100

            bins = np.array([10, 20, 30])
            inds = np.digitize(relative_std_percent, bins).astype(int)

            cm = {0: "navy", 1: "royalblue", 2: "skyblue", 3: "darkgray"}
            colors = pd.Series([cm.get(x, "black") for x in inds], index=relative_std_percent.index)
            color_counts = {color: (colors == color).sum() for color in colors.unique()}
            mask = ~relative_std_percent.isna()

            # intensity vs relative standard deviation
            # TODO show only experiments with full amount of intensities found
            fig, ax = plt.subplots(1, 1, figsize=(14, 7))
            if "log2" in df_to_use:
                pass
                #ax.set_xscale("log", basex=2)
            ax.scatter(y.mean(axis=1)[mask], relative_std_percent[mask], c=colors[mask], marker="o", s=(2* 72./fig.dpi)**2, alpha=0.8)
            ax.set_xlabel(f"Mean {self.intensity_label_names[df_to_use]}")
            ax.set_ylabel("Relative Standard deviation [%]")
            ax.axvline(cutoff, color="black", alpha=0.5)
            xmin, xmax = ax.get_xbound()
            cumulative_count = 0
            for i, bin_ in enumerate(bins):
                cumulative_count += color_counts.get(cm[i], 0)
                ax.axhline(bin_, color=cm[i])
                ax.text(xmin, bin_, cumulative_count)

            res_path = os.path.join(self.file_dir_descriptive,
                                    f"{self.replicate_representation[experiment].replace(' ', '_')}_rel_std" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")
            plt.close(fig)

    @exception_handler
    def plot_pathway_analysis(self, df_to_use: str = "raw"):
        ex_list = list(self.replicates.keys())
        plot_ns = False
        threshhold = 0.05
        for pathway in self.interesting_proteins:
            plt.close("all")
            found_proteins = set(self.interesting_proteins[pathway])
            found_proteins &= set(self.all_intensities_dict[df_to_use].index)
            found_proteins = sorted(list(found_proteins))
            if len(found_proteins) < 1:
                self.logger.warning("Skipping pathway %s in pathway analysis because no proteins were found", pathway)
                continue
            n_rows, n_cols = get_number_rows_cols_for_fig(found_proteins)
            fig, axarr = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * int(len(self.replicates) / 3)))
            fig.suptitle(pathway)
            all_heights = {}
            try:
                axiterator = axarr.flat
            except AttributeError:
                axiterator = [axarr]
            for protein, ax in zip(found_proteins, axiterator):
                ax.set_title(protein)
                if "log2" in df_to_use:
                    pass
                    #ax.set_yscale("log", basex=2)
                heights = []
                for idx, experiment in enumerate(sorted(self.replicates)):
                    protein_intensities = self.intensities_per_experiment_dict[df_to_use][experiment].loc[protein, :]
                    ax.scatter(protein_intensities, [idx] * len(protein_intensities), label=f"{experiment}")
                    # self.logger.debug(f"{protein}: min: {min_i}, max: {max_i}")
                    heights.append(np.max(protein_intensities))
                all_heights[protein] = heights
                ax.set_ylim((-1, len(self.replicates) + 1))

                ax.set_yticks([i for i in range(len(self.replicates))])
                ax.set_yticklabels(sorted(self.replicates))
            #handles, labels = axiterator[0].get_legend_handles_labels()
            #fig.legend(handles, labels, bbox_to_anchor=(1.04, 0.5), loc="center left")
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            res_path = os.path.join(self.file_dir_descriptive, f"pathway_analysis_{pathway}_no_labels" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")
            try:
                axiterator = axarr.flat
            except AttributeError:
                axiterator = [axarr]
            significances = []
            for protein, ax in zip(found_proteins, axiterator):
                n_annotations = 0
                for e1, e2 in combinations(self.replicates, 2):
                    v1 = self.intensities_per_experiment_dict[df_to_use][e1].loc[protein, :].astype(float)
                    v2 = self.intensities_per_experiment_dict[df_to_use][e2].loc[protein, :].astype(float)
                    # only perform the test if more than 3 data points are available for both experiments
                    if len(v1[v1 > 0]) < 3 or len(v2[v2 > 0]) < 3:
                        # self.logger.debug(f"skipping test for: {protein}, {e1} vs {e2}")
                        continue
                    test = stats.ttest_ind(v1, v2, equal_var=self.configs.get("equal_var", False))
                    if plot_ns or test[1] < threshhold:
                        # barplot_annotate_brackets(ax, ex_list.index(e1), ex_list.index(e2), test[1], range(len(ex_list)),
                        #                          all_heights[protein], dh=0.05 + 0.1 * n_annotations)
                        n_annotations += 1
                        significances.append((protein, e1, e2, test[1]))
            df = pd.DataFrame(significances, columns=["protein", "experiment1", "experiment2", "pvalue"])
            df.to_csv(os.path.join(self.file_dir_descriptive, f"pathway_analysis_{pathway}_table.csv"), index=False)

            #handles, labels = axiterator[0].get_legend_handles_labels()
            #fig.legend(handles, labels, bbox_to_anchor=(1.04, 0.5), loc="center left")
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            res_path = os.path.join(self.file_dir_descriptive, f"pathway_analysis_{pathway}" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")

    @exception_handler
    def plot_pathway_timeline(self, df_to_use: str = "raw"):
        group_colors = {
            "SD": "#808080",
            "4W": "#0b8040",
            "6W": "#ec2024",
            "8W": "#4378bb"
        }
        groups = {k: "SD" if "SD" in k else k.split("_")[0] for k in self.replicates}
        x_values = {k: sum([int(s.replace("W", "")) for s in k.split("_") if s.endswith("W")]) for k in self.replicates}
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
                if "log2" in df_to_use:
                    pass
                    #ax.set_yscale("log", basey=2)
                for idx, experiment in enumerate(self.replicates):
                    protein_intensities = self.intensities_per_experiment_dict[df_to_use][experiment].loc[protein, :]
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
                ax.set_ylim(bottom=protein_minimum, top=protein_maximum)
                ax.set_xlim(left=0, right=max_time + 1)
            handles, labels = axiterator[0].get_legend_handles_labels()
            fig.legend(handles, labels, bbox_to_anchor=(1.04, 0.5), loc="center left")
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            res_path = os.path.join(self.file_dir_descriptive, f"pathway_timeline_{pathway}" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")

    @exception_handler
    def plot_experiment_comparison(self, df_to_use: str = "raw"):
        for ex1, ex2 in combinations(self.replicates, 2):
            plt.close("all")
            protein_intensities_ex1 = self.intensities_per_experiment_dict[df_to_use][ex1]
            counts_ex1 = (protein_intensities_ex1 > 0).sum(axis=1) == len(protein_intensities_ex1.columns)
            protein_intensities_ex2 = self.intensities_per_experiment_dict[df_to_use][ex2]
            counts_ex2 = (protein_intensities_ex2 > 0).sum(axis=1) == len(protein_intensities_ex2.columns)
            # combine intersections
            intersection = pd.concat(
                [protein_intensities_ex1.mean(axis=1).rename("ex1"), protein_intensities_ex2.mean(axis=1).rename("ex2")]
                , axis=1)
            intersection = intersection[counts_ex1 & counts_ex2]
            fig6, ax6 = plt.subplots(1, 1, figsize=(7, 7))
            if "log2" in df_to_use:
                pass
                #ax6.set_xscale('log', basex=2)
                #ax6.set_yscale('log', basey=2)
            # plot intersection
            ax6.scatter(intersection["ex1"], intersection["ex2"], s=8, alpha=0.6, marker=".")
            xmin, xmax = ax6.get_xbound()
            ymin, ymax = ax6.get_ybound()
            ax6.set_xlim(min(xmin, ymin), max(xmax, ymax))
            ax6.set_ylim(min(xmin, ymin), max(xmax, ymax))
            ax6.set_xlabel(ex1)
            ax6.set_ylabel(ex2)

            # TODO add r2
            res_path = os.path.join(self.file_dir_descriptive,
                                    f"{self.replicate_representation[ex1].replace(' ', '_')}_vs_{self.replicate_representation[ex2].replace(' ', '_')}_intersection" + FIG_FORMAT)
            fig6.savefig(res_path, dpi=200, bbox_inches="tight")
            plt.close("all")
            # combine total
            protein_intensities_ex1 = protein_intensities_ex1.mean(axis=1).rename("ex1")
            protein_intensities_ex2 = protein_intensities_ex2.mean(axis=1).rename("ex2")
            conc = pd.concat([protein_intensities_ex1, protein_intensities_ex2], axis=1)
            # TODO unique proteins are missing
            conc = conc[conc > 0]
            # plot total
            plt.close("all")
            fig9, ax9 = plt.subplots(1, 1, figsize=(7, 7))
            ax9.scatter(conc["ex1"], conc["ex2"], s=8, alpha=0.6, marker=".")
            #ax9.set_xscale("log")
            #ax9.set_yscale("log")
            xmin, xmax = ax9.get_xbound()
            ymin, ymax = ax9.get_ybound()
            ax9.set_xlim(min(xmin, ymin), max(xmax, ymax))
            ax9.set_ylim(min(xmin, ymin), max(xmax, ymax))
            ax9.set_xlabel(ex1)
            ax9.set_ylabel(ex2)
            # TODO add r2
            res_path = os.path.join(self.file_dir_descriptive,
                                    f"{self.replicate_representation[ex1].replace(' ', '_')}_vs_{self.replicate_representation[ex2].replace(' ', '_')}_total" + FIG_FORMAT)
            plt.savefig(res_path, dpi=200, bbox_inches="tight")

    @exception_handler
    def plot_pathway_proportions(self, df_to_use: str = "raw"):
        plt.close("all")
        experiment_proportion = {
            experiment:
                {
                    pathway: len(set(pproteins) & set(intensities[intensities.mean(axis=1) > 0].index))
                    for pathway, pproteins in self.interesting_proteins.items()
                }
            for experiment, intensities in self.intensities_per_experiment_dict[df_to_use].items()
        }

        df_total_counts = pd.DataFrame(experiment_proportion).T
        df_proportion = df_total_counts / df_total_counts.sum()
        bottom_cumsum = df_proportion.cumsum()

        fig, ax = plt.subplots(1, 1, figsize=(14, 7))
        fig.suptitle("Number of identified proteins per pathway")

        for i, index in enumerate(df_proportion.index):
            if i == 0:
                bottom = [0] * len(df_proportion.loc[index])
            else:
                bottom = bottom_cumsum.iloc[i - 1]
            ax.bar(range(len(df_proportion.loc[index])), df_proportion.loc[index],
                    bottom=bottom, label=df_proportion.index[i], tick_label=df_proportion.columns)
            for j, text_index in enumerate(df_total_counts.columns):
                height = (bottom_cumsum.iloc[i, j] + bottom[j]) / 2
                ax.text(j, height, df_total_counts.iloc[i, j], horizontalalignment="center")
        # fig.xticks(rotation=45)
        fig.legend()
        res_path = os.path.join(self.file_dir_descriptive, "pathway_proportions" + FIG_FORMAT)
        fig.savefig(res_path, dpi=200, bbox_inches="tight")

    def create_report(self):
        pass

    @exception_handler
    def plot_go_analysis(self, df_to_use: str = "raw"):
        self.logger.info("Creating go analysis plots")
        plt.close("all")

        # all proteins that were detected in any replicate
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
            for experiment in self.replicates:
                # create df with intensity means for specific experiment over all replicates
                mean_intensity = self.intensities_per_experiment_dict[df_to_use][experiment].mean(axis=1)
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

        # TODO only show pvalue when sigificant
        # TODO also create table
        # TODO move labels to bars, remove legend
        fig, ax = plt.subplots(1, 1, figsize=(7, int(len(self.replicates) * len(self.go_analysis_gene_names) / 3)))

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
        res_path = os.path.join(self.file_dir_go_analysis, "go_analysis" + FIG_FORMAT)
        fig.savefig(res_path, dpi=200, bbox_inches="tight")

    @exception_handler
    def plot_r_volcano(self, df_to_use: str = "lfq_log2"):
        # import r interface package
        from rpy2.robjects.packages import importr
        from rpy2.robjects import pandas2ri
        # allow conversion of pd objects to r
        pandas2ri.activate()

        # install r packages
        from mqpipeline.Utils import install_r_dependencies
        r_package_names = ("BiocManager", )
        r_bioconducter_package_names = ("limma", )
        install_r_dependencies(r_package_names, r_bioconducter_package_names)

        # import r packages
        limma = importr("limma")

        for g1, g2 in combinations(self.experiment_groups, 2):
            # get groups based on name
            # calculate the mean of the technical replicates as a proxy for the biological replicate
            # then use all biological replicates in a group for the limma
            v1 = {
                exp:
                    self.all_intensities_dict[df_to_use].loc[
                        :, [self.int_mapping[df_to_use] + rep for rep in self.replicates[exp]]
                    ].mean(axis=1)
                for exp in self.experiment_groups[g1]
            }
            v1 = pd.concat(v1, axis=1)
            v2 = {
                exp:
                    self.all_intensities_dict[df_to_use].loc[:, [self.int_mapping[df_to_use] + rep for rep in self.replicates[exp]]].mean(axis=1)
                    for exp in self.experiment_groups[g2]
            }
            v2 = pd.concat(v2, axis=1)
            #v1 = self.all_intensities_dict[df_to_use].loc[
            #    :, [self.int_mapping[df_to_use] + rep for exp in self.experiment_groups[g1] for rep in self.replicates[exp]]]
            #v2 = self.all_intensities_dict[df_to_use].loc[
            #    :, [self.int_mapping[df_to_use] + rep for exp in self.experiment_groups[g2] for rep in self.replicates[exp]]]
            # filter entries with too many nans based on function
            non_na_group_1 = get_number_of_non_na_values(v1.shape[1])
            non_na_group_2 = get_number_of_non_na_values(v2.shape[1])
            mask_1 = (v1 > 0).sum(axis=1) >= non_na_group_1
            mask_2 = (v2 > 0).sum(axis=1) >= non_na_group_2
            mask = np.logical_and(mask_1, mask_2)
            # determine missing
            missing_1 = (v1 > 0).sum(axis=1) == 0
            missing_2 = (v2 > 0).sum(axis=1) == 0
            # determine exclusive
            exclusive_1 = np.logical_and(mask_1, missing_2)
            exclusive_2 = np.logical_and(mask_2, missing_1)

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
            c_matrix = limma.makeContrasts(f"{g2}-{g1}", levels=r_design)
            contrast_fit = limma.contrasts_fit(fit, c_matrix)
            fit_bayes = limma.eBayes(contrast_fit)
            res = limma.topTable(fit_bayes, adjust="BH", number=df.shape[0])
            # transform back to python
            with warnings.catch_warnings():
                # catch a warning from ri2py where DataFrame.from_items is being used
                warnings.simplefilter("ignore", FutureWarning)
                ress = pandas2ri.ri2py(res)
            # extract index
            ress.index = pd.Index([x for x in res.rownames], name="Gene_names")
            col = "adj.P.Val"  # P.Value
            plot_data = ress.loc[:, ["logFC", col]]
            # save all values
            plot_data.to_csv(os.path.join(self.file_dir_volcano, f"volcano_plot_data_{g1}_vs_{g2}_full.csv"))
            # save significant values
            plot_data[plot_data[col] < 0.05].to_csv(
                os.path.join(self.file_dir_volcano, f"volcano_plot_data_{g1}_vs_{g2}_significant.csv"))
            # calculate mean intensity for unique genes for plotting
            unique_g1 = v1[exclusive_1].mean(axis=1).rename(f"{df_to_use} mean intensity")
            unique_g2 = v2[exclusive_2].mean(axis=1).rename(f"{df_to_use} mean intensity")
            # save unique values
            unique_g1.to_csv(os.path.join(self.file_dir_volcano, f"volcano_plot_data_{g1}_vs_{g2}_unique_{g1}.csv"), header=True)
            unique_g2.to_csv(os.path.join(self.file_dir_volcano, f"volcano_plot_data_{g1}_vs_{g2}_unique_{g2}.csv"), header=True)

            def get_volcano_color(fchange, pval):
                if pval > 0.05:
                    return "gray"
                elif fchange >= 0:
                    return "red"
                elif fchange < 0:
                    return "blue"
                else:
                    raise ValueError("heisenbug")
            # plot
            fig, ax = plt.subplots(1, 1, figsize=(7, 7))
            ax_unique = ax.twinx()

            # non sign gray, left side significant blue right side red
            color = [get_volcano_color(log_fold_change, p_val)
                     for log_fold_change, p_val in zip(plot_data["logFC"], plot_data[col])]
            ax.scatter(plot_data["logFC"], -np.log10(plot_data[col]), color=color)
            # add line at significance threshold
            ax.axhline(-np.log10(0.05), linestyle="--", color="black", alpha=0.6)
            xmin, xmax = ax.get_xbound()
            abs_max = max(abs(xmin), abs(xmax))
            # plot unique values with mean intensity at over maximum
            ax_unique.scatter([-(abs_max * 1.05)] * len(unique_g1), unique_g1, color="royalblue")
            ax_unique.scatter([abs_max * 1.05] * len(unique_g2), unique_g2, color="lightcoral")
            # figure stuff
            fig.suptitle(f"{g1} vs {g2}")
            ax.set_xlabel(r"$Log_2$ Fold Change")
            ax.set_ylabel(r"-$Log_{10}$" + " p value")
            ax_unique.set_ylabel(self.intensity_label_names[df_to_use])
            res_path = os.path.join(self.file_dir_volcano, f"volcano_{g1}_{g2}_no_annotation" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")
            significant_upregulated = plot_data[(plot_data["logFC"] < 0) & (plot_data[col] < 0.05)].sort_values(by=[col], ascending=True).head(10)
            significant_downregulated = plot_data[(plot_data["logFC"] > 0) & (plot_data[col] < 0.05)].sort_values(by=[col], ascending=True).head(10)
            significant = pd.concat([significant_upregulated, significant_downregulated])
            texts = []
            for log_fold_change, p_val, gene_name in zip(significant["logFC"], significant[col], significant.index):
                texts.append(ax.text(log_fold_change, -np.log10(p_val), gene_name, ha="center", va="center"))
            adjust_text(texts, arrowprops=dict(width=0.2, headwidth=0, color='black'), ax=ax)
            res_path = os.path.join(self.file_dir_volcano, f"volcano_{g1}_{g2}_annotation" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")
            # TODO scatter plot of significant genes
