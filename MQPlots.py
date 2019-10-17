import pandas as pd
import numpy as np
import os
import sys
from Logger import Logger
from matplotlib_venn import venn3, venn2
import matplotlib.pyplot as plt
from scipy import stats
from itertools import combinations
from collections import defaultdict as ddict
from Utils import barplot_annotate_brackets, get_number_rows_cols_for_fig, venn_names
import logging


# plt.style.use('ggplot')

FIG_FORMAT = ".pdf"
# Venn diagram settings
# TODO figsize
VENN_TITLE_FONT_SIZE = 20
VENN_SET_LABEL_FONT_SIZE = 16
VENN_SUBSET_LABEL_FONT_SIZE = 14
# descriptive plots settings

class MQPlots(Logger):
    def __init__(
        self, start_dir, replicates, configs,
        df_protein_names, df_peptide_names,
        interesting_proteins, interesting_receptors, go_analysis_gene_names,
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
        self.configs = configs
        # data frames
        self.df_protein_names = df_protein_names.set_index(df_protein_names["Gene name fasta"], drop=False)
        if any((col.startswith("LFQ intensity") for col in self.df_protein_names)):
            self.has_lfq = True
        else:
            self.has_lfq = False
            self.logger.warning("LFQ intensities were not found. Raw intensities are used for all plots")
        self.df_peptide_names = df_peptide_names
        # dicts
        self.interesting_proteins = interesting_proteins
        self.interesting_receptors = interesting_receptors
        self.go_analysis_gene_names = go_analysis_gene_names

        # extract all raw intensites from the dataframe
        self.all_intensities_raw = self.df_protein_names[
            [f"Intensity {rep}" for exp in self.replicates for rep in self.replicates[exp]]
        ]
        # filter all rows where all intensities are 0
        mask = (self.all_intensities_raw != 0).sum(axis=1) != 0
        self.all_intensities_raw = self.all_intensities_raw[mask]

        # write all intensites of every experiment to one dict entry
        self.intensites_per_experiment_raw = {
            exp: self.all_intensities_raw[
                [f"Intensity {rep}" for rep in self.replicates[exp]]
            ] for exp in self.replicates
        }

        if self.has_lfq:
            # extract all lfq intensites from the dataframe
            self.all_intensities_lfq = self.df_protein_names[
                [f"LFQ intensity {rep}" for exp in self.replicates for rep in self.replicates[exp]]
            ]
            # filter all rows where all intensities are 0
            mask = (self.all_intensities_lfq != 0).sum(axis=1) != 0
            self.all_intensities_lfq = self.all_intensities_lfq[mask]

            # write all intensites of every experiment to one dict entry
            self.intensites_per_experiment_lfq = {
                exp: self.all_intensities_lfq[
                    [f"LFQ intensity {rep}" for rep in self.replicates[exp]]
                ] for exp in self.replicates
            }

        else:
            self.all_intensities_lfq = self.all_intensities_raw
            self.intensites_per_experiment_lfq = self.intensites_per_experiment_raw

        # create a dict matching the raw and lfq dfs by string
        self.all_intensities_dict = {"lfq": self.all_intensities_lfq, "raw": self.all_intensities_raw}
        self.intensities_per_experiment_dict = {
            "lfq": self.intensites_per_experiment_lfq, "raw": self.intensites_per_experiment_raw
        }

        # create all sets that are required for plotting
        # this could be turned into a property
        self.protein_ids = ddict(dict)
        self.whole_experiment_protein_ids = {}
        for experiment in self.replicates:
            exp_prot_ids = set()
            for rep in self.replicates[experiment]:
                # TODO what to use for venns: raw or lfq?
                idx = self.all_intensities_raw["Intensity " + rep] > 0
                idx = [i for i in idx.index if idx[i]]
                rep_set = set(self.df_protein_names.loc[idx, "Protein name"])
                self.protein_ids[experiment][rep] = rep_set
                exp_prot_ids |= rep_set
            self.whole_experiment_protein_ids[experiment] = exp_prot_ids

        # set all result dirs
        if self.configs.get("venn_results", False):
            # path for venn diagrams
            self.file_dir_venn = os.path.join(self.start_dir, "venn")
            # create file structure and folder
            os.makedirs(self.file_dir_venn, exist_ok=True)
        if self.configs.get("descriptive_plots", False):
            # path for descriptive plots
            self.file_dir_descriptive = os.path.join(self.start_dir, "descriptive")
            # create file structure and folder
            os.makedirs(self.file_dir_descriptive, exist_ok=True)
        if self.configs.get("go_analysis", False):
            self.file_dir_go_analysis = os.path.join(self.start_dir, "go_analysis")
            # create file structure and folder
            os.makedirs(self.file_dir_go_analysis, exist_ok=True)

    @classmethod
    def from_MQInitializer(cls, mqinti_instance, loglevel=logging.DEBUG):
        return cls(
            start_dir = mqinti_instance.start_dir,
            replicates = mqinti_instance.replicates,
            configs = mqinti_instance.configs,
            df_protein_names = mqinti_instance.df_protein_names,
            df_peptide_names = mqinti_instance.df_peptide_names,
            interesting_proteins = mqinti_instance.interesting_proteins,
            interesting_receptors = mqinti_instance.interesting_receptors,
            go_analysis_gene_names = mqinti_instance.go_analysis_gene_names,
            loglevel=loglevel
            )

    @property
    def replicate_representation(self):
        if self._replicate_representation is None:
            d = {}
            for experiment in self.replicates:
                experiment_rep = experiment.rstrip("0123456789_")
                experiment_rep = " ".join(experiment_rep.split("_"))
                d[experiment] = experiment_rep
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
        # create all results according to configs
        if self.configs.get("venn_results", False):
            self.create_venn_results()

        if self.configs.get("descriptive_plots", False):
            self.create_descriptive_plots()

        if self.configs.get("go_analysis", False):
            self.create_go_analysis()

        if self.configs.get("buffer_score", False):
            pass
        # create report summarizing the analysis results
        self.create_report()

    def save_venn(self, ex: str, sets, set_names):
        creates_figure = True
        plt.close("all")
        plt.figure(figsize=(14, 7))
        plt.title(self.replicate_representation[ex], fontsize=VENN_TITLE_FONT_SIZE)

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
                self.file_dir_venn, f"venn_replicate_{self.replicate_representation[ex].replace(' ', '_')}" + FIG_FORMAT
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

    def save_bar_venn(self, ex: str, named_sets):
        plt.close("all")
        # ax1.set_title("size")
        # ax2.set_title("selected samples")

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
        try:
            name = self.replicate_representation[ex]
        except KeyError:
            name = ex
        fig.suptitle(name, fontsize=20)
        # create the bar plot
        ax1.bar(x, heights)
        # add text to the bar plot
        for x_level, height in zip(x, heights):
            ax1.text(x_level, max(heights) / 2, height, verticalalignment='center', horizontalalignment='center')

        # create the line plots
        for x_level, y in zip(x, ys):
            # we just want to draw a straight line every time so we repeat x as often as needed
            ax2.plot([x_level] * len(y), y, linestyle="-", color="gray", marker=".")
        # replace the yticks with the names of the samples
        ax2.set_yticks([i for i in range(len(y_mappings))])
        ax2.set_yticklabels(y_mappings)

        # save the result with the according name
        res_path = os.path.join(
            self.file_dir_venn, f"venn_bar_{name.replace(' ', '_')}" + FIG_FORMAT
        )
        fig.savefig(res_path, dpi=200, bbox_inches="tight")
        plt.close("all")

    def save_venn_names(self, named_sets):
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

    def create_venn_results(self):
        self.logger.info("Creating venn diagrams")

        # create venn diagrams comparing all replicates within an experiment
        for ex in self.protein_ids:
            set_names = self.protein_ids[ex].keys()
            sets = self.protein_ids[ex].values()
            # save the resulting venn diagram
            self.save_venn(ex, sets, set_names)
            # save the sets to txt files
            self.save_venn_names(self.protein_ids[ex])
            # create a mixture of bar and venn diagram
            self.save_bar_venn(ex, self.protein_ids[ex])
        # create venn diagrams comparing all pellets with supernatants  # TODO
        # create venn diagrams comparing all experiments
        # first compare only proteins which are found between all replicates
        experiment_intersection_sets = {
            exp: set.intersection(*(self.protein_ids[exp][rep] for rep in self.protein_ids[exp]))
            for exp in self.protein_ids
        }
        self.save_bar_venn("All experiments intersection", experiment_intersection_sets)
        # then compare all proteins that are found at all
        experiment_union_sets = {
            exp: set.union(*(self.protein_ids[exp][rep] for rep in self.protein_ids[exp]))
            for exp in self.protein_ids
        }
        self.save_bar_venn("All experiments union", experiment_union_sets)

    def plot_detection_counts(self, df_to_use: str = "raw"):
        plt.close("all")
        # determine number of rows and columns in the plot based on the number of experiments
        n_rows_experiment, n_cols_experiment = get_number_rows_cols_for_fig(self.replicates)
        fig, axarr = plt.subplots(n_rows_experiment, n_cols_experiment, sharex=True, squeeze=True,
                                  figsize=(4 * n_rows_experiment, 4 * n_cols_experiment))
        fig.suptitle(f"Detection counts from {df_to_use} intensities")
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

    def plot_number_of_detected_proteins(self, df_to_use: str = "raw"):
        plt.close("all")
        # determine number of rows and columns in the plot based on the number of experiments
        n_rows_experiment, n_cols_experiment = get_number_rows_cols_for_fig(self.replicates)
        fig, axarr = plt.subplots(n_rows_experiment, n_cols_experiment,
                                  figsize=(5 * n_cols_experiment, 3 * n_rows_experiment))
        fig.suptitle(f"Number of detected proteins from {df_to_use} intensity")

        for experiment, ax in zip(self.replicates, axarr.flat):
            plt.close("all")
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
                labels.append(col.replace("Intensity ", ""))
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

    def plot_intensity_histograms(self, df_to_use: str = "raw"):
        plt.close("all")
        n_rows_replicates, n_cols_replicates = get_number_rows_cols_for_fig([x for l in self.replicates.values() for x in l])
        # make a intensity histogram for every replicate
        fig, axarr = plt.subplots(n_rows_replicates, n_cols_replicates,
                                  figsize=(5 * n_cols_replicates, 5 * n_rows_replicates))
        fig_index = 0
        fig.suptitle(f"Replicate {df_to_use} Intensity histograms")
        for experiment in self.replicates:
            intensities = self.intensities_per_experiment_dict[df_to_use][experiment]
            for col in intensities:
                mask = intensities[col] > 0
                ax = axarr.flat[fig_index]
                ax.set_title(f"{col}".replace(self.configs['descriptive_intensity_col'], ""))
                ax.hist(intensities[col][mask],
                        bins=np.logspace(
                            np.log10(intensities[col][mask].min()),
                            np.log10(intensities[col][mask].max()),
                            25
                           )
                        )

                ax.set_xscale('log')
                fig_index += 1

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        res_path = os.path.join(self.file_dir_descriptive, "intensity_histograms" + FIG_FORMAT)
        fig.savefig(res_path, dpi=200, bbox_inches="tight")

    def plot_scatter_replicates(self):
        plt.close("all")
        for experiment in self.replicates:
            # scatter plot of all replicates vs all replicates
            fig, ax = plt.subplots(1, 1, figsize=(7, 7))
            for rep1, rep2 in combinations(self.replicates[experiment], 2):
                x1 = self.df_protein_names[self.configs['descriptive_intensity_col'] + rep1]
                x2 = self.df_protein_names[self.configs['descriptive_intensity_col'] + rep2]
                mask = np.logical_or(x1 != 0, x2 != 0)
                exp = r"$r^{2}$"
                ax.scatter(x1[mask] + 1e2, x2[mask] + 1e2, label=f"{rep1} vs {rep2}, "
                                                      fr"{exp}: {stats.pearsonr(x1[mask], x2[mask])[0] ** 2}",
                            alpha=0.5, marker=".")
                ax.set_xscale('log')
                ax.set_yscale('log')

            fig.legend()

            res_path = os.path.join(self.file_dir_descriptive,
                                    f"{self.replicate_representation[experiment].replace(' ', '_')}_scatter" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")

    def plot_rank(self, df_to_use: str = "raw"):
        plt.close("all")
        all_pathway_proteins = set.union(*(set(x) for x in self.interesting_proteins.values()))
        for experiment in self.replicates:
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
            if rank_identified_proteins:
                median_pathway_rank = int(np.median(rank_identified_proteins))
                median_intensity = m_intensity.iloc[median_pathway_rank]

            fig, ax = plt.subplots(1, 1, figsize=(14, 7))
            # plot the non pathway proteins
            x = [dic[protein][0] for protein in non_pathway_proteins]
            y = [dic[protein][1] for protein in non_pathway_proteins]
            ax.scatter(x, y, c=f"darkgray", s=10, alpha=0.3, marker=".", label="no pathway")
            # plot all proteins of a specific pathway
            for i, (pathway, proteins) in enumerate(self.interesting_proteins.items()):
                proteins = set(proteins) & found_proteins
                x = [dic[protein][0] for protein in proteins]
                y = [dic[protein][1] for protein in proteins]
                ax.scatter(x, y, c=f"C{i}", s=80, alpha=0.6, marker=".", label=pathway)
            #
            ax.set_yscale("log")
            if rank_identified_proteins:
                xmin, xmax = ax.get_xbound()
                xm = (median_pathway_rank + abs(xmin)) / (abs(xmax) + abs(xmin))
                ymin, ymax = ax.get_ybound()
                ym = (np.log10(median_intensity) - np.log10(ymin) ) / (np.log10(ymax) - np.log10(ymin))
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
            ax.set_ylabel(f"{df_to_use} intensity mean")
            fig.legend()

            res_path = os.path.join(self.file_dir_descriptive,
                                     f"{self.replicate_representation[experiment].replace(' ', '_')}_rank" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")

    def plot_relative_std(self, df_to_use: str = "raw"):
        plt.close("all")
        for experiment in self.replicates:
            intensities = self.intensities_per_experiment_dict[df_to_use][experiment]
            y = intensities.copy()
            # below cutoff remove zeros, above threshold include 0
            cutoff = 1e6
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
            fig, ax = plt.subplots(1, 1, figsize=(14, 7))
            ax.scatter(y.mean(axis=1)[mask], relative_std_percent[mask], c=colors[mask], marker=".", s=200, alpha=0.8)
            ax.set_xscale("log")
            ax.set_xlabel("Mean raw intensity")
            ax.set_ylabel("Relative Stadartdeviation [%]")
            ax.axvline(cutoff, color="black", alpha=0.5)
            xmin, xmax = ax.get_xbound()
            cumulative_count = 0
            for i, bin_ in enumerate(bins):
                cumulative_count += color_counts[cm[i]]
                ax.axhline(bin_, color=cm[i])
                ax.text(xmin, bin_, cumulative_count)

            res_path = os.path.join(self.file_dir_descriptive,
                                    f"{self.replicate_representation[experiment].replace(' ', '_')}_rel_std" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")

    def plot_pathway_analysis(self, df_to_use: str = "raw"):
        ex_list = list(self.replicates.keys())
        for pathway in self.interesting_proteins:
            plt.close("all")
            found_proteins = set(self.interesting_proteins[pathway])
            found_proteins &= set(self.all_intensities_dict[df_to_use].index)
            found_proteins = sorted(list(found_proteins))
            n_rows, n_cols = get_number_rows_cols_for_fig(found_proteins)
            fig, axarr = plt.subplots(n_rows, n_cols, figsize=(4 * n_rows, 4 * n_cols))
            fig.suptitle(pathway)
            for protein, ax in zip(found_proteins, axarr.flat):
                ax.set_title(protein)
                heights = []
                for idx, experiment in enumerate(self.replicates):
                    protein_intensities = self.intensites_per_experiment_raw[experiment].loc[protein, :]
                    ax.scatter([idx] * len(protein_intensities), protein_intensities, label=f"{experiment}")
                    # self.logger.debug(f"{protein}: min: {min_i}, max: {max_i}")
                    heights.append(np.max(protein_intensities))
                # ax.set_ylim((0.9 * min_i, 1.1 * max_i))
                # ax.set_yscale("log")
                # ax.legend()
                n_annotations = 0
                for e1, e2 in combinations(self.replicates, 2):
                    v1 = self.intensites_per_experiment_raw[e1].loc[protein, :].astype(float)
                    v2 = self.intensites_per_experiment_raw[e2].loc[protein, :].astype(float)
                    log_transform = False
                    plot_ns = False
                    threshhold = 0.05
                    if log_transform:
                        v1 = np.log2(v1 + 10)
                        v2 = np.log2(v2 + 10)
                    # only perform the test if more than 3 data points are available for both experiments
                    if len(v1[v1 > 0]) < 3 or len(v2[v2 > 0]) < 3:
                        # self.logger.debug(f"skipping test for: {protein}, {e1} vs {e2}")
                        continue
                    test = stats.ttest_ind(v1, v2, equal_var=self.configs.get("equal_var", False))
                    if plot_ns or test[1] < threshhold:
                        barplot_annotate_brackets(ax, ex_list.index(e1), ex_list.index(e2), test[1], range(len(ex_list)),
                                                  heights, dh=0.05 + 0.1 * n_annotations)
                        n_annotations += 1

            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            res_path = os.path.join(self.file_dir_descriptive, f"pathway_analysis_{pathway}" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")

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
            # plot intersection
            ax6.scatter(intersection["ex1"], intersection["ex2"], s=8, alpha=0.6, marker=".")
            ax6.set_xscale('log')
            ax6.set_yscale('log')
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
            ax9.set_xscale("log")
            ax9.set_yscale("log")
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
            # sys.exit(0)

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

    def create_descriptive_plots(self):
        self.logger.info("Creating descriptive plots")
        self.plot_detection_counts()
        self.plot_number_of_detected_proteins()
        self.plot_intensity_histograms()
        self.plot_relative_std("lfq")
        self.plot_rank()
        self.plot_pathway_analysis("lfq")
        self.plot_pathway_proportions()
        self.plot_scatter_replicates()
        self.plot_experiment_comparison("lfq")

    def calculate_buffer_score(self):
        pass

    def create_report(self):
        pass

    def create_go_analysis(self, df_to_use: str = "raw"):
        self.logger.info("Creating go analysis plots")
        plt.close("all")

        # all proteins that were detected in any replicate
        background = set(self.all_intensities_dict[df_to_use].index)
        heights = ddict(list)
        test_results = ddict(list)

        bar_width = 0.25
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
                chi2, p, dof, ex = stats.chi2_contingency(table, correction=True)
                # self.logger.debug(f"{chi2}, {dof}, {ex}")
                # self.logger.debug(f"{compartiment} {experiment}: table: {table} fisher: {pvalue:.4f}, chi2: {p:.4f}")
                test_results[experiment].append(pvalue)

        fig, ax = plt.subplots(1, 1, figsize=(7, 7))

        for i, experiment in enumerate(heights):
            y_pos = np.array([(i + (len(heights) + 1) * x) * bar_width for x in range(len(self.go_analysis_gene_names))])
            ax.barh(y_pos, heights[experiment], height=bar_width, edgecolor='white', label=experiment)
            for x, y, text in zip(heights[experiment], y_pos, test_results[experiment]):
                ax.text(x, y, f" p: {text:.4f}", verticalalignment="center", fontsize="8")

        ax.set_ylabel('compartiment')
        ax.set_xlabel('number of proteins')
        # add yticks on the middle of the group bars
        ax.set_yticks(np.array([((len(heights) * 1.5 - 1.5) * x + 2) * bar_width for x in range(len(self.go_analysis_gene_names))]))
        # replace the y ticks with the compartiment names
        ax.set_yticklabels([x for x in self.go_analysis_gene_names])
        plt.legend()
        res_path = os.path.join(self.file_dir_go_analysis, "go_analysis" + FIG_FORMAT)
        fig.savefig(res_path, dpi=200, bbox_inches="tight")
