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
from Utils import barplot_annotate_brackets, get_number_rows_cols_for_fig


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
        interesting_proteins, interesting_receptors, go_analysis_gene_names
        ):
        super().__init__(self.__class__.__name__)
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

        # extract all lfq intensites from the dataframe
        self.all_intensities_lfq = self.df_protein_names[
            [f"LFQ intensity {rep}" for exp in self.replicates for rep in self.replicates[exp]]
        ]
        # filter all rows where all intensities are 0
        mask = (self.all_intensities_lfq != 0).sum(axis=1) != 0
        self.all_intensities_lfq = self.all_intensities_lfq[mask]

        # write all intensites of every experiment to one dict entry
        self.intensites_per_experiment_lfg = {
            exp: self.all_intensities_lfq[
                [f"LFQ intensity {rep}" for rep in self.replicates[exp]]
            ] for exp in self.replicates
        }

    @classmethod
    def from_MQInitializer(cls, mqinti_instance):
        return cls(
            start_dir = mqinti_instance.start_dir,
            replicates = mqinti_instance.replicates,
            configs = mqinti_instance.configs,
            df_protein_names = mqinti_instance.df_protein_names,
            df_peptide_names = mqinti_instance.df_peptide_names,
            interesting_proteins = mqinti_instance.interesting_proteins,
            interesting_receptors = mqinti_instance.interesting_receptors,
            go_analysis_gene_names = mqinti_instance.go_analysis_gene_names
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

    def create_venn_results(self):
        self.logger.info("Creating venn diagrams")
        file_dir_venn = os.path.join(self.start_dir, "venn")
        # create file structure and folder
        os.makedirs(file_dir_venn, exist_ok=True)

        def save_venn(ex: str, sets, set_names):
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
                    file_dir_venn, f"venn_replicate_{self.replicate_representation[ex].replace(' ', '_')}" + FIG_FORMAT
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

        def save_bar_venn(ex: str, named_sets):
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
                file_dir_venn, f"venn_bar_{name.replace(' ', '_')}" + FIG_FORMAT
            )
            fig.savefig(res_path, dpi=200, bbox_inches="tight")
            plt.close("all")

        def venn_names(named_sets):
            names = set(named_sets)
            for i in range(1, len(named_sets) + 1):
                for to_intersect in combinations(sorted(named_sets), i):
                    others = names.difference(to_intersect)
                    intersected = set.intersection(*(named_sets[k] for k in to_intersect))
                    unioned = set.union(*(named_sets[k] for k in others)) if others else set()
                    yield to_intersect, others, intersected - unioned

        def save_venn_names(named_sets):
            for intersected, unioned, result in venn_names(named_sets):
                # create name based on the intersections and unions that were done
                intersected_name = "&".join(sorted(intersected))
                unioned_name = "-"+"-".join(sorted(unioned)) if unioned else ""
                res_path = os.path.join(
                    file_dir_venn,
                    f"set_{intersected_name}{unioned_name}.txt"
                )
                # write all names line by line into the file
                with open(res_path, "w") as out:
                    for re in result:
                        out.write(re + "\n")

        # create all sets that are required for plotting
        protein_ids = ddict(dict)
        whole_experiment_protein_ids = {}
        for experiment in self.replicates:
            exp_prot_ids = set()
            for rep_name, rep in zip(self.replicates[experiment], self.all_intensities_raw):
                # TODO what to use for venns
                idx = self.all_intensities_raw[rep] > 0
                idx = [i for i in idx.index if idx[i]]
                rep_set = set(self.df_protein_names.loc[idx, "Protein name"])
                protein_ids[experiment][rep_name] = rep_set
                exp_prot_ids |= rep_set
            whole_experiment_protein_ids[experiment] = exp_prot_ids

        # create venn diagrams
        # create venn diagrams comparing all replicates within an experiment
        for ex in protein_ids:
            set_names = protein_ids[ex].keys()
            sets = protein_ids[ex].values()
            # save the resulting venn diagram
            save_venn(ex, sets, set_names)
            # TODO should these results be different options
            # save the sets to txt files
            save_venn_names(protein_ids[ex])
            # create a mixture of bar and venn diagram
            save_bar_venn(ex, protein_ids[ex])
        # create venn diagrams comparing all pellets with supernatants
        # create venn diagrams comparing all experiments
        # first compare only proteins which are found between all replicates
        experiment_sets = {exp: set.intersection(*(protein_ids[exp][rep] for rep in protein_ids[exp])) for exp in protein_ids}
        save_bar_venn("All experiments intersection", experiment_sets)
        # then compare all proteins that are found at all
        experiment_sets = {exp: set.union(*(protein_ids[exp][rep] for rep in protein_ids[exp])) for exp in protein_ids}
        save_bar_venn("All experiments union", experiment_sets)

    def create_descriptive_plots(self):
        self.logger.info("Creating descriptive plots")
        file_dir_descriptive = os.path.join(self.start_dir, "descriptive")
        # create file structure and folder
        os.makedirs(file_dir_descriptive, exist_ok=True)
        plt.close("all")
        # determine number of rows and columns in the plot based on the number of experiments
        n_rows_experiment, n_cols_experiment = get_number_rows_cols_for_fig(self.replicates)

        # plot showing the amount of proteins that were detected x amount of times
        n_rows_replicates, n_cols_replicates = get_number_rows_cols_for_fig([x for l in self.replicates.values() for x in l])

        fig, axarr = plt.subplots(n_rows_experiment, n_cols_experiment, sharex=True, squeeze=True, figsize=(5 * n_rows_experiment, 3 * n_cols_experiment))
        fig.suptitle("Detection counts")
        fig1, axarr1 = plt.subplots(n_rows_experiment, n_cols_experiment, figsize=(14, 7))
        fig1.suptitle("Number of detected proteins")

        # make a intensity histogram for every replicate
        # TODO scale this fig size with number of rows
        fig7, axarr7 = plt.subplots(n_rows_replicates, n_cols_replicates, figsize=(20, 15))
        fig7_index = 0
        fig7.suptitle("Replicate Intensity histograms")

        # make a intensity histogram for every replicate
        fig8, ax8 = plt.subplots(1, 1, figsize=(14, 7))
        fig8.suptitle("Number of identified proteins per pathway")
        pathway_counts = ddict(lambda: ddict(int))
        pathway_colors = {
            'Signaling of Hepatocyte Growth Factor Receptor': "navy",
            'TGF-beta signaling pathway': "royalblue",
            'EGF Signaling Pathway': "skyblue",
            'EPO Signaling Pathway': "darkslateblue",
            'GAS6 Signaling Pathway': "mediumpurple"
        }

        ex_list = list(self.replicates.keys())
        for experiment, ax, ax1 in zip(self.replicates, axarr.flat, axarr1.flat):
            plt.close("all")
            intensities = self.intensites_per_experiment_raw[experiment]

            # plot 1
            # from 0 to number of replicates, how often was each protein detected
            counts = (intensities > 0).sum(axis=1)
            counts = counts[counts > 0]
            ax.set_title(f"{experiment}, total detected: {len(counts)}")
            bins = [x - 0.5 for x in range(1, self.max_number_replicates + 2)]
            ax.hist(counts, bins=bins)
            max_val = max([len(counts[counts == value]) for value in counts.unique()])
            for x, value in zip(bins, sorted(counts.unique())):
                ax.text(x + 0.5, max_val / 2, len(counts[counts == value]),
                        verticalalignment='center', horizontalalignment='center')
            ax.set_xticks(range(1, self.max_number_replicates + 1))
            ax.set_ylabel("Counts")

            # plot 2
            # how many proteins were detected per replicate and in total
            heights = [len(counts), 0]
            labels = ["Total", " "]
            for col in intensities:
                h = len(intensities[col][intensities[col] > 0])
                heights.append(h)
                labels.append(col)
            mean_height = np.mean(heights[2:])
            # self.logger.debug(labels)
            ax1.bar([x for x in range(len(heights))], heights)
            ax1.text(1, mean_height * 1.01, f"{mean_height:.2f}", horizontalalignment='center')
            ax1.set_title(f"{experiment}")
            ax1.axhline(mean_height, linestyle="--", color="black", alpha=0.6)
            # ax1.set_xticklabels(labels, rotation=45)
            ax1.set_ylabel("Counts")
            #ax1.tick_params(rotation=70)

            # plot 3
            # scatter plot of all replicates vs all replicates
            fig3, ax3 = plt.subplots(1, 1, figsize=(7, 7))
            for rep1, rep2 in combinations(self.replicates[experiment], 2):
                x1 = self.df_protein_names[self.configs['descriptive_intensity_col'] + rep1]
                x2 = self.df_protein_names[self.configs['descriptive_intensity_col'] + rep2]
                mask = np.logical_or(x1 != 0, x2 != 0)
                ax3.scatter(x1[mask] + 1e2, x2[mask] + 1e2, label=f"{rep1} vs {rep2}, "
                                                      f"r: {stats.pearsonr(x1[mask], x2[mask])[0] ** 2}",
                            alpha=0.5, marker=".")
                ax3.set_xscale('log')
                ax3.set_yscale('log')

            fig3.legend()

            res3_path = os.path.join(file_dir_descriptive,
                                     f"{self.replicate_representation[experiment].replace(' ', '_')}_scatter" + FIG_FORMAT)
            # TODO add r2
            fig3.savefig(res3_path, dpi=200, bbox_inches="tight")

            # plot 4
            # protein ranks vs intensity
            fig4, ax4 = plt.subplots(1, 1, figsize=(14, 7))
            m_intensity = intensities.mean(axis=1).sort_values(ascending=False)
            m_intensity = m_intensity[m_intensity > 0]

            protein_colors = {protein: pathway_colors.get(pathway, "black") for pathway, proteins in self.interesting_proteins.items() for protein in proteins}

            colors = [protein_colors.get(index, "darkgray") for index in m_intensity.index]
            sizes = [10 if col == "black" or col == "darkgray" else 100 for col in colors]
            colors_pathway = {v: k for k, v in pathway_colors.items()}
            # TODO this is done at mean intensity level now
            for col in colors:
                pathway_counts[colors_pathway.get(col, "Remaining")][experiment] += 1

            ax4.scatter(range(1, len(m_intensity) + 1), m_intensity, c=colors, s=sizes, alpha=0.3, marker=".")
            ax4.set_yscale('log')
            ax4.set_xlabel("Protein rank")
            ax4.set_ylabel("Raw intensity mean")

            res4_path = os.path.join(file_dir_descriptive,
                                     f"{self.replicate_representation[experiment].replace(' ', '_')}_rank" + FIG_FORMAT)
            fig4.savefig(res4_path, dpi=200, bbox_inches="tight")

            # plot 5
            # intensity vs relative standard deviation
            fig5, ax5 = plt.subplots(1, 1, figsize=(14, 7))
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
            mask = ~relative_std_percent.isna()
            ax5.scatter(y.mean(axis=1)[mask], relative_std_percent[mask], c=colors[mask], marker=".", s=200, alpha=0.8)
            ax5.set_xscale('log')
            ax5.set_xlabel("Mean raw intensity")
            ax5.set_ylabel("Relative Stadartdeviation [%]")
            ax5.axvline(cutoff, color="black", alpha=0.5)
            ax5.axhline(10, color=cm[0])
            ax5.axhline(20, color=cm[1])
            ax5.axhline(30, color=cm[2])

            res5_path = os.path.join(file_dir_descriptive,
                                     f"{self.replicate_representation[experiment].replace(' ', '_')}_rel_std" + FIG_FORMAT)
            fig5.savefig(res5_path, dpi=200, bbox_inches="tight")

            for col in intensities:
                ax7 = axarr7.flat[fig7_index]
                mask = intensities[col] > 0
                ax7.set_title(f"{col}".replace(self.configs['descriptive_intensity_col'], ""))
                # ax7.hist(np.log2(intensities[col][mask]))
                ax7.hist(intensities[col][mask],
                         bins=np.logspace(
                             np.log10(intensities[col][mask].min()),
                             np.log10(intensities[col][mask].max()),
                             25
                            )
                         )

                ax7.set_xscale('log')
                fig7_index += 1

        res_path = os.path.join(file_dir_descriptive, f"detected_counts" + FIG_FORMAT)
        fig.savefig(res_path, dpi=200, bbox_inches="tight")

        res1_path = os.path.join(file_dir_descriptive, "detection_per_replicate" + FIG_FORMAT)
        fig1.savefig(res1_path, dpi=200, bbox_inches="tight")

        res7_path = os.path.join(file_dir_descriptive, "intensity_histograms" + FIG_FORMAT)
        fig7.savefig(res7_path, dpi=200, bbox_inches="tight")

        experiment_proportion = ddict(lambda: ddict(int))
        for pathway in pathway_counts:
            counts = np.array(list(pathway_counts[pathway].values()))
            for experiment, count in zip(pathway_counts[pathway], counts):
                experiment_proportion[experiment][pathway] = count

        df_total_counts = pd.DataFrame(experiment_proportion).T
        df_proportion = df_total_counts / df_total_counts.sum()
        bottom_cumsum = df_proportion.cumsum()
        for i, index in enumerate(df_proportion.index):
            if i == 0:
                bottom = [0] * len(df_proportion.loc[index])
            else:
                bottom = bottom_cumsum.iloc[i - 1]
            ax8.bar(range(len(df_proportion.loc[index])), df_proportion.loc[index],
                    bottom=bottom, label=df_proportion.index[i], tick_label=df_proportion.columns)
            for j, text_index in enumerate(df_total_counts.columns):
                height = (bottom_cumsum.iloc[i, j] + bottom[j]) / 2
                ax8.text(j, height, df_total_counts.iloc[i, j], horizontalalignment="center")
        # fig8.xticks(rotation=45)
        fig8.legend()
        res8_path = os.path.join(file_dir_descriptive, "pathway_proportions" + FIG_FORMAT)
        fig8.savefig(res8_path, dpi=200, bbox_inches="tight")

        for pathway in pathway_colors:
            plt.close("all")
            found_proteins = set(self.interesting_proteins[pathway])
            found_proteins &= set(self.all_intensities_raw.index)
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
                    test = stats.ttest_ind(v1, v2)
                    if plot_ns or test[1] < threshhold:
                        barplot_annotate_brackets(ax, ex_list.index(e1), ex_list.index(e2), test[1], range(len(ex_list)),
                                                  heights, dh=0.05 + 0.1 * n_annotations)
                        n_annotations += 1

            res_path = os.path.join(file_dir_descriptive, f"pathway_analysis_{pathway}" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")

        for ex1, ex2 in combinations(self.replicates, 2):
            protein_intensities_ex1 = self.intensites_per_experiment_raw[ex1]
            counts_ex1 = (protein_intensities_ex1 > 0).sum(axis=1) == len(protein_intensities_ex1.columns)
            protein_intensities_ex2 = self.intensites_per_experiment_raw[ex1]
            counts_ex2 = (protein_intensities_ex2 > 0).sum(axis=1) == len(protein_intensities_ex2.columns)
            # combine intersections
            intersection = pd.concat(
                [protein_intensities_ex1.mean(axis=1).rename("ex1"), protein_intensities_ex2.mean(axis=1).rename("ex2")]
                , axis=1)
            intersection = intersection[counts_ex1 & counts_ex2]
            plt.close("all")
            fig6, ax6 = plt.subplots(1, 1, figsize=(7, 7))
            # plot intersection
            ax6.scatter(intersection["ex1"], intersection["ex2"], s=8, alpha=0.6, marker=".")
            ax6.set_xscale('log')
            ax6.set_yscale('log')
            ax6.set_xlabel(ex1)
            ax6.set_ylabel(ex2)

            # TODO add r2
            res_path = os.path.join(file_dir_descriptive,
                                    f"{self.replicate_representation[ex1].replace(' ', '_')}_vs_{self.replicate_representation[ex2].replace(' ', '_')}_intersection" + FIG_FORMAT)
            fig6.savefig(res_path, dpi=200, bbox_inches="tight")
            plt.close("all")
            # combine total
            protein_intensities_ex1 = protein_intensities_ex1.mean(axis=1).rename("ex1")
            protein_intensities_ex2 = protein_intensities_ex2.mean(axis=1).rename("ex2")
            conc = pd.concat([protein_intensities_ex1, protein_intensities_ex2], axis=1)
            conc = conc[conc > 0]
            # plot total
            plt.close("all")
            fig9, ax9 = plt.subplots(1, 1, figsize=(7, 7))
            ax9.scatter(conc["ex1"], conc["ex2"], s=8, alpha=0.6, marker=".")
            ax9.set_xscale('log')
            ax9.set_yscale('log')
            ax9.set_xlabel(ex1)
            ax9.set_ylabel(ex2)
            # TODO add r2
            res_path = os.path.join(file_dir_descriptive,
                                    f"{self.replicate_representation[ex1].replace(' ', '_')}_vs_{self.replicate_representation[ex2].replace(' ', '_')}_total" + FIG_FORMAT)
            plt.savefig(res_path, dpi=200, bbox_inches="tight")
            plt.close("all")
            # sys.exit(0)

    def calculate_buffer_score(self):
        pass

    def create_report(self):
        pass

    def create_go_analysis(self):
        self.logger.info("Creating go analysis plots")
        file_dir_go_analysis = os.path.join(self.start_dir, "go_analysis")
        # create file structure and folder
        os.makedirs(file_dir_go_analysis, exist_ok=True)
        plt.close("all")

        # all proteins that were detected in any replicate
        background = set(self.all_intensities_raw.index)
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
                # create df with intensity means for specific experiment ofer all replicates
                mean_intensity = self.intensites_per_experiment_raw[experiment].mean(axis=1)
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
        res_path = os.path.join(file_dir_go_analysis, "go_analysis" + FIG_FORMAT)
        fig.savefig(res_path, dpi=200, bbox_inches="tight")
