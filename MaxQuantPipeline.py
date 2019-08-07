import pandas as pd
import numpy as np
import os
import sys
from Logger import Logger
import tkinter as tk
from tkinter import filedialog
import argparse
from matplotlib_venn import venn3, venn2
import matplotlib.pyplot as plt
from scipy import stats
from itertools import combinations
from collections import defaultdict as ddict
from ruamel.yaml import YAML
from collections.abc import Sized


# plt.style.use('ggplot')

FIG_FORMAT = ".pdf"
# Venn diagram settings
# TODO figsize
VENN_TITLE_FONT_SIZE = 20
VENN_SET_LABEL_FONT_SIZE = 16
VENN_SUBSET_LABEL_FONT_SIZE = 14
# descriptive plots settings


def get_number_rows_cols_for_fig(obj):
    if isinstance(obj, Sized):
        obj = len(obj)
    n_rows, n_cols = 0, 0
    while n_rows * n_cols < obj:
        if n_cols <= n_rows:
            n_cols += 1
        else:
            n_rows += 1
    return n_rows, n_cols


def barplot_annotate_brackets(ax, num1, num2, data, center, height, yerr=None, dh=.05, barh=.05, fs=None,
                              maxasterix=None):
    """
    From: https://stackoverflow.com/questions/11517986/indicating-the-statistically-significant-difference-in-bar-graph
    Annotate barplot with p-values.

    :param ax: axis of plot to put the annotaion brackets
    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """

    if type(data) is str:
        text = data
    else:
        # * is p < 0.05
        # ** is p < 0.005
        # *** is p < 0.0005
        # etc.
        text = ''
        p = .05

        while data < p:
            text += '*'
            p /= 10.

            if maxasterix and len(text) == maxasterix:
                break

        if len(text) == 0:
            text = 'n. s.'

    lx, ly = center[num1], max(height)
    rx, ry = center[num2], max(height)

    if yerr:
        ly += yerr[num1]
        ry += yerr[num2]

    ax_y0, ax_y1 = ax.get_ylim()
    log = False
    if log:
        dh *= 10 ** (ax_y1 - ax_y0)
        barh *= 10 ** (ax_y1 - ax_y0)
    else:
        dh *= (ax_y1 - ax_y0)
        barh *= (ax_y1 - ax_y0)

    y = max(ly, ry) + dh

    barx = [lx, lx, rx, rx]
    bary = [y, y + barh, y + barh, y]
    mid = ((lx + rx) / 2, y + barh)

    ax.plot(barx, bary, c='black')

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    ax.text(*mid, text, **kwargs)


class MaxQuantPipeline(Logger):
    def __init__(self, dir_: str, file_path_yml: str = None):
        super().__init__(self.__class__.__name__)
        self._replicates = None
        self._replicate_representation = None
        self._min_number_replicates = None
        self._max_number_replicates = None
        self._replicates_representation = None

        self.script_loc = os.path.dirname(os.path.realpath(__file__))
        self.path_pipeline_config = os.path.join(self.script_loc, "config")

        self.yml_file_name_tmp = "config_tmp.yml"
        self.yml_file_name = "config.yml"
        self.default_yml_name = "ms_analysis_default.yml"

        # make sure to be on the right level and set starting dir
        if os.path.split(os.path.split(dir_)[0])[-1] == "txt":
            self.logger.debug("Removing txt ending from path")
            self.start_dir = os.path.join(os.path.split(os.path.split(dir_)[0])[0])
        else:
            self.start_dir = dir_
        self.logger.info(f"Starting dir: {self.start_dir}")

        # if no yml file is passed try to guess it or ask for one
        if file_path_yml is None:
            file_path_yml = self.init_yml_path()
        elif file_path_yml.lower() == "default":
            file_path_yml = self.get_default_yml_path()

        # load the config from the yml file
        self.logger.debug(f"yml file location: {file_path_yml}")
        self.yaml = YAML()
        self.yaml.default_flow_style = False
        self.logger.info("loading yml file")
        with open(file_path_yml) as f:
            self.configs = self.yaml.load(f)
        self.logger.debug(f"Config file contents: {self.configs}")

        self.path_config = os.path.join(self.start_dir, "config")
        os.makedirs(self.path_config, exist_ok=True)
        self.update_config_file()

        # read the proteinGroups.txt and peptides.txt
        self.logger.info("Reading proteinGroups.txt")  # TODO also reading another file
        self.df_protein_names, self.df_peptide_names = self.init_dfs_from_txts()
        self.update_config_file()

        # read all proteins and receptors of interest from the config dir
        self.logger.info("Reading proteins and receptors of interest")
        self.interesting_proteins, self.interesting_receptors, self.go_analysis_gene_names = self.init_interest_from_xlsx()

    @property
    def replicates(self):
        if self._replicates is None:
            self._replicates = self.init_replicates(self.df_protein_names.columns)
        return self._replicates

    def init_replicates(self, df_colnames):
        all_reps = sorted([x.replace('Intensity ', '') for x in df_colnames
                           if x.startswith('Intensity ')], key=len, reverse=True)
        _replicates = ddict(list)
        for experiment in self.configs.get("experiments", False):
            if not experiment:
                raise ValueError("Missing experiments key in config file")
            for rep in all_reps:
                if rep.startswith(experiment):
                    _replicates[experiment].append(rep)
        return _replicates

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

    def get_default_yml_path(self):
        self.logger.debug(f"Loading default yml file from: {self.script_loc}, since no file was selected")
        if self.default_yml_name in os.listdir(self.path_pipeline_config):
            yaml_file = os.path.join(self.path_pipeline_config, self.default_yml_name)
        else:
            raise ValueError("Could not find default yaml file. Please select one.")
        return yaml_file

    def init_yml_path(self) -> str:
        def yml_file_dialog() -> str:
            file_path = filedialog.askopenfilename()
            self.logger.debug(f"selected file path: {file_path}")
            if not file_path:
                yaml_file = self.get_default_yml_path()
            elif not file_path.endswith(".yaml") and not file_path.endswith(".yml"):
                raise ValueError("Please select a yaml / yml file")
            else:
                yaml_file = file_path
            return yaml_file

        if "config" in os.listdir(self.start_dir):
            self.logger.debug("Found config dir")
            config_dir = os.path.join(self.start_dir, "config")
            if self.yml_file_name in os.listdir(config_dir):
                self.logger.debug("Found config.yml file in config dir")
                yaml_file = os.path.join(config_dir, self.yml_file_name)
            else:
                yaml_file = yml_file_dialog()
        else:
            yaml_file = yml_file_dialog()
        return yaml_file

    def init_dfs_from_txts(self):
        file_dir_txt = os.path.join(self.start_dir, "txt")
        if not os.path.isdir(file_dir_txt):
            raise ValueError("Directory does not contain a txt dir")
        file_dir_protein_names = os.path.join(file_dir_txt, "proteinGroups.txt")
        file_dir_peptides_names = os.path.join(file_dir_txt, "peptides.txt")
        # make sure protein groups file exists
        if not os.path.isfile(file_dir_protein_names):
            raise ValueError("txt directory does not contain a proteinGroups.txt file")
        if not os.path.isfile(file_dir_peptides_names):
            raise ValueError("txt directory does not contain a peptides.txt file")
        # read protein groups file
        df_protein_names = pd.read_csv(file_dir_protein_names, sep="\t")
        df_peptide_names = pd.read_csv(file_dir_peptides_names, sep="\t")

        # try to automatically determine experimental setup
        if not self.configs.get("experiments", False):
            self.logger.info("No replicates specified in settings file. Trying to infer.")
            # TODO will there every be more than 9 replicates?
            import difflib

            def get_overlap(s1, s2):
                s = difflib.SequenceMatcher(None, s1, s2)
                pos_a, pos_b, size = s.find_longest_match(0, len(s1), 0, len(s2))
                return s1[pos_a:pos_a + size]

            # TODO can the Intensity column always be expected in the file?
            # TODO will the column names always be the same between Intensity and LFQ intensity?
            all_reps = sorted([x.replace('Intensity ', '') for x in df_protein_names.columns
                               if x.startswith('Intensity ')], key=len, reverse=True)
            # make sure the two files contain the same replicate names
            all_reps_peptides = [x.replace('Intensity ', '') for x in df_protein_names.columns
                                 if x.startswith('Intensity ')]
            experiment_analysis_overlap = [x not in all_reps for x in all_reps_peptides]
            if any(experiment_analysis_overlap):
                unmatched = [x for x in all_reps_peptides if experiment_analysis_overlap]
                raise ValueError("Found replicates in peptides.txt, but not in proteinGroups.txt: " +
                                 ", ".join(unmatched))
            #
            overlap = [[get_overlap(re1, re2) if re1 != re2 else "" for re1 in all_reps] for re2 in all_reps]
            overlap_matrix = pd.DataFrame(overlap, columns=all_reps, index=all_reps)
            unclear_matches = []
            replicates = ddict(list)
            for col in overlap_matrix:
                sorted_matches = sorted(overlap_matrix[col].values, key=len, reverse=True)
                best_match = sorted_matches[0]
                replicates[best_match].append(col)
                # check if any other match with the same length could be found
                if any([len(best_match) == len(match) and best_match != match for match in sorted_matches]):
                    unclear_matches.append(best_match)
            for experiment in replicates:
                if len(replicates[experiment]) == 1:
                    rep = replicates.pop(experiment)
                    replicates[rep[0]] = rep
                elif experiment in unclear_matches:
                    self.logger.debug(f"unclear match for experiment: {experiment}")
            self.logger.info(f"determined experiemnts: {replicates.keys()}")
            self.logger.debug(f"number of replicates per experiment:")
            self.logger.debug("\n".join([f"{ex}: {len(replicates[ex])}" for ex in replicates]))
            self.configs["experiments"] = list(replicates.keys())
            self._replicates = replicates

        # TODO properties seem to make no sense here anymore
        if self._replicates is None:
            self._replicates = self.init_replicates(df_protein_names.columns)

        found_replicates = [rep for l in self.replicates.values() for rep in l]
        for df_cols in [df_peptide_names.columns, df_protein_names.columns]:
            intens_cols = [x.replace('Intensity ', '') for x in df_cols if x.startswith('Intensity ')]
            not_found_replicates = [x not in found_replicates for x in intens_cols]
            if any(not_found_replicates):
                unmatched = [x for x in intens_cols if not_found_replicates]
                raise ValueError("Found replicates in peptides.txt, but not in proteinGroups.txt: " +
                                 ", ".join(unmatched))

        # filter all contaminants by removing all rows where any of the 3 columns contains a +
        not_contaminants = (df_protein_names[
                                ["Only identified by site", "Reverse", "Potential contaminant"]] == "+"
                            ).sum(axis=1) == 0
        df_protein_names = df_protein_names[not_contaminants]
        # split the fasta headers
        # first split all fasta headers that contain multiple entries
        sep_ind = df_protein_names["Fasta headers"].str.contains(";")
        # replace all fasta headers with multiple entries with only the first one
        # TODO will there always be a fasta header?
        df_protein_names["Fasta headers"][sep_ind] = df_protein_names["Fasta headers"][sep_ind].str.split(";").apply(pd.Series)[0]
        # split the fasta headers with the pipe symbol
        fasta_col = df_protein_names["Fasta headers"].str.split("|", n=2).apply(pd.Series)
        fasta_col.columns = ["trash", "protein id", "description"]
        # extract the gene name from the description eg: "GN=abcd"
        gene_names_fasta = fasta_col["description"].str.extract(r"(GN=[A-Za-z0-9]*)")[0].str.split("=").apply(pd.Series)
        gene_names_fasta.columns = ["GN", "Gene name fasta"]
        # concat all important columns with the original dataframe
        df_protein_names = pd.concat([df_protein_names, fasta_col["protein id"], gene_names_fasta["Gene name fasta"]], axis=1)
        # add protein name from fasta description col
        df_protein_names["Protein name"] = fasta_col["description"].str.split("_", expand=True)[0]
        return df_protein_names, df_peptide_names

    def init_interest_from_xlsx(self) -> (dict, dict, dict):
        protein_path = os.path.join(self.path_pipeline_config, "important_protein_names.xlsx")
        receptor_path = os.path.join(self.path_pipeline_config, "important_receptor_names.xlsx")
        go_path = os.path.join(self.path_pipeline_config, "go_analysis_gene_names.xlsx")
        # make sure files exist
        if not os.path.isfile(protein_path):
            raise ValueError("missing important_protein_names.xlsx file")
        # make sure files exist
        if not os.path.isfile(receptor_path):
            raise ValueError("missing important_receptor_names.xlsx file")
        if not os.path.isfile(go_path):
            raise ValueError("missing go_analysis.xlsx file")

        def df_to_dict(df):
            return {col: df[col].dropna().tolist() for col in df}

        df_protein = pd.read_excel(protein_path)
        df_receptor = pd.read_excel(receptor_path)
        df_go = pd.read_excel(go_path)
        return df_to_dict(df_protein), df_to_dict(df_receptor), df_to_dict(df_go)

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
            plt.savefig(res_path, dpi=200, bbox_inches="tight")
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

        self.logger.info("Creating venn diagrams")
        file_dir_venn = os.path.join(self.start_dir, "venn")
        # create file structure and folder
        os.makedirs(file_dir_venn, exist_ok=True)

        # create all sets that are required for plotting
        protein_ids = ddict(dict)
        whole_experiment_protein_ids = {}
        for experiment in self.replicates:
            intensities = self.df_protein_names[
                [f"{self.configs['venn_intensity_col']}{rep}" for rep in self.replicates[experiment]]
            ]
            exp_prot_ids = set()
            for rep_name, rep in zip(self.replicates[experiment], intensities):
                # TODO what to use for venns
                rep_set = set(self.df_protein_names[intensities[rep] > 0]["Protein name"])
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
        experiment_sets = {exp: set.intersection(*(protein_ids[exp][rep] for rep in protein_ids[exp])) for exp in protein_ids}
        save_bar_venn("All experiments intersection", experiment_sets)
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

        all_intensities = self.df_protein_names[
            [f"{self.configs['descriptive_intensity_col']}{rep}"
             for exp in self.replicates for rep in self.replicates[exp]] +["Gene name fasta"]
        ].set_index("Gene name fasta")
        mask = (all_intensities != 0).sum(axis=1) != 0
        all_intensities = all_intensities[mask]
        ex_list = list(self.replicates.keys())
        for experiment, ax, ax1 in zip(self.replicates, axarr.flat, axarr1.flat):
            plt.close("all")
            e_list = [f"{self.configs['descriptive_intensity_col']}{rep}" for rep in self.replicates[experiment]]
            intensities = all_intensities[e_list]

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
            self.logger.debug(labels)
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
            found_proteins &= set(all_intensities.index)
            found_proteins = sorted(list(found_proteins))
            n_rows, n_cols = get_number_rows_cols_for_fig(found_proteins)
            fig, axarr = plt.subplots(n_rows, n_cols, figsize=(4 * n_rows, 4 * n_cols))
            fig.suptitle(pathway)
            for protein, ax in zip(found_proteins, axarr.flat):
                ax.set_title(protein)
                heights = []
                for idx, experiment in enumerate(self.replicates):
                    e_list = [x for x in all_intensities.columns if x.startswith(self.configs['descriptive_intensity_col'] + experiment)]
                    protein_intensities = all_intensities.loc[protein, e_list]
                    ax.scatter([idx] * len(protein_intensities), protein_intensities, label=f"{experiment}")
                    # self.logger.debug(f"{protein}: min: {min_i}, max: {max_i}")
                    heights.append(np.max(protein_intensities))
                # ax.set_ylim((0.9 * min_i, 1.1 * max_i))
                # ax.set_yscale("log")
                # ax.legend()
                n_annotations = 0
                for e1, e2 in combinations(self.replicates, 2):
                    e1_list = [x for x in all_intensities.columns if x.startswith(self.configs['descriptive_intensity_col'] + e1)]
                    e2_list = [x for x in all_intensities.columns if x.startswith(self.configs['descriptive_intensity_col'] + e2)]
                    v1 = all_intensities.loc[protein, e1_list].astype(float)
                    v2 = all_intensities.loc[protein, e2_list].astype(float)
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
            protein_intensities_ex1 = self.df_protein_names[
                [f"{self.configs['descriptive_intensity_col']}{rep}" for rep in self.replicates[ex1]]
            ]
            counts_ex1 = (protein_intensities_ex1 > 0).sum(axis=1) == len(protein_intensities_ex1.columns)
            protein_intensities_ex2 = self.df_protein_names[
                [f"{self.configs['descriptive_intensity_col']}{rep}" for rep in self.replicates[ex2]]
            ]
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

        all_intensities = self.df_protein_names[
            [f"{self.configs['descriptive_intensity_col']}{rep}"
             for exp in self.replicates for rep in self.replicates[exp]] +["Gene name fasta"]
        ].set_index("Gene name fasta")
        # eliminate all rows where intensities are 0
        mask = (all_intensities != 0).sum(axis=1) != 0
        # all proteins that were detected in any replicate
        background = set(all_intensities[mask].index)
        fig, ax = plt.subplots(1, 1, figsize=(7, 7))
        heights = ddict(list)
        test_results = ddict(list)

        bar_width = 0.25
        for compartiment, all_pathway_genes in self.go_analysis_gene_names.items():
            all_pathway_genes = set(all_pathway_genes)
            # 
            pathway_genes = all_pathway_genes & background
            # all genes that are not in the pathway
            not_pathway_genes = background - pathway_genes
            # sanity check
            assert pathway_genes | not_pathway_genes == background
            heights["background"].append(len(pathway_genes))
            for experiment in self.replicates:
                # create df with intensity means for specific experiment ofer all replicates
                abc = all_intensities[[x for x in all_intensities.columns if x.startswith('Intensity ' + experiment)]].mean(axis=1)
                # get all proteins with mean intensity > 0
                experiment_genes = set(abc[abc > 0].index)
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

    def update_config_file(self):
        # store the config file as tmp
        yml_file_loc_tmp = os.path.join(self.path_config, self.yml_file_name_tmp)
        with open(yml_file_loc_tmp, "w") as outfile:
            self.yaml.dump(self.configs, outfile)

        # delete non tmp if exists
        yml_file_loc = os.path.join(self.path_config, self.yml_file_name)
        if self.yml_file_name in os.listdir(self.path_config):
            os.remove(yml_file_loc)

        # rename to non tmp
        os.rename(yml_file_loc_tmp, yml_file_loc)


if __name__ == "__main__":
    # create arg parser to handle different options
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        '--dir',
        dest='dir',
        action='store',
        help="Path to directory of analysis which should be analyzed."
             "If not set the program will open a prompt and ask to select one."
    )
    parser.add_argument(
        '--yml-file',
        dest='yml_file',
        action='store',
        help="Path to yml file which should be used for analysis, or 'default'."
             "If not set the program will open a prompt and ask to select one,"
             "but if canceled the default yml file will be used."
    )
    args = parser.parse_args()
    args_dict = vars(args)
    # args.dir = "/media/b200-linux/Elements/max_quant_txt_results/xaiomeng_combined_buffers_no_isoform/"
    # print(args)
    # disable most ui things
    root = tk.Tk()
    root.withdraw()
    # if no dir was specified ask for one
    if args.dir is None:
        start_dir = filedialog.askdirectory()
        if not start_dir:
            raise ValueError("Please select a directory")
    else:
        start_dir = args.dir
    mxpipeline = MaxQuantPipeline(start_dir, args.yml_file)
    mxpipeline.create_results()
