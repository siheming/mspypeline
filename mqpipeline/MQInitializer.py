import pandas as pd
from pandas.api.types import is_numeric_dtype
import os
from mqpipeline.Logger import Logger
from collections import defaultdict as ddict

from mqpipeline.Utils import get_overlap

try:
    from ruamel_yaml import YAML
except ModuleNotFoundError:
    from ruamel.yaml import YAML
from tkinter import filedialog
import logging


class MQInitializer(Logger):
    # set all file names that are required
    yml_file_name_tmp = "config_tmp.yml"
    yml_file_name = "config.yml"
    default_yml_name = "ms_analysis_default.yml"
    go_path = "go_terms"
    pathway_path = "pathways"
    proteins_txt = "proteinGroups.txt"
    peptides_txt = "peptides.txt"

    def __init__(self, dir_: str, file_path_yml: str = None, loglevel=logging.DEBUG):
        super().__init__(self.__class__.__name__, loglevel=loglevel)
        self.yaml = YAML()
        self.yaml.default_flow_style = False
        self.replicates = None
        self.configs = None
        self.path_config = None
        self.df_protein_names, self.df_peptide_names = None, None
        self.interesting_proteins, self.go_analysis_gene_names = None, None

        self.script_loc = os.path.dirname(os.path.realpath(__file__))
        self.path_pipeline_config = os.path.join(self.script_loc, "config")
        self.logger.debug("Script location %s", self.script_loc)

        self.possible_gos = sorted([x for x in os.listdir(os.path.join(self.path_pipeline_config , MQInitializer.go_path)) if x.endswith(".txt")])
        self.possible_pathways = sorted([x for x in os.listdir(os.path.join(self.path_pipeline_config , MQInitializer.pathway_path)) if x.endswith(".txt")])

        self._start_dir = None
        self.start_dir = dir_

        self._file_path_yaml = None
        self.file_path_yaml = file_path_yml

    @property
    def start_dir(self):
        return self._start_dir

    @start_dir.setter
    def start_dir(self, start_dir):
        start_dir = os.path.normpath(start_dir)
        # make sure to be on the right level and set starting dir
        if os.path.split(start_dir)[1] == "txt":
            self.logger.debug("Removing txt ending from path")
            self._start_dir = os.path.split(start_dir)[0]
        else:
            self._start_dir = start_dir
        self.logger.info(f"Starting dir: {self.start_dir}")
        self.replicates = None
        self.path_config = os.path.join(self.start_dir, "config")

    @property
    def file_path_yaml(self):
        return self._file_path_yaml

    @file_path_yaml.setter
    def file_path_yaml(self, file_path_yml):
        # if no yml file is passed try to guess it or ask for one
        if file_path_yml.lower() == "default":
            self._file_path_yaml = self.get_default_yml_path()
        elif file_path_yml.lower() == "file":
            self._file_path_yaml = os.path.join(self.start_dir, "config", MQInitializer.yml_file_name)
        elif file_path_yml is None:
            self._file_path_yaml = self.init_yml_path()
        self.logger.debug(f"yml file location: {self._file_path_yaml}")

        # load the config from the yml file
        self.logger.info("loading yml file")
        with open(self.file_path_yaml) as f:
            self.configs = self.yaml.load(f)
        self.logger.debug(f"Config file contents: {self.configs}")

    def init_config(self):
        os.makedirs(self.path_config, exist_ok=True)
        self.update_config_file()

    def has_yml_file(self):
        if not self.start_dir:
            return False
        if "config" in os.listdir(self.start_dir):
            self.logger.debug("Found config dir")
            config_dir = os.path.join(self.start_dir, "config")
            if MQInitializer.yml_file_name in os.listdir(config_dir):
                self.logger.debug("Found config.yml file in config dir")
                return True
        return False

    def get_default_yml_path(self) -> str:
        self.logger.debug("Loading default yml file from: %s, since no file was selected", self.script_loc)
        if MQInitializer.default_yml_name in os.listdir(self.path_pipeline_config):
            yaml_file = os.path.join(self.path_pipeline_config, MQInitializer.default_yml_name)
        else:
            raise ValueError("Could not find default yaml file. Please select one.")
        return yaml_file

    def init_yml_path(self) -> str:
        def yml_file_dialog() -> str:
            file_path = filedialog.askopenfilename(filetypes=[(".yaml", ".yml")],
                                                   title='Please select a yml / yaml settings file')
            self.logger.debug(f"selected file path: {file_path}")
            if not file_path:
                yaml_file = self.get_default_yml_path()
            else:
                yaml_file = file_path
            return yaml_file

        if self.has_yml_file():
            yaml_file = os.path.join(self.start_dir, "config", MQInitializer.yml_file_name)
        else:
            yaml_file = yml_file_dialog()
        return yaml_file

    def init_dfs_from_txts(self) -> (pd.DataFrame, pd.DataFrame):
        file_dir_txt = os.path.join(self.start_dir, "txt")
        if not os.path.isdir(file_dir_txt):
            raise ValueError("Directory does not contain a txt dir")
        file_dir_protein_names = os.path.join(file_dir_txt, MQInitializer.proteins_txt)
        file_dir_peptides_names = os.path.join(file_dir_txt, MQInitializer.peptides_txt)
        # make sure protein groups file exists
        if not os.path.isfile(file_dir_protein_names):
            raise ValueError(f"txt directory does not contain a {MQInitializer.proteins_txt} file")
        if os.path.isfile(file_dir_peptides_names):
            df_peptide_names = pd.read_csv(file_dir_peptides_names, sep="\t")
        else:
            df_peptide_names = pd.DataFrame()
        # read protein groups file
        df_protein_names = pd.read_csv(file_dir_protein_names, sep="\t")
        self.logger.debug("%s shape: %s", MQInitializer.proteins_txt, df_protein_names.shape)

        # TODO can the Intensity column always be expected in the file?
        # TODO will the column names always be the same between Intensity and LFQ intensity?
        all_reps = sorted([x.replace('Intensity ', '') for x in df_protein_names.columns
                           if x.startswith('Intensity ')], key=len, reverse=True)
        # make sure the two files contain the same replicate names
        all_reps_peptides = [x.replace('Intensity ', '') for x in df_peptide_names.columns
                             if x.startswith('Intensity ')]
        if df_peptide_names.shape != (0, 0):
            experiment_analysis_overlap = [x not in all_reps for x in all_reps_peptides]
            if any(experiment_analysis_overlap):
                unmatched = [x for x in all_reps_peptides if experiment_analysis_overlap]
                raise ValueError(f"Found replicates in {MQInitializer.proteins_txt}, but not in {MQInitializer.peptides_txt}: " +
                                 ", ".join(unmatched))

        # try to automatically determine experimental setup
        if not self.configs.get("experiments", False):
            self.logger.info("No experiments specified in settings file. Trying to infer.")

            if self.configs.get("has_replicates", True):
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
                        self.logger.warning(f"unclear match for experiment: {experiment}")
                self.logger.info(f"determined experiments: {replicates.keys()}")
                self.logger.debug("number of replicates per experiment:")
                self.logger.debug("\n".join([f"{ex}: {len(replicates[ex])}" for ex in replicates]))
                self.configs["experiments"] = list(replicates.keys())
                self.replicates = replicates
            else:
                self.configs["experiments"] = all_reps
                self.replicates = {rep: [rep] for rep in all_reps}
        else:
            if self.configs.get("has_replicates", True):
                self.logger.info("Using saved experimental setup")
                replicates = ddict(list)
                for rep in all_reps:
                    overlaps = pd.Series(
                        [len(rep.replace(experiment, "")) if len(rep) > len(rep.replace(experiment, "")) else -1
                         for experiment in self.configs.get("experiments")]
                    )
                    overlap_minimum = overlaps[overlaps >= 0].idxmin()
                    replicates[self.configs.get("experiments")[overlap_minimum]].append(rep)
                self.replicates = replicates
            else:
                self.replicates = {rep: [rep] for rep in all_reps}

        # filter all contaminants by removing all rows where any of the 3 columns contains a +
        not_contaminants = (df_protein_names[
                                ["Only identified by site", "Reverse", "Potential contaminant"]] == "+"
                            ).sum(axis=1) == 0
        self.logger.debug("Removing %s rows from %s because they are marked as contaminant",
                          (~not_contaminants).sum(), MQInitializer.proteins_txt)
        df_protein_names = df_protein_names[not_contaminants]
        # split the fasta headers
        colon_start = df_protein_names["Fasta headers"].str.startswith(";")
        df_protein_names.loc[colon_start, "Fasta headers"] = df_protein_names.loc[colon_start, "Fasta headers"].str.lstrip(";")
        # first split all fasta headers that contain multiple entries
        sep_ind = df_protein_names["Fasta headers"].str.contains(";").fillna(False)
        # replace all fasta headers with multiple entries with only the first one
        # TODO will there always be a fasta header?
        df_protein_names.loc[sep_ind, "Fasta headers"] = df_protein_names.loc[sep_ind, "Fasta headers"].str.split(";").apply(pd.Series)[0]
        # split the fasta headers with the pipe symbol
        fasta_col = df_protein_names["Fasta headers"].str.split("|", n=2).apply(pd.Series)
        fasta_col.columns = ["trash", "protein id", "description"]
        # extract the gene name from the description eg: "GN=abcd"
        gene_names_fasta = fasta_col["description"].str.extract(r"(GN=(.*?)(\s|$))")[1]
        gene_names_fasta.name = "Gene name fasta"
        # remove all rows without gene names
        mask = ~pd.isna(gene_names_fasta)
        gene_names_fasta = gene_names_fasta.loc[mask]
        fasta_col = fasta_col.loc[mask]
        df_protein_names = df_protein_names.loc[mask]
        # added upper() function to avoid that non-human gene names are not recognized
        gene_names_fasta = gene_names_fasta.str.upper()
        # concat all important columns with the original dataframe
        df_protein_names = pd.concat([df_protein_names, fasta_col["protein id"], gene_names_fasta], axis=1)
        # add protein name from fasta description col
        df_protein_names["Protein name"] = fasta_col["description"].str.split("_", expand=True)[0]
        # filter all entries with duplicate Gene name fasta
        duplicates = df_protein_names.duplicated(subset="Gene name fasta", keep=False)
        if any(duplicates):
            self.logger.warning("Found duplicates in Gene name fasta column. Dropping all %s duplicates. Duplicate Gene names: %s",
                                duplicates.sum(), ", ".join(df_protein_names[duplicates].loc[:, "Gene name fasta"]))
            df_protein_names = df_protein_names.drop_duplicates(subset="Gene name fasta", keep=False)
        # convert all non numeric intensities
        for col in [col for col in df_protein_names.columns if 'Intensity ' in col or "LFQ " in col or "iBAQ " in col]:
            if not is_numeric_dtype(df_protein_names[col]):
                df_protein_names[col] = df_protein_names[col].apply(lambda x: x.replace(",", ".")).fillna(0)
                df_protein_names[col] = df_protein_names[col].astype("int64")
        self.logger.debug("%s shape after preprocessing: %s", MQInitializer.proteins_txt, df_protein_names.shape)
        return df_protein_names, df_peptide_names

    def init_interest_from_txt(self):
        dict_pathway = {}
        dict_go = {}
        for pathway in self.configs.get("pathways"):
            name, proteins = self.read_config_txt_file(MQInitializer.pathway_path, pathway)
            dict_pathway[name] = proteins

        for go in self.configs.get("go_terms"):
            name, proteins = self.read_config_txt_file(MQInitializer.go_path, go)
            dict_go[name] = proteins
        return dict_pathway, dict_go

    def read_config_txt_file(self, path, file):
        fullpath = os.path.join(self.path_pipeline_config, path, file)
        if path == MQInitializer.pathway_path:
            with open(fullpath) as f:
                name = f.readline().strip()
                f.readline()
                proteins = []
                for line in f:
                    proteins.append(line.strip())
        elif path == MQInitializer.go_path:
            name = file.replace(".txt", "")
            with open(fullpath) as f:
                proteins = []
                for line in f:
                    proteins.append(line.strip())
        else:
            raise ValueError(f"Invalid path: {path}")
        return name, proteins

    def update_config_file(self):
        # store the config file as tmp
        self.logger.debug("Updating yml settings file")
        yml_file_loc_tmp = os.path.join(self.path_config, MQInitializer.yml_file_name_tmp)
        with open(yml_file_loc_tmp, "w") as outfile:
            self.yaml.dump(self.configs, outfile)

        # delete non tmp if exists
        yml_file_loc = os.path.join(self.path_config, MQInitializer.yml_file_name)
        if MQInitializer.yml_file_name in os.listdir(self.path_config):
            os.remove(yml_file_loc)

        # rename to non tmp
        os.rename(yml_file_loc_tmp, yml_file_loc)

    def prepare_stuff(self):
        # read the proteins_txt and peptides_txt
        self.logger.info("Reading %s, and %s", MQInitializer.proteins_txt, MQInitializer.peptides_txt)
        self.df_protein_names, self.df_peptide_names = self.init_dfs_from_txts()

        # read all proteins and receptors of interest from the config dir
        self.logger.info("Reading proteins and receptors of interest")
        self.interesting_proteins, self.go_analysis_gene_names = self.init_interest_from_txt()
        self.update_config_file()
