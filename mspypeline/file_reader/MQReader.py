import os
import pandas as pd
from pandas.api.types import is_numeric_dtype
import logging
from collections import defaultdict as ddict

from mspypeline.helpers import get_overlap, string_similarity_ratio, dict_depth
from mspypeline.file_reader import BaseReader, MissingFilesException


class MQReader(BaseReader):
    proteins_txt = "proteinGroups.txt"
    peptides_txt = "peptides.txt"
    mapping_txt = "sample_mapping.txt"
    summary_txt = "summary.txt"
    parameters_txt = "parameters.txt"
    all_files = [proteins_txt, peptides_txt, summary_txt, parameters_txt]
    name = "mqreader"

    def __init__(self, start_dir: str,
                 reader_config: dict,
                 index_col: str = "Gene name",
                 duplicate_handling: str = "sum",
                 loglevel=logging.DEBUG):
        super().__init__(start_dir, reader_config, loglevel=loglevel)
        # TODO connect this to the configs of the initializer
        self.data_dir = os.path.join(self.start_dir, "txt")
        self.index_col = index_col
        self.duplicate_handling = duplicate_handling

        # start by sampling with all files
        # only save the dataframe if the file was found
        self.sample_data = {}
        for file in MQReader.all_files:
            try:
                file_dir = os.path.join(self.data_dir, file)
                df = pd.read_csv(file_dir, sep="\t", nrows=5)
                self.sample_data[file] = df
            except FileNotFoundError:
                self.logger.info("Did not find file: %s", file)

        # if no files were found throw an error
        if not self.sample_data:
            raise MissingFilesException("Could find none of: " + ", ".join(MQReader.all_files))

        # read the sample name mapping file
        try:
            self.mapping_txt = pd.read_csv(os.path.join(self.start_dir, MQReader.mapping_txt),
                                           sep="\t", header=0, index_col=None)
            if self.mapping_txt.shape[1] != 2:
                raise ValueError(f"{MQReader.mapping_txt} should have two columns")
            duplicated_old = self.mapping_txt.iloc[:, 0].duplicated(keep=False)
            duplicated_new = self.mapping_txt.iloc[:, 1].duplicated(keep=False)
            if any(duplicated_new) or any(duplicated_old):
                raise ValueError(f"{MQReader.mapping_txt} should only contain unique rows, "
                                 f"{self.mapping_txt.iloc[:, 0][duplicated_old]}, {self.mapping_txt.iloc[:, 1][duplicated_new]}")
        except FileNotFoundError:
            self.mapping_txt = None

        # rename all columns based on the mapping
        self.logger.info("Renaming columns using %s", MQReader.mapping_txt)
        self.new_column_names = {file: df.columns for file, df in self.sample_data.items()}
        if self.mapping_txt is not None:
            for file in self.sample_data:
                self.new_column_names[file] = self.rename_df_columns(self.sample_data[file].columns)

        # get columns that should be dropped
        to_drop = self.reader_config.get("drop_columns", [])
        if to_drop:
            if not isinstance(to_drop, list):
                to_drop = [to_drop]

        # subset on all columns that start with intensity
        self.intensity_column_names = {}
        for file, column_name_list in self.new_column_names.items():
            col_names = sorted([x.replace("Intensity ", "") for x in column_name_list if x.startswith('Intensity ')],
                               key=len, reverse=True)
            # drop all columns
            col_names = [x for x in col_names if x not in to_drop]
            self.intensity_column_names[file] = col_names
        self.column_name_sample = self.intensity_column_names[list(self.intensity_column_names.keys())[0]]
        if not self.reader_config.get("all_replicates", False):
            self.reader_config["all_replicates"] = self.column_name_sample
        # check the naming convention
        self.naming_convention = self.check_naming_convention()
        # determine the grouping
        self.analysis_design = self.determine_groupings()
        if not self.reader_config.get("analysis_design", False):
            self.reader_config["analysis_design"] = self.analysis_design
            self.reader_config["levels"] = dict_depth(self.analysis_design)
        # update the config
        #self.update_config_file()  # TODO write this to the config file

        # TODO turning this into a property and only preprocessing when required might be a lot better
        self.full_data = {}
        for file in MQReader.all_files:
            if file in self.sample_data:
                self.logger.debug("Preprocessing %s", file)
                full_data = getattr(self, f"preprocess_{file.replace('.txt', '')}")()
                self.full_data[file] = full_data

    def rename_df_columns(self, col_names: list) -> list:
        mapping_dict = {old_name: new_name for old_name, new_name in zip(self.mapping_txt.iloc[:, 0], self.mapping_txt.iloc[:, 1])}

        def find_key_match(col):
            matches = []
            for key in mapping_dict:
                if key in col:
                    matches.append(key)
            matches = sorted(matches, key=len)
            if len(matches) == 0:
                return ""
            else:
                return matches[-1]
        match_list = [find_key_match(col) for col in col_names]
        # add a replace value for the default string
        mapping_dict[""] = ""
        return [col_names[i].replace(match_list[i], mapping_dict[match_list[i]]) for i in range(len(col_names))]

    def check_naming_convention(self) -> bool:
        # does the name follow the convention
        first_rep_split = len(self.column_name_sample[0].split("_"))
        naming_convention = all(len(rep.split("_")) == first_rep_split for rep in self.column_name_sample)
        # create a mapping template if the convention is not followed and the file doesn't exist
        mapping_template_name = os.path.join(self.start_dir, MQReader.mapping_txt.replace("ing", "ing_template"))

        if not naming_convention:
            self.logger.warning("Naming of experiments does not follow naming convention, "
                                "please consider using a %s file", MQReader.mapping_txt)
            if os.path.isfile(mapping_template_name):
                self.logger.warning("Currently unused %s file in the directory", mapping_template_name)
            else:
                self.logger.warning("Creating %s file. Please provide a mapping that follows "
                                    "the indicated naming convention", mapping_template_name)
                old_sample_names = sorted([x.replace('Intensity ', '') for x in self.sample_data[list(self.sample_data.keys())[0]]
                        if x.startswith('Intensity ')], key=len, reverse=True)
                # dump the names that still correspond to the names in the txt files
                mapping = pd.DataFrame({"old name": old_sample_names,
                                        "new name": ["groupname_experimentname_techrepname"] * len(old_sample_names)})
                mapping.to_csv(mapping_template_name, header=True, index=False, sep="\t")
        return naming_convention

    def determine_groupings(self):
        for file, l in self.intensity_column_names.items():
            unmatched = [x for x in l if x not in self.column_name_sample]
            if len(unmatched) > 0:
                self.logger.warning("mismatch between files found")
                continue
                raise ValueError(
                    f"Found replicates in {MQReader.proteins_txt}, but not in {MQReader.peptides_txt}: " +
                    ", ".join(unmatched))

        # extract the analysis design from the file
        if not self.reader_config.get("analysis_design", False):
            # try to automatically determine experimental setup
            # if the naming convention is followed it is quite easy
            if self.naming_convention:
                factory = lambda: ddict(factory)
                analysis_design = factory()

                def fill_dict(d: dict, s: str, s_split=None):
                    if s_split is None:
                        s_split = s.split("_")
                    if len(s_split) > 1:
                        fill_dict(d[s_split[0]], s, s_split[1:])
                    else:
                        d[s_split[0]] = s

                def default_to_regular(d: dict):
                    if isinstance(d, ddict):
                        d = {k: default_to_regular(v) for k, v in d.items()}
                    return d

                for name in self.column_name_sample:
                    fill_dict(analysis_design, name)

                analysis_design = default_to_regular(analysis_design)
            # otherwise we can just guess grouping
            else:
                analysis_design = self.guess_analysis_design(self.column_name_sample)
        else:
            analysis_design = self.reader_config.get("analysis_design")
        return analysis_design

    def guess_analysis_design(self, all_reps):
        raise NotImplementedError("This is not implemented at the moment. Please stick to the naming convention")
        # TODO update the attempted matching mechanism

    def preprocess_proteinGroups(self):
        file_dir = os.path.join(self.data_dir, MQReader.proteins_txt)
        df_protein_groups = pd.read_csv(file_dir, sep="\t")
        df_protein_groups.columns = self.new_column_names[MQReader.proteins_txt]
        not_contaminants = (df_protein_groups[
                                ["Only identified by site", "Reverse", "Potential contaminant"]] == "+"
                            ).sum(axis=1) == 0
        self.logger.debug("Removing %s rows from %s because they are marked as contaminant",
                          (~not_contaminants).sum(), MQReader.proteins_txt)
        df_protein_groups = df_protein_groups[not_contaminants]
        if any(df_protein_groups["Fasta headers"].isna()):
            self.logger.warning("Missing fasta headers using default columns for information")
            gene_name = df_protein_groups["Gene names"]
            sep_ind = gene_name.str.contains(";").fillna(False)
            gene_name[sep_ind] = gene_name[sep_ind].str.split(";", expand=True)[0]
            concat_df = pd.DataFrame({
                "protein id": df_protein_groups["Protein names"],
                "Gene name": gene_name.str.upper(),
                "Protein name": ["Missing"] * gene_name.shape[0],
            })
        else:
            # split the fasta headers
            colon_start = df_protein_groups["Fasta headers"].str.startswith(";")
            df_protein_groups.loc[colon_start, "Fasta headers"] = df_protein_groups.loc[
                colon_start, "Fasta headers"].str.lstrip(";")
            # first split all fasta headers that contain multiple entries
            sep_ind = df_protein_groups["Fasta headers"].str.contains(";").fillna(False)
            # replace all fasta headers with multiple entries with only the first one
            df_protein_groups.loc[sep_ind, "Fasta headers"] = \
                df_protein_groups.loc[sep_ind, "Fasta headers"].str.split(";", expand=True)[0]
            # split the fasta headers with the pipe symbol
            fasta_col = df_protein_groups["Fasta headers"].str.split("|", n=2).apply(pd.Series)
            fasta_col.columns = ["trash", "protein id", "description"]
            # extract the gene name from the description eg: "GN=abcd"
            gene_names_fasta = fasta_col["description"].str.extract(r"(GN=(.*?)(\s|$))")[1]
            # added upper() function to avoid that non-human gene names are not recognized
            concat_df = pd.DataFrame({
                "protein id": fasta_col["protein id"],
                "Gene name": gene_names_fasta.str.upper(),
                "Protein name": fasta_col["description"].str.split("_", expand=True)[0]
            })
        # concat all important columns with the original dataframe
        df_protein_groups = pd.concat([df_protein_groups, concat_df], axis=1)
        # remove all rows where the column used for indexing is missing
        mask = ~pd.isna(df_protein_groups[self.index_col])
        df_protein_groups = df_protein_groups.loc[mask]
        if ~mask.sum() > 0:
            self.logger.warning("Removing %s rows because the index col information from: %s is missing",
                                ~mask.sum(), self.index_col)
        # set index
        self.logger.info("Setting index of %s to %s", MQReader.proteins_txt, self.index_col)
        df_protein_groups = df_protein_groups.set_index(df_protein_groups[self.index_col], drop=False)
        # convert all non numeric intensities
        for col in [col for col in df_protein_groups.columns if "Intensity " in col or "LFQ " in col or "iBAQ " in col]:
            if not is_numeric_dtype(df_protein_groups[col]):
                df_protein_groups[col] = df_protein_groups[col].apply(lambda x: x.replace(",", ".")).fillna(0)
                df_protein_groups[col] = df_protein_groups[col].astype("int64")
        # handle all rows with duplicated index column
        duplicates = df_protein_groups.duplicated(subset=self.index_col, keep=False)
        if any(duplicates):
            self.logger.warning("Found duplicates in %s column. Duplicate names: %s",
                                self.index_col, ", ".join(df_protein_groups[duplicates].loc[:, self.index_col]))
            if self.duplicate_handling == "drop":
                self.logger.warning("Dropping all %s duplicates.", duplicates.sum())
                df_protein_groups = df_protein_groups.drop_duplicates(subset=self.index_col, keep=False)
            elif self.duplicate_handling == "sum":
                def group_sum(x):
                    x.iloc[0].loc[x.select_dtypes("number").columns] = x.sum(axis=0, numeric_only=True)
                    return x.iloc[0]

                new_index = df_protein_groups.index.drop_duplicates(keep=False)
                duplicate_index = df_protein_groups.index.difference(new_index)
                df_dup = df_protein_groups.loc[duplicate_index, :]
                self.logger.warning("Merging %s rows into %s by summing numerical columns."
                                    "Some information might be incorrect", df_dup.shape[0], duplicate_index.shape[0])
                df_dup = df_dup.groupby(df_dup.index).apply(group_sum)
                df_protein_groups = pd.concat([df_protein_groups.loc[new_index, :], df_dup], axis=0)
        self.logger.debug("%s shape after preprocessing: %s", MQReader.proteins_txt, df_protein_groups.shape)
        return df_protein_groups

    def preprocess_peptides(self):
        file_dir = os.path.join(self.data_dir, MQReader.peptides_txt)
        df_peptides = pd.read_csv(file_dir, sep="\t")
        df_peptides.columns = self.new_column_names[MQReader.peptides_txt]
        not_contaminants = (df_peptides[
                                ["Reverse", "Potential contaminant"]] == "+"
                            ).sum(axis=1) == 0
        df_peptides = df_peptides[not_contaminants]
        self.logger.debug("Removing %s rows from %s because they are marked as contaminant",
                          (~not_contaminants).sum(), MQReader.peptides_txt)

        return df_peptides

    def preprocess_summary(self):
        file_dir = os.path.join(self.data_dir, MQReader.summary_txt)
        df_summary = pd.read_csv(file_dir, sep="\t")
        df_summary.columns = self.new_column_names[MQReader.summary_txt]
        df_summary = df_summary[df_summary["Enzyme"] != ""]
        return df_summary

    def preprocess_parameters(self):
        file_dir = os.path.join(self.data_dir, MQReader.parameters_txt)
        df_parameters = pd.read_csv(file_dir, sep="\t", index_col=[0], squeeze=True)
        df_parameters.columns = self.new_column_names[MQReader.parameters_txt]
        return df_parameters
