import os
from typing import Union
import pandas as pd
from pandas.api.types import is_numeric_dtype
import logging

from mspypeline.helpers import dict_depth, get_analysis_design
from mspypeline.file_reader import BaseReader, MissingFilesException
from mspypeline.core import MaxQuantPlotter


class MQReader(BaseReader):
    proteins_txt = "proteinGroups.txt"
    peptides_txt = "peptides.txt"
    mapping_txt = "sample_mapping.txt"
    summary_txt = "summary.txt"
    parameters_txt = "parameters.txt"
    msms_scans_txt = "msmsScans.txt"
    ms_scans_txt = "msScans.txt"
    evidence_txt = "evidence.txt"
    required_files = [proteins_txt]
    name = "mqreader"
    plotter = MaxQuantPlotter

    def __init__(self, start_dir: str,
                 reader_config: dict,
                 index_col: str = "Gene name",
                 duplicate_handling: str = "sum",
                 drop_columns: Union[list, tuple, str] = None,
                 loglevel=logging.DEBUG):
        super().__init__(start_dir, reader_config, loglevel=loglevel)

        if os.path.split(self.start_dir)[1] != "txt":
            self.data_dir = os.path.join(self.start_dir, "txt")
        else:
            self.data_dir = self.start_dir
        self.index_col = self.reader_config.get("index_col", index_col)
        self.duplicate_handling = self.reader_config.get("duplicate_handling", duplicate_handling)
        to_drop = self.reader_config.get("drop_columns", drop_columns if drop_columns is not None else [])

        # read a sample of all required files. If any required file is missing exit
        # but we need only one file from the max quant results
        try:
            file_dir = os.path.join(self.data_dir, self.proteins_txt)
            df = pd.read_csv(file_dir, sep="\t", nrows=5)
            self.proteins_txt_columns = df.columns
        except FileNotFoundError:
            raise MissingFilesException("Could find all of: " + ", ".join(MQReader.required_files))

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
                                 f"{self.mapping_txt.iloc[:, 0][duplicated_old]}, "
                                 f"{self.mapping_txt.iloc[:, 1][duplicated_new]}")
            self.logger.info("Successfully loaded %s", MQReader.mapping_txt)
        except FileNotFoundError:
            self.mapping_txt = None

        # rename all columns based on the mapping
        self.new_proteins_txt_columns = self.proteins_txt_columns
        if self.mapping_txt is not None:
            self.new_proteins_txt_columns = self.rename_df_columns(self.new_proteins_txt_columns)

        # get columns that should be dropped
        if isinstance(to_drop, str):
            to_drop = [to_drop]

        # subset on all columns that start with intensity
        self.intensity_column_names = sorted([x.replace("Intensity ", "") for x in self.new_proteins_txt_columns
                                              if x.startswith('Intensity ')], key=len, reverse=True)
        self.intensity_column_names = [x for x in self.intensity_column_names if x not in to_drop]
        if not self.reader_config.get("all_replicates", False):
            self.reader_config["all_replicates"] = self.intensity_column_names
        # check the naming convention
        self.naming_convention = self.check_naming_convention()
        # determine the grouping
        self.analysis_design = self.determine_groupings()
        if not self.reader_config.get("analysis_design", False):
            self.reader_config["analysis_design"] = self.analysis_design
            self.reader_config["levels"] = dict_depth(self.analysis_design)
            self.reader_config["level_names"] = [x for x in range(self.reader_config["levels"])]

    def rename_df_columns(self, col_names: list) -> list:
        if self.mapping_txt is None:
            return col_names
        mapping_dict = {old_name: new_name for old_name, new_name
                        in zip(self.mapping_txt.iloc[:, 0], self.mapping_txt.iloc[:, 1])}

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
        first_rep_split = len(self.intensity_column_names[0].split("_"))
        naming_convention = all(len(rep.split("_")) == first_rep_split for rep in self.intensity_column_names)

        # create a mapping template if the convention is not followed and the file doesn't exist
        if not naming_convention:
            mapping_template_name = os.path.join(self.start_dir, MQReader.mapping_txt.replace("ing", "ing_template"))
            self.logger.warning("Naming of experiments does not follow naming convention, "
                                "please consider using a %s file", MQReader.mapping_txt)
            if os.path.isfile(mapping_template_name):
                self.logger.warning("Currently unused %s file in the directory", mapping_template_name)
            else:
                self.logger.warning("Creating %s file. Please provide a mapping that follows "
                                    "the indicated naming convention", mapping_template_name)
                old_sample_names = sorted([x.replace('Intensity ', '') for x in self.proteins_txt_columns
                                           if x.startswith('Intensity ')], key=len, reverse=True)
                # dump the names that still correspond to the names in the txt files
                mapping = pd.DataFrame({"old name": old_sample_names,
                                        "new name": ["groupname_experimentname_techrepname"] * len(old_sample_names)})
                mapping.to_csv(mapping_template_name, header=True, index=False, sep="\t")
        return naming_convention

    def determine_groupings(self):
        # extract the analysis design from the file
        if not self.reader_config.get("analysis_design", False):
            # try to automatically determine experimental setup
            # if the naming convention is followed it is quite easy
            if self.naming_convention:
                analysis_design = get_analysis_design(self.intensity_column_names)
            # otherwise we can just guess grouping
            else:
                analysis_design = self.guess_analysis_design(self.intensity_column_names)
        else:
            analysis_design = self.reader_config.get("analysis_design")
        return analysis_design

    def guess_analysis_design(self, all_reps):
        raise NotImplementedError("This is not implemented at the moment. Please stick to the naming convention")
        # TODO update the attempted matching mechanism

    def preprocess_proteinGroups(self):
        file_dir = os.path.join(self.data_dir, MQReader.proteins_txt)
        df_protein_groups = pd.read_csv(file_dir, sep="\t")
        df_protein_groups.columns = self.rename_df_columns(df_protein_groups.columns)
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
            if sep_ind.sum() > 0:
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
                self.logger.warning("Merging %s rows into %s by summing numerical columns. "
                                    "Some information might be incorrect", df_dup.shape[0], duplicate_index.shape[0])
                df_dup = df_dup.groupby(df_dup.index).apply(group_sum)
                df_protein_groups = pd.concat([df_protein_groups.loc[new_index, :], df_dup], axis=0)
        self.logger.debug("%s shape after preprocessing: %s", MQReader.proteins_txt, df_protein_groups.shape)
        return df_protein_groups

    def preprocess_peptides(self):
        file_dir = os.path.join(self.data_dir, MQReader.peptides_txt)
        df_peptides = pd.read_csv(file_dir, sep="\t")
        df_peptides.columns = self.rename_df_columns(df_peptides.columns)
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
        df_summary.columns = self.rename_df_columns(df_summary.columns)
        df_summary = df_summary[df_summary["Enzyme"].notna()]
        df_summary["Experiment"] = self.rename_df_columns(df_summary["Experiment"])
        return df_summary

    def preprocess_parameters(self):
        file_dir = os.path.join(self.data_dir, MQReader.parameters_txt)
        df_parameters = pd.read_csv(file_dir, sep="\t", index_col=[0], squeeze=True)
        return df_parameters

    def preprocess_evidence(self):
        file_dir = os.path.join(self.data_dir, MQReader.evidence_txt)
        df_evidence = pd.read_csv(file_dir, sep="\t")
        not_contaminants = (df_evidence[["Reverse", "Potential contaminant"]] == "+").sum(axis=1) == 0
        df_evidence = df_evidence[not_contaminants]
        df_evidence.columns = self.rename_df_columns(df_evidence.columns)
        df_evidence["Experiment"] = self.rename_df_columns(df_evidence["Experiment"].tolist())
        return df_evidence

    def preprocess_msScans(self):
        file_dir = os.path.join(self.data_dir, MQReader.ms_scans_txt)
        df_msscans = pd.read_csv(file_dir, sep="\t", index_col=[0],
                                 usecols=["Raw file", "Total ion current", "Retention time"])
        return df_msscans

    def preprocess_msmsScans(self):
        file_dir = os.path.join(self.data_dir, MQReader.msms_scans_txt)
        df_msmsscans = pd.read_csv(file_dir, sep="\t", index_col=[0],
                                   usecols=["Raw file", "Total ion current", "Retention time"])
        return df_msmsscans
