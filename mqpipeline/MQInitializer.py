import os
try:
    from ruamel_yaml import YAML
except ModuleNotFoundError:
    from ruamel.yaml import YAML
import logging

from mqpipeline.helpers import get_logger
from mqpipeline.file_reader import MQReader


class MQInitializer:
    # set all file names that are required
    yml_file_name_tmp = "config_tmp.yml"
    yml_file_name = "config.yml"
    default_yml_name = "ms_analysis_default.yml"
    go_path = "go_terms"
    pathway_path = "pathways"
    script_loc = os.path.dirname(os.path.realpath(__file__))
    path_pipeline_config = os.path.join(script_loc, "config")
    possible_gos = sorted([x for x in os.listdir(os.path.join(path_pipeline_config, go_path)) if x.endswith(".txt")])
    possible_pathways = sorted([x for x in os.listdir(os.path.join(path_pipeline_config, pathway_path)) if x.endswith(".txt")])

    def __init__(self, dir_: str, file_path_yml: str = "default", loglevel=logging.DEBUG):
        self.logger = get_logger(self.__class__.__name__, loglevel=loglevel)
        # create a yaml file reader
        self.yaml = YAML()
        # self.yaml.indent(mapping=2, sequence=4, offset=2)
        self.yaml.indent(offset=2)
        self.yaml.default_flow_style = False

        # attributes that change upon changing the starting dir
        self.analysis_design = None
        self.all_replicates = None
        self.configs = None
        self.path_config = None
        self.naming_convention = None

        self.df_protein_names, self.df_peptide_names = None, None
        self.interesting_proteins, self.go_analysis_gene_names = None, None

        self.logger.debug("Script location %s", MQInitializer.script_loc)
        # properties
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
        self.analysis_design = None
        self.all_replicates = None
        self.naming_convention = None
        self.configs = None
        self.file_path_yaml = "file"
        self.path_config = os.path.join(self.start_dir, "config")

    @property
    def file_path_yaml(self):
        return self._file_path_yaml

    @file_path_yaml.setter
    def file_path_yaml(self, file_path_yml: str):
        """

        Parameters
        ----------
        file_path_yml

        Raises
        ------
        ValueError
            if no valid value was provided
        FileNotFoundError
            if the file specified by the file_path_yml was not found

        """
        if file_path_yml.lower() == "default":
            self._file_path_yaml = self.get_default_yml_path()
        elif file_path_yml.lower() == "file":
            if self.has_yml_file():
                self._file_path_yaml = os.path.join(self.start_dir, "config", MQInitializer.yml_file_name)
            else:
                self._file_path_yaml = self.get_default_yml_path()
        elif file_path_yml.lower().endswith(('.yml', '.yaml')):
            self._file_path_yaml = os.path.normpath(file_path_yml)
        else:
            raise ValueError(f"Invalid value provided for yaml file: {file_path_yml}")
        self.logger.debug("yml file location: %s", self._file_path_yaml)

        # load the config from the yml file
        self.logger.info("loading yml file")
        with open(self.file_path_yaml) as f:
            self.configs = self.yaml.load(f)
        self.logger.debug(f"Config file contents: {self.configs}")

    def init_config(self):
        os.makedirs(self.path_config, exist_ok=True)
        self.update_config_file()

    def has_yml_file(self) -> bool:
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
        self.logger.debug("Loading default yml file from: %s, since no file was selected", MQInitializer.script_loc)
        if MQInitializer.default_yml_name in os.listdir(MQInitializer.path_pipeline_config):
            yaml_file = os.path.join(MQInitializer.path_pipeline_config, MQInitializer.default_yml_name)
        else:
            raise ValueError("Could not find default yaml file. Please select one.")
        return yaml_file

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
        fullpath = os.path.join(MQInitializer.path_pipeline_config, path, file)
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
        self.logger.info("Reading %s, and %s", MQReader.proteins_txt, MQReader.peptides_txt)
        # TODO make this work properly
        mqreader = MQReader(self.start_dir, self.configs)
        self.df_protein_names = mqreader.full_data[MQReader.proteins_txt]
        self.df_peptide_names = mqreader.full_data[MQReader.peptides_txt]
        self.all_replicates = mqreader.all_samples
        self.analysis_design = mqreader.analysis_design
        self.configs = mqreader.configs

        # read all proteins and receptors of interest from the config dir
        self.logger.info("Reading proteins and receptors of interest")
        self.interesting_proteins, self.go_analysis_gene_names = self.init_interest_from_txt()
        self.update_config_file()
