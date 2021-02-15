# setup variables required for package
import os
path_package = os.path.dirname(os.path.realpath(__file__))
path_package_config = os.path.join(path_package, "config")


# this function is located here to host the flask sever
def create_app(test_config=None):
    from mspypeline.flask_scripts import create_app_helper
    return create_app_helper(test_config)


# flatten package imports for the core package
from .version import __version__
from .plotting_backend import plotly_plots, matplotlib_plots
from .modules import *
from .core import *
from .file_reader import *
from .file_reader.MQReader import MQReader
# import for "from package import *"
__all__ = [
    "create_app",
    "path_package",
    "path_package_config",
    "__version__",
    "plotly_plots",
    "matplotlib_plots",
    "MQReader"
]
__all__.extend(core.__all__)
__all__.extend(modules.__all__)
__all__.extend(file_reader.__all__)


def load_example_dataset(path_dataset: str = None, configs: dict = None):
    import pandas as pd

    if configs is None:
        configs = {}
    if path_dataset is None:
        path_dataset = os.path.join(os.path.split(path_package)[0], "mspypeline_public_dataset")
    path_dataset_txt = os.path.join(path_dataset, "txt")
    init = MSPInitializer(path_dataset)
    init.init_config()
    init.file_path_yaml = "default"
    init.configs.update(configs)
    init.read_data()

    init.reader_data["mqreader"]["contaminants"] = pd.read_pickle(os.path.join(path_dataset_txt, "contaminants.gzip"))
    init.reader_data["mqreader"]["evidence"] = pd.read_pickle(os.path.join(path_dataset_txt, "evidence.gzip"))
    init.reader_data["mqreader"]["msmsScans"] = pd.read_pickle(os.path.join(path_dataset_txt, "msmsScans.gzip"))
    init.reader_data["mqreader"]["msScans"] = pd.read_pickle(os.path.join(path_dataset_txt, "msScans.gzip"))
    init.reader_data["mqreader"]["parameters"] = pd.read_pickle(os.path.join(path_dataset_txt, "parameters.gzip"))
    init.reader_data["mqreader"]["peptides"] = pd.read_pickle(os.path.join(path_dataset_txt, "peptides.gzip"))
    init.reader_data["mqreader"]["proteinGroups"] = pd.read_pickle(os.path.join(path_dataset_txt, "proteinGroups.gzip"))
    init.reader_data["mqreader"]["summary"] = pd.read_pickle(os.path.join(path_dataset_txt, "summary.gzip"))
    return init


if __name__ == "__main__":
    msparser = MSPParser()
    UIHandler(**msparser.args_dict)
