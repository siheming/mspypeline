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
from .plotting_backend import matplotlib_plots
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
    "matplotlib_plots",
    "MQReader"
]
__all__.extend(core.__all__)
__all__.extend(modules.__all__)
__all__.extend(file_reader.__all__)


def load_example_dataset(
        path_dataset: str = None, configs: dict = None, download: bool = False, remove_on_fail: bool = True
) -> MSPInitializer:
    """
    Loads the example dataset into a MSPInitializer. If download is true the data is downloaded from github.

    Parameters
    ----------
    path_dataset
        path to where the dataset should be read from and / or downloaded to
    configs
        configs that will be passed to the initializer configs
    download
        should the data be downloaded
    remove_on_fail
        should the passed path_dataset directory be removed if the download fails

    Returns
    -------
    An initializer with all data loaded into memory.

    """
    import pandas as pd

    if configs is None:
        configs = {}
    if path_dataset is None:
        path_dataset = os.path.join(os.path.split(path_package)[0], "mspypeline_public_dataset")
    path_dataset_txt = os.path.join(path_dataset, "txt")

    file_list = ["contaminants", "evidence", "msmsScans", "msScans", "parameters", "peptides", "proteinGroups",
                 "summary"]
    if download:
        assert int(pd.__version__.split(".")[0]) >= 1, "Pandas Version must be at least 1.0 to download"
        try:
            url_base = "https://raw.githubusercontent.com/siheming/mspypeline_public_dataset/main/txt/"
            os.makedirs(path_dataset_txt, exist_ok=True)
            for file in file_list:
                pd.read_pickle(url_base + file + ".gzip").to_pickle(os.path.join(path_dataset_txt, file + ".gzip"))
            pd.read_csv(url_base + "proteinGroups.txt", sep="\t"
                        ).to_csv(os.path.join(path_dataset_txt, "proteinGroups.txt"), sep="\t")
        except Exception as e:
            if remove_on_fail:
                import shutil
                shutil.rmtree(path_dataset_txt)
            raise e

    init = MSPInitializer(path_dataset)
    init.init_config()
    init.file_path_yaml = "default"
    init.configs.update(configs)
    init.read_data()

    try:
        for file in file_list:
            init.reader_data["mqreader"][file] = pd.read_pickle(os.path.join(path_dataset_txt, file + ".gzip"))
    except FileNotFoundError:
        raise FileNotFoundError("The files are not under the specified directory. Try specifying download=True")
    return init


if __name__ == "__main__":
    msparser = MSPParser()
    UIHandler(**msparser.args_dict)
