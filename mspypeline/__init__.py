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

if __name__ == "__main__":
    msparser = MSPParser()
    UIHandler(**msparser.args_dict)
