# setup variables required for package
import os
path_package = script_loc = os.path.dirname(os.path.realpath(__file__))
path_package_config = os.path.join(path_package, "config")


# this function is located here to host the flask sever
def create_app(test_config=None):
    from mspypeline.flask_scripts import create_app_helper
    return create_app_helper(test_config)

# flatten package imports for the core package
from .version import __version__
from .core.MSPInitializer import MSPInitializer
from .core.MSPPlots import MSPPlots
from .core.MSPypeline import MSPUI, MSPParser, UIHandler
from .plotter import plotly_plots
# import for "from package import *"
__all__ = [
    "create_app",
    "path_package",
    "path_package_config",
    "MSPInitializer",
    "MSPUI",
    "MSPParser",
    "UIHandler",
    "MSPPlots"
]
# make sure everything that was imported is a class
import inspect
assert inspect.isclass(MSPPlots)
assert inspect.isclass(MSPInitializer)
assert inspect.isclass(MSPUI)
assert inspect.isclass(MSPParser)
assert inspect.isclass(UIHandler)

if __name__ == "__main__":
    msparser = MSPParser()
    UIHandler(**msparser.args_dict)
