# setup variables required for package
import os
path_package = script_loc = os.path.dirname(os.path.realpath(__file__))
path_package_config = os.path.join(path_package, "config")


# this function is located here to host the flask sever
def create_app(test_config=None):
    from mqpipeline.flask_scripts import create_app_helper
    return create_app_helper(test_config)

# flatten package imports for the core package
from .version import __version__
from .core import MQUI, MQParser, MQPlots, UIHandler, MQInitializer
from .plotter import plotly_plots
# import for "from package import *"
__all__ = [
    "create_app",
    "path_package",
    "path_package_config",
    "MQInitializer",
    "MQUI",
    "MQParser",
    "UIHandler",
    "MQPlots"
]
# make sure everything that was imported is a class
import inspect
assert inspect.isclass(MQPlots)
assert inspect.isclass(MQInitializer)
assert inspect.isclass(MQUI)
assert inspect.isclass(MQParser)
assert inspect.isclass(UIHandler)

if __name__ == "__main__":
    mqparser = MQParser()
    UIHandler(**mqparser.args_dict)
