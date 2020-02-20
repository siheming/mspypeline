from .version import __version__
import os
# import for "from package import *"
__all__ = [
    "create_app",
    "path_package",
    "path_package_config",
]

path_package = script_loc = os.path.dirname(os.path.realpath(__file__))
path_package_config = os.path.join(path_package, "config")


def create_app(test_config=None):
    from mqpipeline.flask_scripts import create_app_helper
    return create_app_helper(test_config)


if __name__ == "__main__":
    from mqpipeline.core import MQParser, UIHandler
    mqparser = MQParser()
    UIHandler(**mqparser.args_dict)
