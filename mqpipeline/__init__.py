from .version import __version__
from .MQInitializer import MQInitializer
from .MQPlots import MQPlots
from .MQPipeline import MQParser, MQUI
from .DataStructure import DataNode, DataTree

# import for "from package import *"
__all__ = [
    "MQPlots",
    "MQInitializer",
    "MQParser",
    "MQUI",
    "DataNode",
    "DataTree"
]

if __name__ == "__main__":
    mqparser = MQParser()
    gui = MQUI(**mqparser.args_dict)
