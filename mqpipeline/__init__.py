from .version import __version__
from .MQInitializer import MQInitializer
from .MQPlots import MQPlots
from .MQPipeline import MQParser, MQUI

# import for "from package import *"
__all__ = [
    "MQPlots",
    "MQInitializer",
    "MQParser",
    "MQUI"
]

if __name__ == "__main__":
    mqparser = MQParser()
    gui = MQUI(**mqparser.args_dict)
