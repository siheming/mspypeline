from .version import __version__
from .MSPInitializer import MSPInitializer
from .MSPPlots import MSPPlots
from .MSPypeline import MSPParser, MSPUI

# import for "from package import *"
__all__ = [
    "MSPPlots",
    "MSPInitializer",
    "MSPParser",
    "MSPUI"
]

if __name__ == "__main__":
    mqparser = MSPParser()
    gui = MSPUI(**mqparser.args_dict)
