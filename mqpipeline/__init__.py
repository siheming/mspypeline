from .version import __version__
from .MQInitializer import MQInitializer
from .MQPlots import MQPlots
from .MQPipeline import MQParser, MQGUI

# import for "from package import *"
__all__ = [
    "MQPlots",
    "MQInitializer",
    "MQParser",
    "MQGUI"
]

if __name__ == "__main__":
    import logging
    mqparser = MQParser()

    # determine logging level
    try:
        loglevel = getattr(logging, mqparser.args.loglevel.upper())
    except AttributeError:
        try:
            loglevel = int(mqparser.args.loglevel)
        except ValueError:
            loglevel = logging.DEBUG

    gui = MQGUI(mqparser.args, loglevel=loglevel)
