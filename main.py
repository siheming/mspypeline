import logging
from mqpipeline import MQParser, MQGUI

if __name__ == "__main__":
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