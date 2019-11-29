import logging
import datetime

loglevel = logging.DEBUG


class Logger:
    def __init__(self, name=None, loglevel=loglevel):
        if name is None:
            name = __name__

        self.logger = logging.getLogger(name)

        if not self.logger.handlers:
            ch = logging.StreamHandler()
            ch.setLevel(loglevel)
            formatter = logging.Formatter(
                "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            )
            ch.setFormatter(formatter)

            self.logger.setLevel(loglevel)
            self.logger.addHandler(ch)

