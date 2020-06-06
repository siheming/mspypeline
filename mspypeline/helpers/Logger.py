import logging
import warnings
from typing import Union


def get_logger(name: str = None, loglevel: Union[int, str] = logging.DEBUG) -> logging.Logger:
    """

    Parameters
    ----------
    name
        name of the logger, if none the filename will be used instead
    loglevel
        loglevel of the logger, either an int or the str level names of the logging module

    Returns
    -------
    logging.Logger
        A logger with respective level and formatted output

    """
    if name is None:
        name = __name__
    logger = logging.getLogger(name)

    # determine logging level
    try:
        loglevel = getattr(logging, loglevel.upper())
    except AttributeError:
        try:
            loglevel = int(loglevel)
        except ValueError:
            warnings.warn(f"Did not recognize loglevel={loglevel}, setting loglevel to DEBUG", RuntimeWarning)
            loglevel = logging.DEBUG

    if not logger.handlers:
        ch = logging.StreamHandler()
        ch.setLevel(loglevel)
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        ch.setFormatter(formatter)

        logger.setLevel(loglevel)
        logger.addHandler(ch)
    return logger
