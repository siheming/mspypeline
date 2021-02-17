import logging
from logging import Logger
from abc import abstractmethod, ABC
from typing import Iterable
try:
    from ruamel_yaml import YAML
except ModuleNotFoundError:
    from ruamel.yaml import YAML

from mspypeline.helpers import get_logger, DataDict


class BaseReader(ABC):
    """
    | Base reader to provide a data dictionary with keys to the data. Data stored on system hardware, is thus only
      loaded on demand. This is the parent class of any file reader that will be use to preprocess data to the
      internal format.

    Example
    ---------
    >>> # example for a new custom reader
    >>> class CustomReader(BaseReader):
    ...     name = "reader"  # this is the name of the reader in the yaml file
    ...     required_files = []  # this is a list of strings of all files that should be parsed
    ...     plotter = BasePlotter
    ...
    ...     def __init__(self, start_dir, reader_config, loglevel):
    ...         super().__init__(start_dir, reader_config, loglevel)
    ...         for file in Reader.required_files:
    ...             self.full_data[file] = [0, 0, 10]  # this should be the data from the file
    >>> r = Reader("", {}, 10)

    """
    def __init__(self, start_dir: str, reader_config: dict, loglevel: int = logging.DEBUG):
        """

        Parameters
        ----------
        start_dir
            location where the directory/txt folder to the data can be found.
        reader_config
            mapping of the file reader configuration (as e.g. given in the config.yml file)
        loglevel
            level of the logger

        """
        self.full_data: DataDict = DataDict(data_source=self)
        self.start_dir: str = start_dir
        self.reader_config: dict = reader_config
        self.logger: Logger = get_logger(self.__class__.__name__, loglevel)

        # log which files will be read
        self.logger.info("Required files: %s", self.required_files)

        if not reader_config:
            self.logger.warning("Empty configs")
        else:
            self.logger.debug("Got configs: %s", self.reader_config)
        if start_dir is None:
            raise ValueError("Invalid starting dir")

    @property
    @classmethod
    @abstractmethod
    def name(cls) -> str:
        raise NotImplementedError

    @property
    @classmethod
    @abstractmethod
    def required_files(cls) -> Iterable[str]:
        raise NotImplementedError

    @property
    @classmethod
    @abstractmethod
    def plotter(cls):  # -> Type[BasePlotter]
        raise NotImplementedError


class MissingFilesException(Exception):
    pass


# the decorators add to the docstring so it is removed here
BaseReader.name.__doc__ = ""
BaseReader.required_files.__doc__ = ""
BaseReader.plotter.__doc__ = ""

if __name__ == "__main__":
    # minimal example of a class implementing the BaseReader
    class Reader(BaseReader):
        name = "reader"  # this is the name of the reader in the yaml file
        required_files = []  # this is a list of strings of all files that should be parsed
        #plotter = BasePlotter

        def __init__(self, start_dir, reader_config, loglevel):
            super().__init__(start_dir, reader_config, loglevel)
            for file in Reader.required_files:
                self.full_data[file] = [0, 0, 10]  # this should be the data from the file


    r = Reader("", {}, 10)
