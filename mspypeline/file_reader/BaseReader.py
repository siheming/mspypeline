import logging
from abc import abstractmethod, ABC, abstractstaticmethod, abstractproperty
from mspypeline.helpers import get_logger


class DataDict(dict):
    def __init__(self, data_source, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.data_source = data_source

    def __missing__(self, key):
        try:
            self.data_source.logger.debug("Reading %s from disk", key)
            data = getattr(self.data_source, f"preprocess_{key}")()
            self[key] = data
            return data
        except FileNotFoundError as e:
            raise KeyError("Missing file:", key, e)
        except AttributeError as e:
            raise KeyError("Missing function to load:", key, e)


class BaseReader(ABC):
    def __init__(self, start_dir: str, reader_config: dict, loglevel=logging.DEBUG):
        self.full_data = DataDict(data_source=self)
        self.start_dir = start_dir
        self.reader_config = reader_config
        self.logger = get_logger(self.__class__.__name__, loglevel)

        # log which files will be read
        self.logger.info("Reading files: %s", self.all_files)
        self.logger.debug("Got configs: %s", self.reader_config)
        if start_dir is None:
            raise ValueError("Invalid starting dir")
        if not reader_config:
            self.logger.warning("Empty configs")

    @property
    @classmethod
    @abstractmethod
    def name(cls):
        raise NotImplementedError

    @property
    @classmethod
    @abstractmethod
    def all_files(cls):
        raise NotImplementedError

    @property
    @classmethod
    @abstractmethod
    def plotter(cls):
        raise NotImplementedError


class MissingFilesException(Exception):
    pass


if __name__ == "__main__":
    # minimal example of a class implementing the BaseReader
    class Reader(BaseReader):
        name = "reader"  # this is the name of the reader in the yaml file
        all_files = []  # this is a list of strings of all files that should be parsed

        def __init__(self, start_dir, reader_config, loglevel):
            super().__init__(start_dir, reader_config, loglevel)
            for file in Reader.all_files:
                self.full_data[file] = [0, 0, 10]  # this should be the data from the file


    r = Reader("", {}, 10)
