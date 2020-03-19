from abc import abstractmethod, ABC, abstractstaticmethod, abstractproperty
from mspypeline.helpers import get_logger


class BaseReader(ABC):
    def __init__(self, start_dir: str, reader_config: dict, loglevel):
        self.full_data = {}
        self.start_dir = start_dir
        self.reader_config = reader_config
        self.logger = get_logger(self.get_class_name(), loglevel)
        if start_dir is None:
            raise ValueError("Invalid starting dir")
        if not reader_config:
            # raise ValueError("Empty configs")
            pass

    @classmethod
    def get_class_name(cls):
        return cls.__name__

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

    def get_full_data(self):
        if not self.full_data:
            raise NotImplementedError("full_data needs to be filled with the file information")
        return self.full_data


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
