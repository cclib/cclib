import abc
from abc import ABC


class base_parser(ABC):
    @staticmethod
    @abc.abstractmethod
    def parse(file_handler, program, ccdata):
        return
