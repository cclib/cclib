from abc import ABC
import abc

class base_parser(ABC):

    @staticmethod
    @abc.abstractmethod
    def parse(file_handler, program, ccdata):
        return

