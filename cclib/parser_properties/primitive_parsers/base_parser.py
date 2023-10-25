from abc import ABC


class base_parser(ABC):
    @staticmethod
    @abstractmethod
    def parse(file_handler, program, ccdata):
        return
