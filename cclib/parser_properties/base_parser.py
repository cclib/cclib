import abc
from abc import ABC


class base_parser(ABC):
    @staticmethod
    @abc.abstractmethod
    def parse(file_handler, program, ccdata):
        return

    @staticmethod
    def check_dependencies(dependency_list, ccdata, current_property):
        for i in dependency_list:
            if getattr(ccdata, i) == None:
                return False
        return True
