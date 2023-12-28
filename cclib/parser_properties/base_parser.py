import abc
from abc import ABC


class base_parser(ABC):
    @staticmethod
    @abc.abstractmethod
    def parse(file_handler, program, ccdata):
        return

    @staticmethod
    def check_dependencies(dependency_list, ccdata, current_property):
        if len(dependency_list) > 0:
            for i in dependency_list:
                if not hasattr(ccdata, i):
                    raise Warning(
                        f"to parse {current_property}, ccdata required dependency of {i} which is not present"
                    )
                    return False
        return True
