# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
import abc


class base_parser(abc.ABC):
    @staticmethod
    @abc.abstractmethod
    def parse(file_handler, program, ccdata):
        return

    @staticmethod
    def check_dependencies(dependency_list, ccdata, current_property):
        for i in dependency_list:
            prop = ccdata.__getattr__(i)
            if prop is None:
                return False
        return True
