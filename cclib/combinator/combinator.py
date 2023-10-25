# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from dataclasses import dataclass

from cclib.parser_properties import scfenergies


@dataclass
class combinator:
    name: str
    job_list: list


@dataclass
class sp_combinator(combinator):
    def __init__(self):
        self.name = "single_point"
        self.job_list = [[scfenergies]]
