# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from dataclasses import dataclass

from primitive_parsers import sp_parser


@dataclass
class combinator:
    name: str = "base"
    job_list: list = []


@dataclass
class sp_combinator(combinator):
    name: str = "single_point"
    job_list: list = [[sp_energies]]
