from dataclasses import dataclass
from primitive_parsers import sp_parser

@dataclass
class combinator():
    name: str = 'base'
    job_list: list = []


@dataclass
class sp_combinator(combinator):
    name: str = 'single_point'
    job_list: list = [[sp_energies]]

