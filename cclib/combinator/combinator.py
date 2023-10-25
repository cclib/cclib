from dataclasses import dataclass
from cclib.parser_properties import scfenergies

@dataclass
class combinator():
    name: str
    job_list: list


@dataclass
class sp_combinator(combinator):
    def __init__(self):
        self.name = 'single_point'
        self.job_list =  [[scfenergies]]

