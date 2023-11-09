from dataclasses import dataclass
from cclib.parser_properties import scfenergies
from cclib.parser_properties import atomcoords
from cclib.parser_properties import atomnos
from cclib.parser_properties import atommasses
from cclib.parser_properties import nbasis

@dataclass
class combinator():
    name: str
    job_list: list


@dataclass
class sp_combinator(combinator):
    def __init__(self):
        self.name = 'single_point'
        self.job_list =  [[scfenergies,atomcoords,atomnos,atommasses,nbasis]]

