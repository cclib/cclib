from dataclasses import dataclass
import cclib.parser_properties as cpp

@dataclass
class combinator():
    name: str
    job_list: list


@dataclass
class sp_combinator(combinator):
    def __init__(self):
        self.name = 'single_point'
        self.job_list =  [[cpp.scfenergies,cpp.atomcoords,cpp.atomnos,cpp.atommasses,cpp.nbasis,cpp.nmo,cpp.charge,cpp.mult]]

