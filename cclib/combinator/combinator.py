from dataclasses import dataclass

import cclib.parser_properties as cprops


@dataclass
class combinator:
    name: str
    job_list: list


@dataclass
class sp_combinator(combinator):
    def __init__(self):
        self.name = "single_point"
        self.job_list = [
            [
                cprops.scfenergies,
                cprops.atomcoords,
                cprops.atombasis,
                cprops.atomnos,
                cprops.charge,
                cprops.mult,
                cprops.nbasis,
                cprops.atommasses,
                cprops.mocoeffs,
                cprops.mosyms,
                cprops.nmo,
                cprops.natom,
                cprops.gbasis,  # after natom due to dependency inpsi4
            ]
        ]
