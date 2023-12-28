from cclib.parser_properties import utils
from cclib.parser_properties.base_parser import base_parser

import numpy as np


class atombasis(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> list | None:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        dependency_list = ['nmo','nbasis']
        line = file_handler.last_line
        constructed_data = None
        if (
            line[5:35] == "Molecular Orbital Coefficients"
            or line[5:41] == "Alpha Molecular Orbital Coefficients"
            or line[5:40] == "Beta Molecular Orbital Coefficients"
        ):
                constructed_data=[]
                if not base_parser.check_dependencies(dependency_list, ccdata, 'atombasis'):
                    return None
                beta=False
                if line[5:40] == "Beta Molecular Orbital Coefficients":
                    beta = True
                colNames = file_handler.virtual_next()
                symmetries = file_handler.virtual_next()
                eigenvalues = file_handler.virtual_next()
                base = 0 
                atombasis=[]
                for base in range(0, ccdata.nmo, 5):
                    for i in range(ccdata.nbasis):
                        line = file_handler.virtual_next()
                        if i == 0:
                            # Find location of the start of the basis function name
                            start_of_basis_fn_name = line.find(line.split()[3]) - 1
                        if base == 0 and not beta:  # Just do this the first time 'round
                            parts = line[:start_of_basis_fn_name].split()
                            if len(parts) > 1:  # New atom
                                if i > 0:
                                    constructed_data.append(atombasis)
                                atombasis = []
                            atombasis.append(i)

                        atombasis.append(i)
                return constructed_data
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> list | None:
        constructed_data = None
        if program in atombasis.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(atombasis, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
