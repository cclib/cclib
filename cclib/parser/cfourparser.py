# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for CFOUR output files"""

import datetime
import re

from cclib.parser import data, logfileparser, utils
from cclib.parser.logfileparser import StopParsing

import numpy as np


class CFOUR(logfileparser.Logfile):
    """A CFOUR v2.0 log file."""

    def __init__(self, *args, **kwargs):
        super().__init__(logname="CFOUR", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return f"CFOUR log file {self.filename}"

    def __repr__(self):
        """Return a representation of the object."""
        return f'CFOUR("{self.filename}")'

    def normalisesym(self, label):
        """Use standard symmetry labels instead of Gaussian labels.

        To normalise:
        (1) If label is one of [SG, PI, PHI, DE], replace by [sigma, pi, phi, delta]
        (2) replace any G or U by their lowercase equivalent
        """
        # note: DLT must come after DLTA
        greek = [("SG", "sigma"), ("PI", "pi"), ("PHI", "phi"), ("DE", "delta")]
        for k, v in greek:
            if label.startswith(k):
                tmp = label[len(k) :]
                label = v
                if tmp:
                    label = f"{v}.{tmp}"
        ans = label.replace("U", "u").replace("G", "g")
        return ans

    def before_parsing(self):
        pass

    def after_parsing(self):
        pass

    def extract(self, inputfile, line):
        if "Version" in line:
            self.metadata["pavkage_version"]=line.split()[1]
        if 'Coordinates used in calculation (QCOMP)' in line:
            line=next(inputfile)
            line=next(inputfile)
            line=next(inputfile)
            line=next(inputfile)
            line=next(inputfile)
            atomnos=[]
            while not '----------------------------------------------------------------' in line:
                atomnos.append(int(line.strip().split()[1]))
                line=next(inputfile)
            self.set_attribute("atomnos",np.array(atomnos))


