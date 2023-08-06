# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Analyses related to orbitals."""

import logging
from typing import Optional

import numpy

from cclib.method.calculationmethod import Method
from cclib.progress import Progress


class Orbitals(Method):
    """A class for orbital related methods."""

    def __init__(self, data, progress: Optional[Progress] = None, \
                 loglevel: int = logging.INFO, logname: str = "Log") -> None:
        super().__init__(data, progress, loglevel, logname)

        self.required_attrs = ('mocoeffs','moenergies','homos')
        self.fragresults = None

    def __str__(self) -> str:
        """Return a string representation of the object."""
        return "Orbitals"

    def __repr__(self) -> str:
        """Return a representation of the object."""
        return "Orbitals"

    def closed_shell(self, precision: float = 10e-6) -> bool:
        """Return Boolean indicating if system is closed shell."""

        # If there are beta orbitals, we can assume the system is closed
        # shell if the orbital energies are identical within numerical accuracy.
        if len(self.data.mocoeffs) == 2:
            return numpy.allclose(*self.data.moenergies, atol=precision)

        # Restricted open shell will have one set of MOs but two HOMO indices,
        # and the indices should be different (otherwise it's still closed shell).
        if len(self.data.homos) == 2 and self.data.homos[0] != self.data.homos[1]:
            return False

        return True
