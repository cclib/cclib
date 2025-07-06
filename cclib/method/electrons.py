# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Calculate properties for electrons."""

import logging
from typing import TYPE_CHECKING, Optional

from cclib.method.calculationmethod import Method

if TYPE_CHECKING:
    from cclib.parser.data import ccData
    from cclib.progress import Progress


class Electrons(Method):
    """A container for methods pertaining to electrons."""

    def __init__(
        self,
        data: "ccData",
        progress: Optional["Progress"] = None,
        loglevel: int = logging.INFO,
        logname: str = "Log",
    ):
        super().__init__(data, progress, loglevel, logname)
        self.required_attrs = ("atomnos", "charge", "coreelectrons", "homos")

    def __str__(self) -> str:
        """Returns a string representation of the object."""
        return "Electrons"

    def __repr__(self) -> str:
        """Returns a representation of the object."""
        return "Electrons"

    def alpha(self) -> int:
        """Number of alpha electrons"""
        return self.data.homos[0] + 1

    def beta(self) -> int:
        """Number of beta electrons"""
        return self.data.homos[-1] + 1

    def count(self, core: bool = False) -> int:
        """Returns the electron count in system.

        Normally returns electrons used in calculation, but will include
        core electrons in pseudopotentials if core is True.
        """
        nelectrons = sum(self.data.atomnos) - self.data.charge
        if core:
            nelectrons += sum(self.data.coreelectrons)
        return nelectrons
