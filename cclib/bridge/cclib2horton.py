# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in horton (http://theochem.github.io/horton)."""

import numpy
from cclib.parser.data import ccData
from cclib.parser.utils import find_package

_found_iodata = find_package("iodata")

# Detect whether iodata (part of horton 3) is present or not
# Horton 3 is divided into smaller (sub)packages that each take different functionalities.
if _found_iodata:
    from iodata import IOData
    from iodata.orbitals import MolecularOrbitals


def check_horton():
    if not _found_iodata:
        raise ImportError("You must install `horton` to use this function.")


def makehorton(data):
    """ Create horton IOData object from ccData object """

    check_horton()
    attributes = {}

    # In horton 3 IOData, inputs are type-checked -- thus the bridge also verifies the types.
    # if coordinates are known; numpy.ndarray (size natom-by-3)
    if hasattr(data, "atomcoords"):
        attributes["atcoords"] = numpy.asanyarray(data.atomcoords)[-1]
    # if atomic numbers are known; numpy.ndarray (size natom)
    if hasattr(data, "atomnos"):
        attributes["atnums"] = numpy.asanyarray(data.atomnos)
    # if orbital coefficients known; iodata.orbitals.MolecularOrbitals
    if hasattr(data, "mocoeffs"):
        # First build a dictionary of inputs that will be passed to the constructor of
        # horton3's MolecularOrbitals object.
        moattr = {
            "kind": "restricted",
            "norba": data.nbasis,
            "norbb": None,
            # In horton 3, occupation in MOs are represented as 1's (occupied) and
            # 0's (unoccupied). Beta orbitals follow directly after alpha orbitals, forming
            # 1D array.
            "occs": numpy.concatenate(
                numpy.ones(data.homos[0]), numpy.zeros(data.nbasis - data.homos[0])
            ),
            "coeffs": data.mocoeffs[0],
            "energies": None,
            "irreps": None,
        }
        # and if unrestricted:
        if len(mocoeffs) == 2:
            moattr["kind"] = "unrestricted"
            moattr["norbb"] = data.nbasis
            moattr["coeffs"].append(data.mocoeffs[1].T)
            moattr["occs"].append(
                numpy.concatenate(
                    numpy.ones(data.homos[1]),
                    numpy.zeros(data.nbasis - data.homos[1]),
                )
            )
        # Then construct MolecularOrbitals object
        attributes["mo"] = MolecularOrbitals(**moattr)
    # if multiplicity known; float / should not be set when mocoeffs present
    # Refer to IOData code:
    # https://github.com/theochem/iodata/blob/b36513d162f99b57264005583701c6987037839c/iodata/iodata.py#L174
    if (
        hasattr(data, "mult")
        and not hasattr(data, "mocoeffs")
    ):
        attributes["spinpol"] = data.mult - 1  # horton has 2S+1, iodata has 2S
    # if pseudopotentials exist; numpy.ndarray (size natom)
    if hasattr(data, "coreelectrons"):
        attributes["atcorenums"] = data.atomnos - numpy.asanyarray(data.coreelectrons)
    # if mulliken charges are known; dict of numpy.ndarrays (size natom)
    if hasattr(data, "atomcharges"):
        attributes["atcharges"] = data.atomcharges

    return IOData(**attributes)  # Pass collected attributes into IOData constructor


def makecclib(iodat):
    """ Create cclib ccData object from horton IOData object """

    check_horton()
    attributes = {}

    # Horton 3 IOData class uses attr and does not have __dict__.
    # In horton 3, some attributes have a default value of None.
    # Therefore, second hasattr statement is needed for mo attribute.
    if hasattr(iodat, "atcoords"):
        # cclib parses the whole history of coordinates in the list, horton keeps the last one.
        attributes["atomcoords"] = [iodat.atcoords]
    if hasattr(iodat, "mo") and hasattr(iodat.mo, "norba"):
        # MO coefficient should be transposed to match the dimensions.
        attributes["mocoeffs"] = [iodat.mo.coeffs[: iodat.mo.norba].T]
        if iodat.mo.kind == "unrestricted":
            attributes["mocoeffs"].append(iodat.mo.coeffs[iodat.mo.norba :].T)
    if hasattr(iodat, "spinpol") and isinstance(iodat.spinpol, int):
        # IOData stores 2S, ccData stores 2S+1.
        attributes["mult"] = iodat.spinpol + 1
    if hasattr(iodat, "atnums"):
        attributes["atnums"] = numpy.asanyarray(iodat.atnums)
    if hasattr(iodat, "atcorenums") and isinstance(iodat.atnums, numpy.ndarray):
        # cclib stores num of electrons screened out by pseudopotential
        # horton stores num of electrons after applying pseudopotential
        attributes["coreelectrons"] = numpy.asanyarray(iodat.atnums) - numpy.asanyarray(
            iodat.atcorenums
        )
    if hasattr(iodat, "atcharges"):
        attributes["atomcharges"] = iodat.atcharges

    return ccData(attributes)
