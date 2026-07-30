# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
"""Bridge for using cclib data in Chemfiles (https://github.com/chemfiles/chemfiles.py)."""

from collections import defaultdict
from typing import Optional

from cclib.parser.data import ccData
from cclib.parser.utils import PeriodicTable, find_package


class MissingAttributeError(Exception):
    pass


_found_chemfiles = find_package("chemfiles")
if _found_chemfiles:
    from chemfiles import Atom, Frame


def _check_chemfiles(found_chemfiles: bool) -> None:
    if not found_chemfiles:
        raise ImportError("You must install `chemfiles` to use this function")


def makechemfiles(data: ccData, charge: Optional[str] = None) -> "Frame":
    """Create a Chemfiles frame from a ccData object.

    Inputs:
        data - ccData object
        charge - str or None, optional.
            Key in `data.atomcharges` selecting a partial-charge scheme (e.g. "mulliken").
            If None, no charges are assigned.
    """
    _check_chemfiles(_found_chemfiles)

    required_attrs = {"atomcoords", "atomnos"}
    missing = [x for x in required_attrs if not hasattr(data, x)]

    if missing:
        missing = " ".join(missing)
        raise MissingAttributeError(
            f"Could not create chemfiles frame due to missing attribute: {missing}"
        )

    masses = getattr(data, "atommasses", None)
    charge_arr = None

    if charge is not None:
        charges = getattr(data, "atomcharges", None)
        if charges is None:
            raise MissingAttributeError(
                "Charge method requested but ccData has no atomcharges attribute"
            )
        if charge not in charges:
            raise ValueError(
                f"Charge type '{charge}' not available. Available: {', '.join(charges.keys())}"
            )
        charge_arr = charges[charge]

    frame = Frame()
    pt = PeriodicTable()
    counts = defaultdict(int)
    for i, (atomno, coords) in enumerate(zip(data.atomnos, data.atomcoords[-1])):
        symbol = pt.element[atomno]
        counts[symbol] += 1
        # Chemfiles atom names should be unique within a frame; generate
        # element-based sequential names (e.g., C1, C2, H1, H2).
        name = f"{symbol}{counts[symbol]}"

        atom = Atom(name, type=symbol)

        if masses is not None:
            atom.mass = masses[i]

        if charge_arr is not None:
            atom.charge = charge_arr[i]

        frame.add_atom(atom, coords)

    return frame
