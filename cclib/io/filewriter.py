# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Generic file writer and related tools"""

try:
    from cclib.bridge import makeopenbabel
    import openbabel as ob
    import pybel as pb
    _has_openbabel = True
except ImportError:
    _has_openbabel = False

from math import sqrt
from collections import Iterable

from cclib.parser.utils import PeriodicTable


class MissingAttributeError(Exception):
    pass


class Writer(object):
    """Abstract class for writer objects."""

    required_attrs = ()

    def __init__(self, ccdata, jobfilename=None, indices=None, terse=False,
                 *args, **kwargs):
        """Initialize the Writer object.

        This should be called by a subclass in its own __init__ method.

        Inputs:
          ccdata - An instance of ccData, parsed from a logfile.
          jobfilename - The filename of the parsed logfile.
          indices - One or more indices for extracting specific geometries/etc. (zero-based)
          terse - Whether to print the terse version of the output file - currently limited to cjson/json formats
        """

        self.ccdata = ccdata
        self.jobfilename = jobfilename
        self.indices = indices
        self.terse = terse

        self.pt = PeriodicTable()

        self._check_required_attributes()

        # Open Babel isn't necessarily present.
        if _has_openbabel:
            # Generate the Open Babel/Pybel representation of the molecule.
            # Used for calculating SMILES/InChI, formula, MW, etc.
            self.obmol, self.pbmol = self._make_openbabel_from_ccdata()
            self.bond_connectivities = self._make_bond_connectivity_from_openbabel(self.obmol)

        self._fix_indices()

    def generate_repr(self):
        """Generate the written representation of the logfile data.

        This should be overriden by all the subclasses inheriting from
        Writer.
        """
        raise NotImplementedError(
            'generate_repr is not implemented for ' + str(type(self)))

    def _calculate_total_dipole_moment(self):
        """Calculate the total dipole moment."""
        return sqrt(sum(self.ccdata.moments[1] ** 2))

    def _check_required_attributes(self):
        """Check if required attributes are present in ccdata."""
        missing = [x for x in self.required_attrs
                   if not hasattr(self.ccdata, x)]
        if missing:
            missing = ' '.join(missing)
            raise MissingAttributeError(
                'Could not parse required attributes to write file: ' + missing)

    def _make_openbabel_from_ccdata(self):
        """Create Open Babel and Pybel molecules from ccData.
        """
        obmol = makeopenbabel(self.ccdata.atomcoords,
                              self.ccdata.atomnos,
                              charge=self.ccdata.charge,
                              mult=self.ccdata.mult)
        if self.jobfilename is not None:
            obmol.SetTitle(self.jobfilename)
        return (obmol, pb.Molecule(obmol))

    def _make_bond_connectivity_from_openbabel(self, obmol):
        """Based upon the Open Babel/Pybel molecule, create a list of tuples
        to represent bonding information, where the three integers are
        the index of the starting atom, the index of the ending atom,
        and the bond order.
        """
        bond_connectivities = []
        for obbond in ob.OBMolBondIter(obmol):
            bond_connectivities.append((obbond.GetBeginAtom().GetIndex(),
                                        obbond.GetEndAtom().GetIndex(),
                                        obbond.GetBondOrder()))
        return bond_connectivities

    def _fix_indices(self):
        """Clean up the index container type and remove zero-based indices to
        prevent duplicate structures and incorrect ordering when
        indices are later sorted.
        """
        if not self.indices:
            self.indices = set()
        elif not isinstance(self.indices, Iterable):
            self.indices = set([self.indices])
        # This is the most likely place to get the number of
        # geometries from.
        if hasattr(self.ccdata, 'atomcoords'):
            lencoords = len(self.ccdata.atomcoords)
            indices = set()
            for i in self.indices:
                if i < 0:
                    i += lencoords
                indices.add(i)
            self.indices = indices
        return
