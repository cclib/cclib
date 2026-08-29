# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

from cclib.bridge import cclib2chemfiles
from cclib.parser.utils import find_package

from numpy.testing import assert_array_almost_equal

from ..test_data import getdatafile


class ChemfilesTest:
    """Tests for the chemfiles bridge in cclib."""

    def setup_method(self) -> None:
        if not find_package("chemfiles"):
            raise ImportError("Must install chemfiles to run this test")

        self.data, self.logfile = getdatafile("Gaussian", "basicGaussian16", ["dvb_sp.out"])

    def test_makechemfiles(self) -> None:
        """Check that geometry is transferred correctly from ccData to Frame"""
        from chemfiles.frame import Frame

        frame = cclib2chemfiles.makechemfiles(self.data)

        assert isinstance(frame, Frame)
        assert len(frame.atoms) == len(self.data.atomnos)

        assert_array_almost_equal(frame.positions, self.data.atomcoords[-1])

    def test_makechemfiles_masses(self) -> None:
        """Verify atomic masses are copied when present."""
        frame = cclib2chemfiles.makechemfiles(self.data)
        assert hasattr(self.data, "atommasses")
        frame_masses = [atom.mass for atom in frame.atoms]
        assert_array_almost_equal(frame_masses, self.data.atommasses)

    def test_makechemfiles_charges(self) -> None:
        """Verify selected charge scheme is transferred correctly."""
        charges = getattr(self.data, "atomcharges", None)

        assert isinstance(charges, dict)
        assert "mulliken" in charges

        frame = cclib2chemfiles.makechemfiles(self.data, "mulliken")
        frame_charges = [atom.charge for atom in frame.atoms]
        assert_array_almost_equal(frame_charges, charges["mulliken"])
