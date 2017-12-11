# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for Molden writer."""

import os
import unittest

import cclib
from cclib.io.filewriter import MissingAttributeError
from cclib.io.moldenwriter import MoldenReformatter
from cclib.io.moldenwriter import round_molden

__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")
__testdir__ = __filedir__

class MOLDENTest(unittest.TestCase):

    def test_missing_attribute_error(self):
        """Check if MissingAttributeError is raised as expected."""
        fpath = os.path.join(__datadir__,
                             "data/GAMESS/basicGAMESS-US2014/dvb_un_sp.out")
        required_attrs = ['atomcoords', 'atomnos', 'natom']
        for attr in required_attrs:
            data = cclib.io.ccread(fpath)
            delattr(data, attr)

            # Molden files cannot be wriiten if required attrs are missing.
            with self.assertRaises(MissingAttributeError):
                cclib.io.moldenwriter.MOLDEN(data)

    def test_atoms_section_size(self):
        """Check if size of Atoms section is equal to expected."""
        fpath = os.path.join(__datadir__,
                             "data/GAMESS/basicGAMESS-US2014/dvb_un_sp.out")
        data = cclib.io.ccread(fpath)
        writer = cclib.io.moldenwriter.MOLDEN(data)
        # Check size of Atoms section.
        self.assertEqual(len(writer._coords_from_ccdata(-1)), data.natom)

    def test_gto_section_size(self):
        """Check if size of GTO section is equal to expected."""
        fpath = os.path.join(__datadir__,
                             "data/GAMESS/basicGAMESS-US2014/dvb_un_sp.out")
        data = cclib.io.ccread(fpath)
        writer = cclib.io.moldenwriter.MOLDEN(data)
        # Check size of GTO section.
        size_gto_ccdata = 0
        for atom in data.gbasis:
            size_gto_ccdata += 1
            for prims in atom:
                size_gto_ccdata += len(prims[1]) + 1
        # Filter blank lines.
        size_gto_writer = len(list(filter(None, writer._gto_from_ccdata())))
        self.assertEqual(size_gto_writer, size_gto_ccdata)

    def test_mo_section_size(self):
        """Check if size of MO section is equal to expected."""
        fpath = os.path.join(__datadir__,
                             "data/GAMESS/basicGAMESS-US2014/dvb_un_sp.out")
        data = cclib.io.ccread(fpath)
        writer = cclib.io.moldenwriter.MOLDEN(data)
        # Check size of MO section.
        size_mo_ccdata = 0
        extra = 4 if hasattr(data, 'mosyms') else 3
        for i in range(data.mult):
            size_mo_ccdata += len(data.moenergies[i]) *\
                                (len(data.mocoeffs[i][0]) + extra)
        # Filter blank lines.
        size_mo_writer = len(list(filter(None, writer._mo_from_ccdata())))
        self.assertEqual(size_mo_writer, size_mo_ccdata)

    def test_round_molden(self):
        """Check if Molden Style number rounding works as expected."""
        # If the 6th digit after dot is greater than 5, but is not 7,
        # round the number upto 6th place.
        # Else truncate at 6th digit after dot.
        self.assertEqual(round_molden(1), 1)
        self.assertEqual(round_molden(-1), -1)
        self.assertEqual(round_molden(0.999995789), 0.999995)
        self.assertEqual(round_molden(-0.999995789), -0.999995)
        self.assertEqual(round_molden(0.999996789), 0.999997)
        self.assertEqual(round_molden(-0.999997789), -0.999997)
        self.assertEqual(round_molden(0.999997789), 0.999997)
        self.assertEqual(round_molden(-0.999998789), -0.999999)
        self.assertEqual(round_molden(-0.999999999), -1.0)

    def test_molden_cclib_diff(self):
        """Check if file written by cclib matched file written by Molden."""
        filenames = ['dvb_un_sp', 'C_bigbasis', 'water_mp2']
        for fn in filenames:
            fpath = os.path.join(__datadir__,
                                 "data/GAMESS/basicGAMESS-US2014/"+fn+".out")
            data = cclib.io.ccread(fpath)
            cclib_out = cclib.io.moldenwriter.MOLDEN(data).generate_repr()
            # Reformat cclib's output to remove extra spaces.
            cclib_out_formatted = MoldenReformatter(cclib_out).reformat()
            fpath = os.path.join(__testdir__, "data/molden5.7_"+fn+".molden")
            with open(fpath) as handle:
                molden_out = handle.read()
            # Reformat Molden's output to remove extra spaces,
            # and fix number formatting.
            molden_out_formatted = MoldenReformatter(molden_out).reformat()
            # Assert if reformatted files from both writers are same.
            self.assertMultiLineEqual(molden_out_formatted,
                                      cclib_out_formatted)


if __name__ == "__main__":
    unittest.main()
