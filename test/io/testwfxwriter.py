# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for writer WFXwriter module."""

import os
import unittest

import cclib

from cclib.io.filewriter import MissingAttributeError

__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")


class WFXTest(unittest.TestCase):

    def test_missing_attribute_error(self):
        """Check if MissingAttributeError is raised as expected."""
        fpath = os.path.join(__datadir__,
                             "data/GAMESS/basicGAMESS-US2014/C_bigbasis.out")
        required_attrs = ('atomcoords', 'atomnos', 'gbasis', 'charge',
                          'homos', 'mult')
        for attr in required_attrs:
            data = cclib.io.ccopen(fpath).parse()
            delattr(data, attr)

            # WFX files cannot be written if required attrs are missing.
            with self.assertRaises(MissingAttributeError):
                cclib.io.wfxwriter.WFXWriter(data)

    def test_no_of_prims(self):
        """Check if number of primitives are calculated correctly."""
        num_orb = {'s':1, 'p':3, 'd':6, 'f':10, 'g':15, 'h':21}
        gamessdir = os.path.join(__datadir__,
                                 "data/GAMESS/basicGAMESS-US2014")
        filenames = ["C_bigbasis.out", "dvb_un_sp.out"]
        filepaths = [os.path.join(gamessdir, fn) for fn in filenames]

        for fpath in filepaths:
            data = cclib.io.ccopen(fpath).parse()
            wfx = cclib.io.wfxwriter.WFXWriter(data)

            no_prims_writer = wfx._no_of_prims()
            no_prims_ccdata = 0
            for atom in data.gbasis:
                for prims in atom:
                    no_prims_ccdata += num_orb[prims[0].lower()]\
                                        * len(prims[1])

            self.assertEqual(no_prims_writer, no_prims_ccdata)

    def test_no_of_electrons(self):
        """Check if number of electrons are being calculated correctly."""
        gamessdir = os.path.join(__datadir__,
                                 "data/GAMESS/basicGAMESS-US2014")
        filenames = ["C_bigbasis.out", "dvb_un_sp.out"]
        filepaths = [os.path.join(gamessdir, fn) for fn in filenames]

        for fpath in filepaths:
            data = cclib.io.ccopen(fpath).parse()
            wfx = cclib.io.wfxwriter.WFXWriter(data)

            nelectrons = wfx._no_electrons()
            nalpha = wfx._no_alpha_electrons()
            nbeta = wfx._no_beta_electrons()
            maxelectrons = data.homos[0]*2 + 2
            # Sum of alpha and beta electrons be the number of electrons.
            self.assertEqual(nelectrons, nalpha + nbeta)
            # Number of electrons must be either homos*2, or (homos*2 - 1).
            self.assertTrue(maxelectrons - 1 <= nelectrons <= maxelectrons)

    def test_mo_normalization(self):
        """Check if MO section is printed correctly."""
        fpath = os.path.join(__datadir__,
                             "data/GAMESS/basicGAMESS-US2014/C_bigbasis.out")
        data = cclib.io.ccopen(fpath).parse()
        wfx = cclib.io.wfxwriter.WFXWriter(data)

        normalized_mocoeffs = wfx._normalized_mocoeffs()
        if len(data.homos) > 1:
            self.assertEqual(len(normalized_mocoeffs), wfx._no_electrons())
        else:
            self.assertEqual(len(normalized_mocoeffs), wfx._no_of_mos())
        self.assertEqual(len(normalized_mocoeffs[0]), wfx._no_of_prims())


if __name__ == "__main__":
    unittest.main()
