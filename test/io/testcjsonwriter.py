# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for the CJSON writer."""

import json
import os
import unittest
from math import sqrt

import cclib

__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")


class CJSONTest(unittest.TestCase):
    """Unit tests for the CJSON writer."""

    def test_init(self):
        """Does the class initialize correctly?"""
        fpath = os.path.join(__datadir__, "data/ADF/basicADF2007.01/dvb_gopt.adfout")
        data = cclib.io.ccread(fpath)
        cjson = cclib.io.cjsonwriter.CJSON(data)

        # The object should keep the ccData instance passed to its constructor.
        self.assertEqual(cjson.ccdata, data)

    def test_cjson_generation(self):
        """Does the CJSON format get generated properly?"""
        fpath = os.path.join(__datadir__, "data/ADF/basicADF2007.01/NH3.adfout")
        data = cclib.io.ccread(fpath)

        cjson = cclib.io.cjsonwriter.CJSON(data).generate_repr()

        # The data available in the cjson and ccdata objects should be equal.
        json_data = json.loads(cjson)
        number_of_atoms = json_data['properties']['number of atoms']
        self.assertEqual(number_of_atoms, data.natom)

        dipole_moment = json_data['properties']['total dipole moment']
        self.assertAlmostEqual(
            dipole_moment,
            sqrt(sum(data.moments[1] ** 2))
        )

        # Ensure the bond connectivity index starts from 0
        bonds = json_data.get('bonds', None)
        self.assertIsNotNone(bonds)
        indices = bonds['connections']['index']
        self.assertEqual(min(indices), 0)
        self.assertTrue(max(indices) < number_of_atoms)

    def test_zero_dipole_moment(self):
        """Does the CJSON writer handle zero dipole moment correctly?"""
        fpath = os.path.join(__datadir__, "data/GAMESS/basicGAMESS-US2017/C_bigbasis.out")
        data = cclib.io.ccopen(fpath).parse()

        cjson = cclib.io.cjsonwriter.CJSON(data).generate_repr()

        json_data = json.loads(cjson)
        self.assertAlmostEqual(json_data["properties"]['total dipole moment'], 0.0)

    def test_missing_dipole_moment(self):
        """Does the CJSON writer handle missing properties correctly?"""
        fpath = os.path.join(__datadir__, "data/GAMESS/basicGAMESS-US2017/C_bigbasis.out")
        data = cclib.io.ccopen(fpath).parse()
        del data.moments

        cjson = cclib.io.cjsonwriter.CJSON(data).generate_repr()

        json_data = json.loads(cjson)
        self.assertFalse("total dipole moment" in json_data["properties"])


if __name__ == "__main__":
    unittest.main()
