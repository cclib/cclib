# -*- coding: utf-8 -*-
#
# Copyright (c) 2016, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for writer cjsonwriter module."""

import os
import unittest
import json

import cclib


__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")


class CJSONTest(unittest.TestCase):

    def setUp(self):
        self.CJSON = cclib.io.CJSONWriter

    def test_init(self):
        """Does the class initialize correctly?"""
        fpath = os.path.join(__datadir__, "data/ADF/basicADF2007.01/dvb_gopt.adfout")
        data = cclib.io.ccopen(fpath).parse()
        cjson = cclib.io.cjsonwriter.CJSON(data)

        # The object should keep the ccData instance passed to its constructor.
        self.assertEqual(cjson.ccdata, data)

    def test_cjson_generation(self):
        """Does the CJSON format get generated properly?"""
        fpath = os.path.join(__datadir__, "data/ADF/basicADF2007.01/dvb_gopt.adfout")
        data = cclib.io.ccopen(fpath).parse()

        cjson = cclib.io.cjsonwriter.CJSON(data).generate_repr()

        # The data available in the cjson and ccdata objects should be equal.
        json_data = json.loads(cjson)
        number_of_atoms = json_data['properties']['number of atoms']
        self.assertEqual(number_of_atoms, data.natom)

if __name__ == "__main__":
    unittest.main()
