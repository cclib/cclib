#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for general file writer."""

import inspect
import os
import unittest

import cclib


__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")


class FileWriterTest(unittest.TestCase):

    def test_init(self):
        """Does the class initialize properly?"""

        # You cannot instantiate a class with abstract methods.
        fpath = os.path.join(__datadir__, "data/ADF/basicADF2007.01/dvb_gopt.adfout")
        data = cclib.io.ccread(fpath)
        with self.assertRaises(TypeError):
            cclib.io.filewriter.Writer(data)


if __name__ == "__main__":
    unittest.main()
