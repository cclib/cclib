# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2016, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Test BOMD logfiles in cclib"""

import os
import unittest

# import numpy


__filedir__ = os.path.realpath(os.path.dirname(__file__))

class GenericBOMDTest(unittest.TestCase):
    """Generic Born-Oppenheimer molecular dynamics unittest"""

    nsteps = 35
    natoms = 20

    def testdimscfenergies(self):
        """"""
        self.assertEquals(self.data.scfenergies.shape, (self.nsteps, ))

    def testdimatomcoords(self):
        """"""
        self.assertEquals(self.data.atomcoords.shape, (self.nsteps, self.natoms, 3))

    def testdimtime(self):
        """"""
        print(self.data.time)
        self.assertEquals(self.data.time.shape, (self.nsteps, ))


if __name__=="__main__":

    import sys
    sys.path.append(os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['BOMD'])
    suite.testall()
