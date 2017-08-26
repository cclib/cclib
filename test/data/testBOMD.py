# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test Born-Oppenheimer molecular dynamics (BOMD) logfiles in cclib"""

import os
import unittest


__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericBOMDTest(unittest.TestCase):
    """Generic Born-Oppenheimer molecular dynamics unittest"""

    nsteps = 35
    natoms = 20

    def testdimscfenergies(self):
        """Are the number of parsed energies consistent with the number of MD
        steps?
        """
        self.assertEquals(self.data.scfenergies.shape, (self.nsteps, ))

    def testdimatomcoords(self):
        """Are the number of parsed geometries consistent with the number of
        MD steps?
        """
        self.assertEquals(self.data.atomcoords.shape, (self.nsteps, self.natoms, 3))

    def testdimtime(self):
        """Are the number of time points consistent with the number of MD
        steps?
        """
        print(self.data.time)
        self.assertEquals(self.data.time.shape, (self.nsteps, ))


if __name__=="__main__":

    import sys
    sys.path.append(os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['BOMD'])
    suite.testall()
