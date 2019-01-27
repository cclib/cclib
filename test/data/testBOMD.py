# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test Born-Oppenheimer molecular dynamics (BOMD) logfiles in cclib."""

import os
import unittest


__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericBOMDTest(unittest.TestCase):
    """Generic Born-Oppenheimer molecular dynamics unittest"""

    # To calculate the initial set of forces/velocities, programs
    # first converge the energy at the input geometry, so there is
    # always one more geometry/energy than MD step.
    nsteps = 35
    nenergies = 36

    def testdimscfenergies(self):
        """Are the number of parsed energies consistent with the number of MD
        steps?
        """
        self.assertEqual(self.data.scfenergies.shape, (self.nenergies, ))

    def testdimatomcoords(self):
        """Are the number of parsed geometries consistent with the number of
        MD steps?
        """
        self.assertEqual(self.data.atomcoords.shape, (self.nenergies, 20, 3))

    def testdimtime(self):
        """Are the number of time points consistent with the number of MD
        steps?
        """
        self.assertEqual(self.data.time.shape, (self.nsteps, ))


class GaussianBOMDTest(GenericBOMDTest):
    """Customized Born-Oppenheimer molecular dynamics unittest"""

    # This may have something to do with our corrections for
    # extrapolation step rejection.
    nenergies = 35


if __name__=="__main__":

    import sys
    sys.path.append(os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['BOMD'])
    suite.testall()
