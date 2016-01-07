# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2007,2008,2012,2014,2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Test gradients of Herzberg-Teller (HT) transition dipole moments of 
electronic transitions obtained from a numerical TD/CIS frequency calculation. 
"""

import os
import unittest

import numpy


__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericHTTest(unittest.TestCase):
    """Generic TD/CIS num freq HT unittest"""

    number = 3

    def testeteltrdipgradsshape(self):
        """Is eteltrdips the right shape?"""
        self.assertEqual(numpy.shape(self.data.eteltrdipgrads),
                         (self.number, self.data.natom, 3, 3))


class GaussianHTTest(GenericHTTest):
    """Customized TD/CIS num freq HT unittest"""

    def testetveleltrdipgradsshape(self):
        """Is etveleltrdipgrads the right shape?"""
        self.assertEqual(numpy.shape(self.data.etveleltrdipgrads),
                         (self.number, self.data.natom, 3, 3))

    def testetmagtrdipgradsshape(self):
        """Is etmagtrdipgrads the right shape?"""
        self.assertEqual(numpy.shape(self.data.etmagtrdipgrads),
                         (self.number, self.data.natom, 3, 3))

if __name__ == "__main__":

    import sys
    sys.path.append(os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['HT'])
    suite.testall()
