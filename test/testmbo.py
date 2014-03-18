# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2007-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

from __future__ import print_function
import os
import logging
import unittest

# If numpy is not installed, try to import Numeric instead.
try:
    import numpy
except ImportError:
    import Numeric as numpy

from testall import getfile
from cclib.method import MBO
from cclib.parser import Gaussian


def main(log=True):
    data, logfile = getfile(Gaussian, "basicGaussian09", "dvb_sp.out")
    mbo = MBO(data)
    if not log:
        mbo.logger.setLevel(logging.ERROR)
    mbo.calculate()

    return mbo


class MBOTest(unittest.TestCase):
    def test_mbo_sp(self):
        """Testing Mayer bond orders for restricted single point."""
        mbo = main(log=False)
        e_mbo = numpy.loadtxt("dvb_sp.mbo")
        self.assertTrue(numpy.all(mbo.fragresults >= e_mbo - 0.25))
        self.assertTrue(numpy.all(mbo.fragresults <= e_mbo + 0.25))


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(unittest.makeSuite(MBOTest))
