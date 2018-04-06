# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the CDA method in cclib"""

from __future__ import print_function

import sys
import os
import logging
import unittest

import numpy

sys.path.insert(1, "..")

from ..test_data import getdatafile
from cclib.method import CDA
from cclib.parser import Gaussian


def main(log=True):
    data1, logfile1 = getdatafile(Gaussian, "CDA", ["BH3CO-sp.log"])
    data2, logfile2 = getdatafile(Gaussian, "CDA", ["BH3.log"])
    data3, logfile3 = getdatafile(Gaussian, "CDA", ["CO.log"])
    fa = CDA(data1)
    if not log:
        fa.logger.setLevel(logging.ERROR)
    fa.calculate([data2, data3])

    return fa


def printResults():
    fa = main()

    print("       d       b       r")
    print("---------------------------")

    spin = 0
    for i in range(len(fa.donations[0])):

        print("%2i: %7.3f %7.3f %7.3f" % (i,
                                            fa.donations[spin][i],
                                            fa.bdonations[spin][i],
                                            fa.repulsions[spin][i]))


    print("---------------------------")
    print("T:  %7.3f %7.3f %7.3f" % (fa.donations[0].sum(),
                                        fa.bdonations[0].sum(),
                                        fa.repulsions[0].sum()))
    print("\n\n")


class CDATest(unittest.TestCase):

    def runTest(self):
        """Testing CDA results against Frenking's code"""
        fa = main(log=False)

        donation = fa.donations[0].sum()
        bdonation = fa.bdonations[0].sum()
        repulsion = fa.repulsions[0].sum()

        self.assertAlmostEqual(donation, 0.181, 3)
        self.assertAlmostEqual(bdonation, 0.471, 3)
        self.assertAlmostEqual(repulsion, -0.334, 3)


if __name__ == "__main__":
    printResults()
    unittest.TextTestRunner(verbosity=2).run(unittest.makeSuite(CDATest))
