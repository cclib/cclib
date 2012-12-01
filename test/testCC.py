# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2007, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

__revision__ = "$Revision$"

import numpy

import bettertest


class GenericCCTest(bettertest.TestCase):
    """Coupled cluster unittest."""

    def testsign(self):
        corrections = self.data.ccenergies - self.data.scfenergies
        self.failUnless(numpy.alltrue(corrections < 0.0))


class GenericCCDTest(GenericCCTest):
    """CCD unittest."""
    
    def testsign(self):
        """CCD: Are the Coupled cluster corrections negative?"""
        super(GenericCCDTest, self).testsign()


class GenericCCSDTest(GenericCCTest):
    """CCSD unittest."""
    
    def testsign(self):
        """CCSD: Are the Coupled cluster corrections negative?"""
        super(GenericCCSDTest, self).testsign()


class GenericCCSDTTest(GenericCCTest):
    """CCSD(T) unittest."""
    
    def testsign(self):
        """CCSD(T): Are the Coupled cluster correction negative?"""
        super(GenericCCSDTTest, self).testsign()


class GAMESSUSCCDTest(GenericCCDTest):
    """GAMESS-US CCD unittest."""

    old_tests = ["GAMESS/GAMESS-US/water_ccd_2005.06.27.r3.out.gz"]


class GAMESSUSCCSDTest(GenericCCSDTest):
    """GAMESS-US CCSD unittest."""

    old_tests = ["GAMESS/GAMESS-US/water_ccsd_2005.06.27.r3.out.gz"]

    
class GAMESSUSCCSDTTest(GenericCCSDTTest):
    """GAMESS-US CCSD(T) unittest."""

    old_tests = ["GAMESS/GAMESS-US/water_ccsd(t)_2005.06.27.r3.out.gz"]

    
class GaussianCCDTest(GenericCCDTest):
    """Gaussian CCD unittest."""

    
class GaussianCCSDTest(GenericCCSDTest):
    """Gaussian CCSD unittest."""

    
class GaussianCCSDTTest(GenericCCSDTTest):
    """Gaussian CCSD(T) unittest."""

    
class MolproCCDTest(GenericCCDTest):
    """Molpro CCD unittest."""

    
class MolproCCSDTest(GenericCCSDTest):
    """Molpro CCSD unittest."""

    
class MolproCCSDTTest(GenericCCSDTTest):
    """Molpro CCSD(T) unittest."""


if __name__ == "__main__":

    from testall import testall
    testall(modules=["CC"])
