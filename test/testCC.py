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
    """Coupled-Cluster unittest."""

    def testsign(self):
        corrections = self.data.ccenergies - self.data.scfenergies
        self.failUnless(numpy.alltrue(corrections < 0.0))


class GenericCCDTest(GenericCCTest):
    """CCD unittest."""
    
    def testsign(self):
        """CCD: Are the Coupled-Cluster correction negative?"""
        super(GenericCCDTest, self).testsign()


class GenericCCSDTest(GenericCCTest):
    """CCSD unittest."""
    
    def testsign(self):
        """CCSD: Are the Coupled-Cluster correction negative?"""
        super(GenericCCSDest, self).testsign()


class GenericCCSDTest(GenericCCTest):
    """CCSD(T) unittest."""
    
    def testsign(self):
        """CCSD(T): Are the Coupled-Cluster correction negative?"""
        super(GenericCCSDTest, self).testsign()


class GAMESSUSCCDTest(GenericCCDTest):
    """GAMESS-US CCD unittest."""


class GAMESSUSCCSDTest(GenericCCSDTest):
    """GAMESS-US CCSD unittest."""

    
class GAMESSUSCCSDTest(GenericCCSDTest):
    """GAMESS-US CCSD(T) unittest."""

    
class GaussianCCDTest(GenericCCDTest):
    """Gaussian CCD unittest."""

    
class GaussianCCSDTest(GenericCCSDTest):
    """Gaussian CCSD unittest."""

    
class GaussianCCSDTest(GenericCCSDTest):
    """Gaussian CCSD(T) unittest."""

    
class MolproCCDTest(GenericCCDTest):
    """Molpro CCD unittest."""

    
class MolproCCSDTest(GenericCCSDTest):
    """Molpro CCSD unittest."""

    
class MolproCCSDTest(GenericCCSDTest):
    """Molpro CCSD(T) unittest."""


if __name__ == "__main__":

    from testall import testall
    testall(modules=["CC"])
