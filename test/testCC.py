__revision__ = "$Revision$"

# If numpy is not installed, try to import Numeric instead.
try:
    import numpy
except ImportError:
    import Numeric as numpy

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


if __name__=="__main__":

    from testall import testmodule
    testmodule("CC")
