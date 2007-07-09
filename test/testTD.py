__revision__ = "$Revision$"

# If numpy is not installed, try to import Numeric instead.
try:
    import numpy
except ImportError:
    import Numeric as numpy

import bettertest


class GenericTDTest(bettertest.TestCase):
    """Time-dependent HF/DFT unittest."""

    def testenergiesnumber(self):
        """Is the length of etenergies correct?"""
        self.assertEqual(len(self.data.etenergies), self.number)
    
    def testoscsnumber(self):
        """Is the length of eotscs correct?"""
        self.assertEqual(len(self.data.etoscs), self.number)

    def testrotatsnumber(self):
        """Is the length of etrotats correct?"""
        self.assertEqual(len(self.data.etrotats), self.number)

    def testsecsnumber(self):
        """Is the length of etsecs correct?"""
        self.assertEqual(len(self.data.etsecs), self.number)

    def testsymsnumber(self):
        """Is the length of etsyms correct?"""
        self.assertEqual(len(self.data.etsyms), self.number)

class GaussianTDDFTTest(GenericTDTest):
    """Gaussian time-dependent HF/DFT unittest."""

    number = 24
    
    def testsyms(self):
        """Is etsyms populated by singlets and triplets 50/50?"""
        singlets = [sym for sym in self.data.etsyms if "Singlet" in sym]
        triplets = [sym for sym in self.data.etsyms if "Triplet" in sym]
        self.assertEqual(len(singlets), self.number/2)
        self.assertEqual(len(triplets), self.number/2)


if __name__=="__main__":

    from testall import testmodule
    testmodule("TD")
