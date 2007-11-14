__revision__ = "$Revision$"

import numpy

import bettertest


class GenericTDTest(bettertest.TestCase):
    """Time-dependent HF/DFT unittest."""

    def testenergies(self):
        """Is the l_max reasonable?"""
        self.assertEqual(len(self.data.etenergies), self.number)
        idx_lambdamax = [i for i, x in enumerate(self.data.etoscs)
                         if x==max(self.data.etoscs)][0]
        self.assertInside(self.data.etenergies[idx_lambdamax], 41000, 5000)
    
    def testoscs(self):
        """Is the maximum of eotscs in the right range?"""
        self.assertEqual(len(self.data.etoscs), self.number)
        self.assertInside(max(self.data.etoscs), 0.67, 0.1)

    def testsecs(self):
        """Is the sum of etsecs close to 1?"""
        self.assertEqual(len(self.data.etsecs), self.number)
        lowestEtrans = self.data.etsecs[1]
        sumofsec = sum([z*z for (x, y, z) in lowestEtrans])
        self.assertInside(sumofsec, 1.0, 0.16)
        
    def testsymsnumber(self):
        """Is the length of etsyms correct?"""
        self.assertEqual(len(self.data.etsyms), self.number)

class GaussianTDDFTTest(GenericTDTest):
    """Gaussian time-dependent HF/DFT unittest."""
    number = 5

    def testrotatsnumber(self):
        """Is the length of etrotats correct?"""
        self.assertEqual(len(self.data.etrotats), self.number)

class GAMESSUSTDDFTTest(GenericTDTest):
    """GAMESS time-dependent HF/DFT unittest."""
    number = 10

class PCGamessTDDFTTest(GenericTDTest):
    """PC-GAMESS time-dependent HF/DFT unittest."""
    number = 5

class OrcaTDDFTTest(GenericTDTest):
    """ORCA time-dependent HF/DFT unittest."""
    number = 24

if __name__=="__main__":

    from testall import testall
    testall(modules=["TD"])
