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

    def testsecs_transition(self):
        """Is the lowest E transition from the HOMO or to the LUMO?"""
        idx_minenergy = [i for i, x in enumerate(self.data.etenergies)
                         if x==min(self.data.etenergies)][0]
        sec = self.data.etsecs[idx_minenergy]
        t = [(c*c, s, e) for (s, e, c) in sec]
        t.sort()
        t.reverse()        
        self.assert_(t[0][1][0]==self.data.homos[0] or
                    t[0][2][0]==self.data.homos[0]+1, t[0])
        
    def testsymsnumber(self):
        """Is the length of etsyms correct?"""
        self.assertEqual(len(self.data.etsyms), self.number)


class ADFTDDFTTest(GenericTDTest):
    """ADF time-dependent DFT unittest."""
    number = 5

    def testsecs(self):
        """Is the sum of etsecs close to 1?"""
        self.assertEqual(len(self.data.etsecs), self.number)
        lowestEtrans = self.data.etsecs[1]

        #ADF squares the etsecs
        sumofsec = sum([z for (x, y, z) in lowestEtrans])
        self.assertInside(sumofsec, 1.0, 0.16)


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
    number = 10
    def testoscs(self):
        """Is the maximum of eotscs in the right range?"""
        self.assertEqual(len(self.data.etoscs), self.number)
        self.assertInside(max(self.data.etoscs), 1.1, 0.1)


class GenericTDTesttrp(GenericTDTest):
    """Time-dependent HF/DFT (triplet) unittest."""

    def testenergies(self):
        """Is the l_max reasonable?"""
        self.assertEqual(len(self.data.etenergies), self.number)
        idx_lambdamax = [i for i, x in enumerate(self.data.etoscs)
                         if x==max(self.data.etoscs)][0]
        self.assertInside(self.data.etenergies[idx_lambdamax], 24500, 100)
    
    def testoscs(self):
        """Triplet excitations should be disallowed."""
        self.assertEqual(len(self.data.etoscs), self.number)
        self.assertInside(max(self.data.etoscs), 0.0, 0.01)


class GAMESSUSTDDFTtrpTest(GenericTDTesttrp):
    """GAMESS TD DFT (restricted) triplet unittest."""
    number = 5
    def testsymsnumber(self):
        """Is the length of etsyms correct? PASS"""
        pass


class PCGamessTDDFTtrpTest(GenericTDTesttrp):
    """PC-GAMESS TD DFT (restricted) triplet unittest."""
    number = 5


if __name__=="__main__":

    from testall import testall
    testall(modules=["TD"])
