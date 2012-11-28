# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

__revision__ = "$Revision$"

import bettertest


class GenericIRTest(bettertest.TestCase):
    """Generic vibrational frequency unittest."""

    def testvibdisps(self):
        """Are the dimensions of vibdisps consistent with 3N-6 x N x 3"""
        numvib = 3*len(self.data.atomnos) - 6
        self.assertEqual(self.data.vibdisps.shape,
                        (numvib, len(self.data.atomnos), 3))

    def testlengths(self):
        """Are the lengths of vibfreqs and vibirs correct?"""
        numvib = 3*len(self.data.atomnos) - 6
        self.assertEqual(len(self.data.vibfreqs), numvib)
        self.assertEqual(len(self.data.vibirs), numvib)

    def testfreqval(self):
        """Is the highest freq value 3630 +/- 200 cm-1?"""
        self.assertInside(self.data.vibfreqs[-1], 3630, 200)

    def testirintens(self):
        """Is the maximum IR intensity 100 +/- 10 km mol-1?"""
        self.assertInside(max(self.data.vibirs), 100, 10)


class GenericIRimgTest(bettertest.TestCase):
    """Generic imaginary vibrational frequency unittest."""

    def testvibdisps(self):
        """Are the dimensions of vibdisps consistent with 3N-6 x N x 3"""
        numvib = 3*len(self.data.atomnos) - 6
        self.assertEqual(self.data.vibdisps.shape,
                        (numvib, len(self.data.atomnos), 3))

    def testlengths(self):
        """Are the lengths of vibfreqs and vibirs correct?"""
        numvib = 3*len(self.data.atomnos) - 6
        self.assertEqual(len(self.data.vibfreqs), numvib)
        self.assertEqual(len(self.data.vibirs), numvib)

    def testfreqval(self):
        """Is the lowest freq value negative?"""
        self.assertTrue(self.data.vibfreqs[0] < 0)


##    def testmaxvibdisps(self):
##        """What is the maximum value of displacement for a H vs a C?"""
##        Cvibdisps = compress(self.data.atomnos==6, self.data.vibdisps, 1)
##        Hvibdisps = compress(self.data.atomnos==1, self.data.vibdisps, 1)
##        self.assertEqual(max(abs(Cvibdisps).flat), 1.0)
        

class ADFIRTest(GenericIRTest):
    """ADF vibrational frequency unittest."""

    
class GamessUKIRTest(GenericIRTest):
    """GAMeSS-UK vibrational frequency unittest."""


class GamessUSIRTest(GenericIRTest):
    """GAMESS-US vibrational frequency unittest."""


class GaussianIRTest(GenericIRTest):
    """Gaussian vibrational frequency unittest."""

    def testvibsyms(self):
        """Is the length of vibsyms correct?"""
        numvib = 3*len(self.data.atomnos) - 6        
        self.assertEqual(len(self.data.vibsyms), numvib)

       
class JaguarIRTest(GenericIRTest):
    """Jaguar vibrational frequency unittest."""

    def testvibsyms(self):
            """Is the length of vibsyms correct?"""
            numvib = 3*len(self.data.atomnos) - 6        
            self.assertEqual(len(self.data.vibsyms), numvib)


class MolproIRTest(GenericIRTest):
    """Molpro vibrational frequency unittest."""


class OrcaIRTest(GenericIRTest):
    """ORCA vibrational frequency unittest."""


class PCGamessIRTest(GenericIRTest):
    """PC-GAMESS vibrational frequency unittest."""

    def testirintens(self):
        """Is the maximum IR intensity 135 +/- 5 km mol-1?"""
        self.assertInside(max(self.data.vibirs), 135, 5)     


class GenericRamanTest(bettertest.TestCase):
    """Generic Raman unittest."""

    def testlengths(self):
        """Is the length of vibramans correct?"""
        numvib = 3*len(self.data.atomnos) - 6
        self.assertEqual(len(self.data.vibramans), numvib)

    def testramanintens(self):
        """Is the maximum Raman intensity 575 +/- 5 A**4/amu?"""
        self.assertInside(max(self.data.vibramans), 575, 5)


class GamessUKRamanTest(GenericRamanTest):
    """GAMESS-UK Raman unittest."""


class GaussianRamanTest(GenericRamanTest):
    """Gaussian Raman unittest."""

    def testramanintens(self):
        """Is the maximum Raman intensity 1066 +/- 5 A**4/amu?"""
        self.assertInside(max(self.data.vibramans), 1066, 5)


class MolproRamanTest(GenericRamanTest):
    """Molpro Raman unittest."""


class OrcaRamanTest(GenericRamanTest):
    """ORCA Raman unittest."""

    
class PCGamessRamanTest(GenericRamanTest):
    """PC-GAMESS Raman unittest."""


class GamessUSIRimgTest(GenericIRimgTest):
    """GAMESS-US imaginary vibrational frequency unittest."""


if __name__=="__main__":

    from testall import testall
    testall(modules=["vib"])
