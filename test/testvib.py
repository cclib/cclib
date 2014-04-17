# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

import bettertest


class GenericIRTest(bettertest.TestCase):
    """Generic vibrational frequency unittest."""

    # Unit tests should normally give this value for the largest IR intensity.
    max_IR_intensity = 100

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
        self.assertInside(max(self.data.vibirs), self.max_IR_intensity, 10)


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

    # We have not been able to determine why ORCA gets such a different
    # maximum IR intensity. The coordinates are exactly the same, and
    # the basis set seems close enough to other programs. It would be nice
    # to determine whether this difference is algorithmic in nature,
    # but in the meanwhile we will expect to parse this value.
    max_IR_intensity = 215

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

    # The tolerance for this number has been increased, since ORCA
    # failed to make it inside +/-5, but it would be nice in the future
    # to determine is it's not too much work whether this is due to
    # algorithmic differences, or to differences in the input basis set
    # or coordinates. The first would be OK, but in the second case the
    # unit test jobs should be made more comparable. With cclib, we first
    # of all want to succeed in parsing, but would also like to remain
    # as comparable between programs as possible (for these tests).
    # Note also that this value is adjusted for Gaussian - why?
    def testramanintens(self):
        """Is the maximum Raman intensity 575 +/- 8 A**4/amu?"""
        self.assertInside(max(self.data.vibramans), 575, 8)

        # We used to test this, but it seems to vary wildly between
        # programs... perhaps we could use it if we knew why...
        #self.assertInside(self.data.vibramans[1], 2.6872, 0.0001)

    def testvibdisps(self):
        """Is the length and value of vibdisps correct?"""
        assert hasattr(self.data, "vibdisps")
        assert len(self.data.vibdisps) == 54

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
