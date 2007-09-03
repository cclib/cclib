__revision__ = "$Revision$"

import numpy

import bettertest


class GenericCITest(bettertest.TestCase):
    """CI unittest."""
    
    def testnumberofstates(self):
        """ Are there nstate elements in etenergies/etsecs/etsyms/etoscs?"""
        self.assertEqual(len(self.data.etenergies), self.nstates)
        self.assertEqual(len(self.data.etsecs), self.nstates)
        self.assertEqual(len(self.data.etsyms), self.nstates)
        self.assertEqual(len(self.data.etoscs), self.nstates)
    
    def testetenergies(self):
        """ Are transition energies positive and rising?"""
        self.failUnless(numpy.alltrue(self.data.etenergies > 0.0))
        changes = self.data.etenergies[1:] - self.data.etenergies[:-1]
        self.failUnless(numpy.alltrue(changes > 0.0))

class GenericCISTest(GenericCITest):
    """CIS unittest."""
    
class GenericCISWaterTest(GenericCISTest):
    """CIS(RHF)/STO-3G water unittest."""

    # First four singlet/triplet state excitation energies [cm-1].
    # Based on output in GAMESS test.
    etenergies0 = numpy.array([98614.56, 114906.59, 127948.12, 146480.64])
    etenergies1 = numpy.array([82085.34,  98999.11, 104077.89, 113978.37])

    # First four singlet/triplet state excitation orbitals and coefficients.
    # Tuples: (from MO, to MO, coefficient) - don't need spin indices.
    # Based on output in GAMESS test (using the "dets" algorithm).
    # Note that only coefficients larger than 0.1 are included here, as
    #   the Gaussian test does not contain smaller ones.
    # The simple example of water should yield the same first 4 states in all programs.
    # Note that the GAMESS test "water_cis_dets" calculates also triplet states,
    #  but the resulting transition dipoles and oscillator strengths are not correct.
    etsecs0 = [ [(4, 5, -0.70710678)],
                [(4, 6, -0.70710678)],
                [(3, 5,  0.68368723)],
                [(2, 5,  0.31163855), (3, 6, -0.63471970)] ]
    etsecs1 = [ [(4, 5,  0.70710678)],
                [(2, 6, -0.16333506), (3, 6, -0.68422741)],
                [(4, 6,  0.70710678)],
                [(2, 5, -0.37899667), (3, 6, -0.59602876)] ]
                
    def testetenergiesvalues(self):
        """ Are etenergies within 50cm-1 of the correct values?"""
        indices0 = [i for i in range(self.nstates) if self.data.etsyms[i][0] == "S"]
        indices1 = [i for i in range(self.nstates) if self.data.etsyms[i][0] == "T"]
        singlets = [self.data.etenergies[i] for i in indices0]
        triplets = [self.data.etenergies[i] for i in indices1]
        # All programs do singlets.
        singletdiff = singlets[:4] - self.etenergies0
        self.failUnless(numpy.alltrue(singletdiff < 50))
        # Not all programs do triplets (i.e. Jaguar).
        if len(triplets) >= 4:
            tripletdiff = triplets[:4] - self.etenergies1
            self.failUnless(numpy.alltrue(tripletdiff < 50))

    def testetsecsvalues(self):
        """ Are etsecs correct and coefficients within 0.0005 of the correct values?"""
        indices0 = [i for i in range(self.nstates) if self.data.etsyms[i][0] == "S"]
        indices1 = [i for i in range(self.nstates) if self.data.etsyms[i][0] == "T"]
        singlets = [self.data.etsecs[i] for i in indices0]
        triplets = [self.data.etsecs[i] for i in indices1]
        # All programs do singlets.
        found = False
        for i in range(4):
            for exc in self.etsecs0[i]:
                for s in singlets[i]:
                    if s[0][0] == exc[0] and s[1][0] == exc[1]:
                        found = True
                        self.assertInside(abs(s[2]), abs(exc[2]), 0.0005)
        if not found:
            self.fail("Excitation %i->%s not found (singlet state %i)" %(exc[0], exc[1], i))
        # Not all programs do triplets (i.e. Jaguar).
        if len(triplets) >= 4:
            found = False
            for i in range(4):
                for exc in self.etsecs1[i]:
                    for s in triplets[i]:
                        if s[0][0] == exc[0] and s[1][0] == exc[1]:
                            found = True
                            self.assertInside(abs(s[2]), abs(exc[2]), 0.0005)
            if not found:
                self.fail("Excitation %i->%s not found (triplet state %i)" %(exc[0], exc[1], i))

class GAMESSUSCISTest(GenericCISWaterTest):
    """GAMESS-US CIS(RHF)/STO-3G water unittest."""

    nstates = 5
        
    def testnocoeffs(self):
        """(MP2) Are Natural Orbital coefficients the right size?"""
        self.assertEquals(self.data.nocoeffs.shape, (self.data.nmo, self.data.nbasis))

class GaussianCISTest(GenericCISWaterTest):
    """Gaussian CIS(RHF)/STO-3G water unittest."""

    nstates = 10

    def testnocoeffs(self):
        """(MP2) Are Natural Orbital coefficients the right size?"""
        self.assertEquals(self.data.nocoeffs.shape, (self.data.nmo, self.data.nbasis))

class Jaguar65CISTest(GenericCISWaterTest):
    """Jaguar 6.5 CIS(RHF)/STO-3G water unittest."""

    nstates = 5

             
if __name__=="__main__":

    from testall import testall
    testall(modules=["CI"])
