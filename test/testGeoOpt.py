import os
import math
import Numeric
import unittest
import bettertest

from testall import getfile
from cclib.parser import ADF, GAMESS, Gaussian, Jaguar, GAMESSUK

class GenericGeoOptTest(bettertest.TestCase):
    def testhomos(self):
        """Is the index of the homo equal to 34?"""
        self.assertArrayEquals(self.data.homos,Numeric.array([34],"i"),"%s != array([34],'i')" % Numeric.array_repr(self.data.homos))

    def testatomcoords(self):
        """Are atomcoords consistent with natom and Angstroms?"""
        coords = self.data.atomcoords
        self.assertEquals(self.data.natom,len(coords[0]),"natom  is %d but len(atomcoords[0]) is %d" % (self.data.natom,len(coords[0])))

        # Find the minimum distance between two C atoms
        mindist = 999
        for i in range(self.data.natom-1):
            if self.data.atomnos[i]==6:
                for j in range(i+1,self.data.natom):
                    if self.data.atomnos[j]==6:
                        # Find the distance in the final iteration
                        dist = math.sqrt(sum((coords[-1][i]-coords[-1][j])**2))
                        mindist = min(mindist,dist)
        self.assert_(abs(mindist-1.34)<0.03,"Mindist is %f (not 1.34)" % mindist)

    def testatomnos(self):
        """Are the atomnos correct?"""
        self.assertEquals(self.data.atomnos.shape, (20,) )
        self.assertEquals(sum(self.data.atomnos==6) + sum(self.data.atomnos==1),
                          20)        

    def testatomcoords_more(self):
        """Are atomcoords consistent with geovalues?"""
        coords = self.data.atomcoords
        self.assertEquals(len(self.data.geovalues),len(coords),"len(atomcoords) is %d but len(geovalues) is %d" % (len(coords),len(self.data.geovalues)))
        
    def testnatom(self):
        """Is the number of atoms equal to 20?"""
        self.assertEquals(self.data.natom,20)

    def testnbasis(self):
        """Is the number of basis set functions equal to 60?"""
        self.assertEquals(self.data.nbasis,60)

    def testscfenergy(self):
        """Is the SCF energy within 40eV of -10365"""
        self.assertInside(self.data.scfenergies[-1],-10365,40,"Final scf energy: %f not -10365+-40eV" % self.data.scfenergies[-1])

    def testnormalisesym(self):
        """Did this subclass overwrite normalisesym?"""
        self.assertNotEquals(self.data.normalisesym("A"),"ERROR: This should be overwritten by this subclass")

    def testlengthmoenergies(self):
        """Is the number of evalues equal to 60?"""
        self.assertEquals(60,len(self.data.moenergies[0]))

    def testtypemoenergies(self):
        """Is moenergies a list containing one Numeric array?"""
        self.assertEquals(type(self.data.moenergies), type([]))
        self.assertEquals(type(self.data.moenergies[0]), type(Numeric.array([])))

    def testsymlabels(self):
        """Are all the symmetry labels either Ag/u or Bg/u?"""
        sumwronglabels = sum([x not in ['Ag','Bu','Au','Bg'] for x in self.data.mosyms[0]])
        self.assertEquals(sumwronglabels,0)

    def testscfvaluetype(self):
        """Do the scf values have the right type?"""
        self.assertEquals(type(self.data.scfvalues),type([]))
        self.assertEquals(type(self.data.scfvalues[0]),type(Numeric.array([])))

    def testgeotargets(self):
        """Do the geo targets have the right dimensions?"""
        self.assertEquals(self.data.geotargets.shape,(len(self.data.geovalues[0]),))
    
    def testscfvaluedim(self):
        """Do the scf values have the right dimensions?"""
        self.assertEquals(len(self.data.scfvalues),len(self.data.geovalues))

    def testscftargetdim(self):
        """Do the scf targets have the right dimensions?"""
        self.assertEquals(self.data.scftargets.shape,(len(self.data.scfvalues),len(self.data.scfvalues[0][0])))

class GaussianGeoOptTest(GenericGeoOptTest):
    def setUp(self):
        self.data = data[0]

class PCGamessGeoOptTest(GenericGeoOptTest):
    def setUp(self):
        self.data = data[1]

class GamessUSGeoOptTest(GenericGeoOptTest):
    def setUp(self):
        self.data = data[2]

class ADFGeoOptTest(GenericGeoOptTest):
    def setUp(self):
        self.data = data[3]

    def testscfvaluedim(self):
        """Do the scf values have the right dimensions?"""
        # ADF calculations have one more SCF cycle after the geometry is converged
        # with a tighter SCF convergence cutoff
        self.assertEquals(len(self.data.scfvalues),len(self.data.geovalues)+1)

    def testatomcoords_more(self):
        """Are atomcoords consistent with geovalues?"""
        coords = self.data.atomcoords
        self.assertEquals(len(self.data.geovalues),len(coords)-1,"len(atomcoords)-1 is %d but len(geovalues) is %d" % (len(coords)-1,len(self.data.geovalues)))
    
    def testscfenergy(self):
        """Is the SCF energy within 1eV of -140eV"""
        self.assertInside(self.data.scfenergies[-1],-140,1,"Final scf energy: %f not -140+-1eV" % self.data.scfenergies[-1])

class Jaguar42GeoOptTest(GenericGeoOptTest):
    def setUp(self):
        self.data = data[4]

class Jaguar65GeoOptTest(GenericGeoOptTest):
    def setUp(self):
        self.data = data[5]
        
class GamessUKGeoOptTest(GenericGeoOptTest):
    def setUp(self):
        self.data = data[6]
        
names = [ "Gaussian", "PCGamess", "GAMESS", "ADF", "Jaguar 4.2",
          "Jaguar 6.5", "GAMESS UK" ]
tests = [ GaussianGeoOptTest, PCGamessGeoOptTest,
          GamessUSGeoOptTest, ADFGeoOptTest,
          Jaguar42GeoOptTest, Jaguar65GeoOptTest,
          GamessUKGeoOptTest ]
data = [ getfile(Gaussian,"basicGaussian03","dvb_gopt.out"),
         getfile(GAMESS,"basicPCGAMESS","dvb_gopt_b.out"),
         getfile(GAMESS,"basicGAMESS-US","dvb_gopt_a.out"),
         getfile(ADF,"basicADF2004.01","dvb_gopt_b.adfout"),
         getfile(Jaguar,"basicJaguar4.2","dvb_gopt_b.out"),
         getfile(Jaguar,"basicJaguar6.5","dvb_gopt_b.out"),
         getfile(GAMESSUK,"basicGAMESS-UK","dvb_gopt_d.out")]         

if __name__=="__main__":
    total = errors = failures = 0
    for name,test in zip(names,tests):
        print "\n**** Testing %s Geo Opt ****" % name
        myunittest = unittest.makeSuite(test)
        a = unittest.TextTestRunner(verbosity=2).run(myunittest)
        total += a.testsRun
        errors += len(a.errors)
        failures += len(a.failures)
    
    print "\n\n********* SUMMARY OF Geo Opt **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)
