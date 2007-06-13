__revision__ = "$Revision$"

import os
import unittest

# If numpy is not installed, try to import Numeric instead.
try:
    import numpy
except ImportError:
    import Numeric as numpy

from cclib.parser.utils import PeriodicTable

import bettertest
from testall import gettestdata


class GenericCoreTest(bettertest.TestCase):
    """Core electrons"""
    coredict = {'Mo': 28, 'O':0, 'Cl':10}
    def testcorrect(self):
        """Is coreelectrons equal to what it should be?"""
        pt = PeriodicTable()
        ans = []
        for x in self.data.atomnos:
            ans.append(self.coredict[pt.element[x]])
        ans = numpy.array(ans, "i")
        self.assertArrayEquals(self.data.coreelectrons, ans)

class ADFCoreTest(GenericCoreTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]
        self.coredict = {'Mo': 28, 'O':0, 'Cl':0}

class GAMESSUKCoreTest(GenericCoreTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GAMESSUSCoreTest(GenericCoreTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GaussianCoreTest(GenericCoreTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]


 # Load test data using information in file.
testdata = gettestdata(module="Core")
           
if __name__=="__main__":
    total = errors = failures = 0
    for test in testdata:
        module = testdata[test]["module"]
        print "\n**** test%s for %s ****" %(module, '/'.join(testdata[test]["location"]))
        test = eval(test)
        myunittest = unittest.makeSuite(test)
        a = unittest.TextTestRunner(verbosity=2).run(myunittest)
        total += a.testsRun
        errors += len(a.errors)
        failures += len(a.failures)

    print "\n\n********* SUMMARY OF CORE TEST **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)
