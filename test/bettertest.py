__revision__ = "$Revision$"

import unittest

import numpy

from cclib.parser import ADF, GAMESS, GAMESSUK, Gaussian, Jaguar, Molpro


class TestCase(unittest.TestCase):
    """Create a class with extra 'asserts' for testing numerical data,
        and a special run() method for loading cclib test files on run-time.

    It is not possible to test equality of numpy arrays using assertEquals().
    Instead, use assertArrayEquals() as defined below. For the original solution see:
    http://mail.python.org/pipermail/python-list/2005-November/311235.html

    Also, for testing near equality of floats use assertInside.
    (Taken from Python Cookbook 2nd Ed. Recipe 8.11)
    """

    def assertInside(self, first, second, error, msg=None):
        """Fail if the second number isn't within a certain error of the first."""
        if not (second-error) < first < (second+error):
            raise self.failureException, (msg or '%r != %r (+-%r)' % (first,second,error))

    def assertArrayEquals(self, first, second, msg=None):
        """Fails unless two numpy arrays are identical."""
        
        errormsg = None
        if not first.shape == second.shape:
            errormsg = "Shapes are different: %s != %s" % (first.shape, second.shape)
        # If NumPy was not imported, assume Numeric was.
        try:
            type1 = first.dtype
            type2 = second.dtype
        except:
            type1 = first.typecode()
            type2 = second.typecode()
        if not type1 == type2:
            errormsg = "Array types are different: %s != %s" % (type1, type2)
        if not numpy.alltrue(first == second):
            errormsg = "Not equal: %s != %s" % (first, second)
        if errormsg:
            raise self.failureException, (msg or errormsg)

    def run(self, result=None):
        """Custom run method for cclib."""
        
        # Skip the actuall call to run() if we want to skip the test,
        #   for any reason (feature not implemented, etc.).
        # Still, increment the total, and append to a new "skipped" attribute.
        # In Python 2.4, the document string is in _TestCase__testMethodDoc.
        # In Python 2.5, the document string is in _testMethodDoc.
        if hasattr(self, "_TestCase__testMethodDoc"):
            if "PASS" in self._TestCase__testMethodDoc:
                print self._TestCase__testMethodDoc
                result.testsRun += 1
                if not hasattr(result, "skipped"):
                    result.skipped = []
                result.skipped.append(self._testMethodDoc)
                return
        if hasattr(self, "_testMethodDoc"):
            if "PASS" in self._testMethodDoc:
                print self._testMethodDoc
                result.testsRun += 1
                if not hasattr(result, "skipped"):
                    result.skipped = []
                result.skipped.append(self._testMethodDoc)
                return

        # If test not skipped, run the test now.
        unittest.TestCase.run(self, result=result)
