import Numeric
import unittest

class TestCase(unittest.TestCase):
    """Create a class with extra 'asserts' for testing numerical data.

    It is not possible to test equality of Numeric arrays using assertEquals().
    Instead use assertArrayEquals() defined below.
    (For the original solution see:
    http://mail.python.org/pipermail/python-list/2005-November/311235.html)

    Also, for testing near equality of floats use assertInside.
    (Taken from Python Cookbook 2nd Ed. Recipe 8.11)
    """
    def assertInside(self,first,second,error,msg=None):
        """Fail if the second number isn't within a certain error of the first."""
        if not (second-error) < first < (second+error):
            raise self.failureException, (msg or '%r != %r (+-%r)' % (first,second,error))

    def assertArrayEquals(self,first,second,msg=None):
        """Fails unless two Numeric arrays are identical."""
        errormsg = None
        if not first.shape==second.shape:
            errormsg = "Shapes are different: %s != %s" % (first.shape, second.shape)
        if not first.typecode()==second.typecode():
            errormsg = "Typecodes are differnts: %s != %s" % (first.typecode(), second.typecode())
        if not Numeric.alltrue(first==second):
            errormsg = "Not equal: %s != %s" % (first, second)
        if errormsg:
            raise self.failureException, (msg or errormsg)
