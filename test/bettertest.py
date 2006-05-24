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
        if not (first.shape==second.shape and
                first.typecode()==second.typecode() and
                Numeric.alltrue(first==second)):
            raise self.failureException, (msg or 'These two arrays are not identical')
