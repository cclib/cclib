__revision__ = "$Revision$"

import numpy

from cclib.parser.utils import PeriodicTable
import bettertest


class GenericCoreTest(bettertest.TestCase):
    """Core electrons unittest."""

    def testcorrect(self):
        """Is coreelectrons equal to what it should be?"""
        pt = PeriodicTable()
        ans = []
        for x in self.data.atomnos:
            ans.append(self.coredict[pt.element[x]])
        ans = numpy.array(ans, "i")
        self.assertArrayEquals(self.data.coreelectrons, ans)

class ADFCoreTest(GenericCoreTest):
    """ADF core electrons unittest."""

    coredict = {'Mo': 28, 'O':0, 'Cl':0}

class GAMESSUKCoreTest(GenericCoreTest):
    """GAMESS-UK core electrons unittest."""

    coredict = {'Mo': 28, 'O':0, 'Cl':10}

class GAMESSUSCoreTest(GenericCoreTest):
    """GAMESS-US core electrons unittest."""

    coredict = {'Mo': 28, 'O':0, 'Cl':10}

class GaussianCoreTest(GenericCoreTest):
    """Gaussian core electrons unittest."""

    coredict = {'Mo': 28, 'O':0, 'Cl':10}

           
if __name__=="__main__":

    from testall import testall
    testall(modules=["Core"])
