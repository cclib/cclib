"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 1 $"

import Numeric
import os
from testall import getfile
from cclib.parser import Gaussian
from cclib.method import CDA

parser1 = getfile(Gaussian, os.path.join("CDA","BH3CO-sp.log"))
parser2 = getfile(Gaussian, os.path.join("CDA","BH3.log"))
parser3 = getfile(Gaussian, os.path.join("CDA","CO.log"))

fa = CDA(parser1)
fa.calculate([parser2, parser3])

def printResults():
    print "       d       b       r"
    print "---------------------------"

    spin = 0
    for i in range(len(fa.donations[0])):

        print "%2i: %7.3f %7.3f %7.3f"%(i,fa.donations[spin][i], fa.bdonations[spin][i], \
                                        fa.repulsions[spin][i])
            

    print "---------------------------"
    print "T:  %7.3f %7.3f %7.3f"%(reduce(Numeric.add, fa.donations[0]), \
                reduce(Numeric.add, fa.bdonations[0]), reduce(Numeric.add, fa.repulsions[0]))

donation = reduce(Numeric.add, fa.donations[0])
bdonation = reduce(Numeric.add, fa.bdonations[0])
repulsion = reduce(Numeric.add, fa.repulsions[0])

print donation, bdonation, repulsion

assert "%.3f" % donation == "0.181"
assert "%.3f" % bdonation == "0.471"
assert "%.3f" % repulsion == "-0.334"


