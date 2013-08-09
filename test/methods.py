# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2007, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

__revision__ = "$Rev$"

import sys
import unittest


testmodules = ['testpopulation', 'testcda']


def importname(modulename, name):
    """Import from a module whose name is determined at runtime.

    (Python Cookbook 2nd ed.)
    """
    module = __import__(modulename, globals(), locals(), [name])
    if not module:
        raise ImportError
    return getattr(module, name)


if __name__=="__main__":
    fullsuite = unittest.TestSuite()
    for testmodule in testmodules:
        try:
            tests = importname(testmodule, "tests")
            for test in tests:
                test = unittest.TestLoader().loadTestsFromTestCase(test)
                fullsuite.addTest(test)
        except ImportError:
            print("%s failed!" % testmodule)
    unittest.TextTestRunner(verbosity=2).run(fullsuite)
        
