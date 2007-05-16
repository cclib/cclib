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
            print "%s failed!" % testmodule
    unittest.TextTestRunner(verbosity=2).run(fullsuite)
        
