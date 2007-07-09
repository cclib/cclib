__revision__ = "$Revision$"

import os
import unittest

from cclib.parser import ADF, GAMESS, GAMESSUK, Gaussian, Jaguar, Molpro


test_modules = [ "SP", "SPun", "GeoOpt", "Basis", "Core",
                 "MP", "CC", "CI", "TD",
                 "vib" ]


def getfile(parser, *location):
    """Returns a parsed logfile."""
    if parser.__name__ in ["ADF", "GAMESS", "Gaussian", "Jaguar", "Molpro"]:
        fullpath = ("..","data",parser.__name__) + location
    elif parser.__name__=="GAMESSUK":
        fullpath = ("..","data","GAMESS-UK") + location
    logfile = parser(os.path.join(*fullpath))
    logfile.logger.setLevel(0)
    logfile.parse()
    return logfile

def gettestdata(module=None):
    """Returns a dict of test files from a given module."""
    lines = open('testdata').readlines()
    lines = [line.split() for line in lines if line.strip()]
    if module:
        lines = [line for line in lines if line[0] == module]
    testdata = {}
    for line in lines:
        test = {}
        test["module"] = line[0]
        test["parser"] = line[1]
        test["location"] = line[3:]
        # The file is parsed bettertest.TestCase.run().
        #test["data"] = getfile(eval(line[1]), *line[3:])
        testclass = line[2]
        testdata[testclass] = test
    return testdata

def visualtests():
    """These are not formal tests -- but they should be eyeballed."""
    logfiles = [ getfile(Gaussian,"basicGaussian03","dvb_gopt.out"),
                 getfile(GAMESS,"basicPCGAMESS","dvb_gopt_a.out"),
                 getfile(GAMESS,"basicGAMESS-US","dvb_gopt_a.out"),
                 getfile(ADF,"basicADF2004.01","dvb_gopt.adfout"),
                 getfile(Jaguar,"basicJaguar4.2", "dvb_gopt.out"),
                 getfile(Jaguar,"basicJaguar6.5", "dvb_gopt.out"),
                 ]

    print "\n\nMO energies of optimised dvb"
    print "    ","".join(["%-10s" % x for x in ['Gaussian','PC-GAMESS','GAMESS-US','ADF','Jaguar4.2','Jaguar6.5']])
    print "HOMO", "   ".join(["%+7.4f" % x.moenergies[0][x.homos[0]] for x in logfiles])
    print "LUMO", "   ".join(["%+7.4f" % x.moenergies[0][x.homos[0]+1] for x in logfiles])
    print "H-L ", "   ".join(["%7.4f" % (x.moenergies[0][x.homos[0]+1]-x.moenergies[0][x.homos[0]],) for x in logfiles])

def importName(modulename, name):
    """Import from a module whose name is determined at run-time.

    Taken from Python Cookbook 2nd ed O'Reilly Recipe 16.3
    """
    try:
        module = __import__(modulename, globals(), locals(), [name])
    except ImportError:
        return None
    return getattr(module, name)

def testmodule(module):
    """Run all the unittests in a module."""

    testdata = gettestdata(module=module)    
    total = errors = failures = 0
    for test in testdata:

        print "\n**** test%s for %s ****" %(module, '/'.join(testdata[test]["location"]))
        test = importName("test%s" %module, test)
        parser = testdata[test.__name__]["parser"]
        location = testdata[test.__name__]["location"]
        test.data = getfile(eval(parser), *location)
        myunittest = unittest.makeSuite(test)

        a = unittest.TextTestRunner(verbosity=2).run(myunittest)
        total += a.testsRun
        errors += len(a.errors)
        failures += len(a.failures)

    print "\n\n********* SUMMARY OF %s TEST **************" %module.upper()
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)


def testall():
    """Run all unittests in all modules."""

    # Make sure we are in the test directory of this script,
    #   so that getfile() can access the data files.
    curdir = os.path.abspath(os.curdir)
    destdir = os.path.dirname(__file__)
    if destdir:
        os.chdir(destdir)

    perpackage = {}
    errors = []
    for module in test_modules:

        testdata = gettestdata(module)
        for name in testdata:

            path = '/'.join(testdata[name]["location"])
            program = testdata[name]["location"][0][5:]
            print "\n**** test%s for %s ****" %(module, path)

            try:
                test = importName("test%s" %module, name)
            except:
                error.append("ERROR: could not import %s from %s." %(name, module))
            else:
                parser = testdata[test.__name__]["parser"]
                location = testdata[test.__name__]["location"]
                test.data = getfile(eval(parser), *location)
                myunittest = unittest.makeSuite(test)
                a = unittest.TextTestRunner(verbosity=2).run(myunittest)
                l = perpackage.setdefault(program, [0, 0, 0])
                l[0] += a.testsRun
                l[1] += len(a.errors)
                l[2] += len(a.failures)

    print "\n\n********* SUMMARY PER PACKAGE ****************"
    names = perpackage.keys()
    names.sort()
    total = [0, 0, 0]
    print " "*14, "\t".join(["Total", "Passed", "Failed", "Errors"])
    for name in names:
        l = perpackage[name]
        print name.ljust(15), "%3d\t%3d\t%3d\t%3d" % (l[0], l[0]-l[1]-l[2], l[2], l[1])
        for i in range(3):
            total[i] += l[i]

    print "\n\n********* SUMMARY OF EVERYTHING **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total[0],total[0]-(total[1]+total[2]),total[2],total[1])

    if errors:
        print "\n".join(errors)

    print "\n\n*** Visual tests ***"
    visualtests()
    
    # Return to the directory we started from.
    if destdir:
        os.chdir(curdir)


if __name__ == "__main__":
    testall()
