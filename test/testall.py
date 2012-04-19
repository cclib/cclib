# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

__revision__ = "$Revision$"

import os
import sys
import unittest

from cclib.parser import ADF, GAMESS, GAMESSUK, Gaussian, Jaguar, Molpro, ORCA


# The modules to be included in the global test testall().
test_modules = [ "SP", "SPun", "GeoOpt", "Basis", "Core",   # Basic calculations.
                 "MP", "CC", "CI", "TD", "TDun",            # Post-SCF calculations.
                 "vib" ]                                    # Other property calculations.


def getfile(parser, *location):
    """Returns a parsed logfile.
    
    Inputs:
        parser - a logfile parser class (subclass of LogFile)
        *location - subdirectory and data filename(s)
    
    Outputs:
        data - the resulting data object
        logfile - the parser object used for parsing
    """

    # GAMESS-UK files are in the GAMESS path.
    if parser.__name__=="GAMESSUK":
        location = ("..", "data", "GAMESS-UK") + location
    else:
        location = ("..", "data", parser.__name__) + location

    # Now construct the proper full path(s).
    # Multiple paths will be in a list only if more than one data file given.
    # Presently, location contains only one subdirectory (basic*),
    #   so this is easy since there are normally 5 elements in location.
    if len(location) == 5:
        filename = os.path.join(*location)
    else:
        filename = [os.path.join(*(location[:4]+location[n:n+1])) for n in range(4,len(location))]

    logfile = parser(filename)
    logfile.logger.setLevel(0)
    data = logfile.parse()
    
    return data, logfile

def gettestdata(module=None):
    """Returns a dict of test files for a given module."""

    lines = open('testdata').readlines()

    # Remove blank lines and those starting with '#'.
    lines = [line.split() for line in lines if (line.strip() and line[0] != '#')]
    
    # Filter for lines only for the given module.
    if module:
        lines = [line for line in lines if line[0] == module]

    testdata = {}
    for line in lines:
        test = {}
        test["module"] = line[0]
        test["parser"] = line[1]
        test["location"] = line[3:]
        testclass = line[2]
        testdata[testclass] = test

    return testdata

def visualtests():
    """These are not formal tests -- but they should be eyeballed."""
    
    output = [ getfile(Gaussian,"basicGaussian03","dvb_gopt.out")[0],
               getfile(GAMESS,"basicPCGAMESS","dvb_gopt_a.out")[0],
               getfile(GAMESS,"basicGAMESS-US","dvb_gopt_a.out")[0],
               getfile(ADF,"basicADF2007.01","dvb_gopt.adfout")[0],
               getfile(Jaguar,"basicJaguar7.0", "dvb_gopt.out")[0],
               getfile(Molpro,"basicMolpro2006", "dvb_gopt.out", "dvb_gopt.out")[0],
             ]

    print "\n\nMO energies of optimised dvb"
    print "      ", "".join(["%-12s" % x for x in ['Gaussian03','PC-GAMESS','GAMESS-US','ADF2007.01','Jaguar7.0','Molpro2006']])
    print "HOMO", "   ".join(["%+9.4f" % x.moenergies[0][x.homos[0]] for x in output])
    print "LUMO", "   ".join(["%+9.4f" % x.moenergies[0][x.homos[0]+1] for x in output])
    print "H-L ", "   ".join(["%9.4f" % (x.moenergies[0][x.homos[0]+1]-x.moenergies[0][x.homos[0]],) for x in output])

def importName(modulename, name):
    """Import from a module whose name is determined at run-time.

    Taken from Python Cookbook 2nd ed O'Reilly Recipe 16.3.
    Additionally, also returns None if module does not habe attribute name.
    
    Inputs:
        modulename - name of the module
        name - name to be imported
    """

    try:
        module = __import__(modulename, globals(), locals(), [name])
    except ImportError:
        return None

    return getattr(module, name, None)

def testall(parserchoice=None, modules=test_modules):
    """Run all unittests in all modules."""

    # Make sure we are in the test directory of this script,
    #   so that getfile() can access the data files.
    curdir = os.path.abspath(os.curdir)
    destdir = os.path.dirname(__file__)
    if destdir:
        os.chdir(destdir)

    perpackage = {}
    errors = []
    
    for module in modules:

        testdata = gettestdata(module)
        
        if parserchoice:
            testdata = dict([ (x,y) for x,y in testdata.iteritems()
                              if y['parser']==parserchoice ])
                
        testnames = testdata.keys()
        testnames.sort()
        for name in testnames:

            path = '/'.join(testdata[name]["location"])
            program = testdata[name]["location"][0][5:]

            try:
                test = importName("test%s" %module, name)
            except:
                errors.append("ERROR: could not import %s from %s." %(name, module))
            else:
                print "\n**** test%s: %s ****" %(module, test.__doc__)
                parser = testdata[name]["parser"]
                location = testdata[name]["location"]
                test.data, test.logfile = getfile(eval(parser), *location)
                myunittest = unittest.makeSuite(test)
                a = unittest.TextTestRunner(verbosity=2).run(myunittest)
                l = perpackage.setdefault(program, [0, 0, 0, 0])
                l[0] += a.testsRun
                l[1] += len(a.errors)
                l[2] += len(a.failures)
                if hasattr(a, "skipped"):
                    l[3] += len(a.skipped)

    print "\n\n********* SUMMARY PER PACKAGE ****************"
    names = perpackage.keys()
    names.sort()
    total = [0, 0, 0, 0]
    print " "*14, "\t".join(["Total", "Passed", "Failed", "Errors", "Skipped"])
    for name in names:
        l = perpackage[name]
        print name.ljust(15), "%3d\t%3d\t%3d\t%3d\t%3d" % (l[0], l[0]-l[1]-l[2]-l[3], l[2], l[1], l[3])
        for i in range(4):
            total[i] += l[i]

    print "\n\n********* SUMMARY OF EVERYTHING **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d\tSKIPPED: %d" \
            %(total[0], total[0]-(total[1]+total[2]+total[3]), total[2], total[1], total[3])

    if errors:
        print "\n".join(errors)

    print "\n\n*** Visual tests ***"
    visualtests()
    
    # Return to the directory we started from.
    if destdir:
        os.chdir(curdir)


if __name__ == "__main__":
    parser = None
    if len(sys.argv)==2:
        parser = sys.argv[1]
    testall(parser)
