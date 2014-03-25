# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

from __future__ import print_function
import os
import sys
import unittest

from cclib.parser import ADF, GAMESS, GAMESSUK, Gaussian, Jaguar, Molpro, ORCA


# All supported parsers.
parsers = [ "ADF", "GAMESS", "GAMESSUK", "Gaussian", "Jaguar", "Molpro", "ORCA" ]

# The modules to be included in the global test testall().
test_modules = [ "SP", "SPun", "GeoOpt", "Basis", "Core",   # Basic calculations.
                 "MP", "CC", "CI", "TD", "TDun",            # Post-SCF calculations.
                 "vib" ]                                    # Other property calculations.


def get_program_dir(parser_name):
    """In one case the directory is names differently than the parser."""
    if parser_name == "GAMESSUK":
        return "GAMESS-UK"
    else:
        return parser_name

def getfile(parser, *location, **kwds):
    """Returns a parsed logfile.
    
    Inputs:
        parser - a logfile parser class (subclass of LogFile)
        *location - subdirectory and data filename(s)
        **kwds - accepts 'stream' keyword argument
    
    Outputs:
        data - the resulting data object
        logfile - the parser object used for parsing
    """

    location = os.path.join(("..", "data", get_program_dir(parser.__name__)) + location)
    stream = kwds.get('stream', sys.stdout)

    # Now construct the proper full path(s).
    # Multiple paths will be in a list only if more than one data file given.
    # Presently, location contains only one subdirectory (basic*),
    #   so this is easy since there are normally 5 elements in location.
    if len(location) == 5:
        filename = os.path.join(*location)
    else:
        filename = [os.path.join(*(location[:4]+location[n:n+1])) for n in range(4,len(location))]

    logfile = parser(filename, logstream=stream)
    logfile.logger.setLevel(0)
    data = logfile.parse()
    
    return data, logfile


def gettestdata(module=None):
    """Returns a dict of test files for a given module."""

    testdatadir = os.path.dirname(os.path.realpath(sys.argv[0]))
    lines = open(testdatadir+'/testdata').readlines()

    # Remove blank lines and those starting with '#'.
    lines = [line.split() for line in lines if (line.strip() and line[0] != '#')]
    
    # Filter for lines only for the given module.
    if module:
        lines = [line for line in lines if line[0] == module]

    # This dictionary contains information needed to run unit tests, with one entry
    # for each unit test class used. Since there may in principle be multiple tests
    # for a test class (files from different program versions, for example), the values
    # in this dictionary are themselves lists of dictionaries.
    testdata = {}
    for line in lines:
        test = {}
        test["module"] = line[0]
        test["parser"] = line[1]
        test["location"] = line[3:]
        testclass = line[2]
        if testclass not in testdata:
            testdata[testclass] = []
        testdata[testclass].append(test)

    return testdata


def visualtests(stream=sys.stdout):
    """These are not formal tests -- but they should be eyeballed."""
    
    output = [ getfile(Gaussian,"basicGaussian03","dvb_gopt.out")[0],
               getfile(GAMESS,"basicPCGAMESS","dvb_gopt_a.out")[0],
               getfile(GAMESS,"basicGAMESS-US2012","dvb_gopt_a.out")[0],
               getfile(ADF,"basicADF2007.01","dvb_gopt.adfout")[0],
               getfile(Jaguar,"basicJaguar7.0", "dvb_gopt.out")[0],
               getfile(Molpro,"basicMolpro2006", "dvb_gopt.out", "dvb_gopt.out")[0],
             ]

    print("\n\nMO energies of optimised dvb", file=stream)
    print("      ", "".join(["%-12s" % x for x in ['Gaussian03','PC-GAMESS','GAMESS-US','ADF2007.01','Jaguar7.0','Molpro2006']]), file=stream)
    print("HOMO", "   ".join(["%+9.4f" % x.moenergies[0][x.homos[0]] for x in output]), file=stream)
    print("LUMO", "   ".join(["%+9.4f" % x.moenergies[0][x.homos[0]+1] for x in output]), file=stream)
    print("H-L ", "   ".join(["%9.4f" % (x.moenergies[0][x.homos[0]+1]-x.moenergies[0][x.homos[0]],) for x in output]), file=stream)


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


def testall(parserchoice=parsers, modules=test_modules, stream=sys.stdout):
    """Run all unittests in all modules."""

    # Make sure we are in the test directory of this script,
    #   so that getfile() can access the data files.
    curdir = os.path.abspath(os.curdir)
    destdir = os.path.dirname(__file__)
    if destdir:
        os.chdir(destdir)

    perpackage = {}
    errors = []
    
    alltests = []
    for module in modules:

        testdata = gettestdata(module)

        # Filter the test data to use if a list of requested parsers is passed
        # to this function, assuming that all unit tests for a given test class
        # use the same parser (use the first unit test to check).
        if parserchoice:
            testdata = dict([(x,y) for x,y in testdata.items() if y[0]['parser'] in parserchoice ])
                
        testnames = sorted(testdata.keys())
        for name in testnames:

            for test_instance in testdata[name]:

                path = '/'.join(test_instance["location"])
                program = test_instance["location"][0][5:]

                try:
                    test = importName("test%s" %module, name)
                except:
                    errors.append("ERROR: could not import %s from %s." %(name, module))
                else:
                    print("\n**** test%s (%s): %s ****" %(module, program, test.__doc__), file=stream)
                    parser = test_instance["parser"]
                    location = test_instance["location"]
                    test.data, test.logfile = getfile(eval(parser), *location, stream=stream)
                    myunittest = unittest.makeSuite(test)
                    a = unittest.TextTestRunner(stream=stream, verbosity=2).run(myunittest)
                    l = perpackage.setdefault(program, [0, 0, 0, 0])
                    l[0] += a.testsRun
                    l[1] += len(a.errors)
                    l[2] += len(a.failures)
                    if hasattr(a, "skipped"):
                        l[3] += len(a.skipped)
                    alltests.append(test)

    print("\n\n********* SUMMARY PER PACKAGE ****************", file=stream)
    names = sorted(perpackage.keys())
    total = [0, 0, 0, 0]
    print(" "*14, "\t".join(["Total", "Passed", "Failed", "Errors", "Skipped"]), file=stream)
    for name in names:
        l = perpackage[name]
        print(name.ljust(15), "%3d\t%3d\t%3d\t%3d\t%3d" % (l[0], l[0]-l[1]-l[2]-l[3], l[2], l[1], l[3]), file=stream)
        for i in range(4):
            total[i] += l[i]

    print("\n\n********* SUMMARY OF EVERYTHING **************", file=stream)
    print("TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d\tSKIPPED: %d" \
            %(total[0], total[0]-(total[1]+total[2]+total[3]), total[2], total[1], total[3]), file=stream)

    if errors:
        print("\n".join(errors))

    print("\n\n*** Visual tests ***", file=stream)
    visualtests(stream=stream)
    
    # Return to the directory we started from.
    if destdir:
        os.chdir(curdir)

    return alltests


if __name__ == "__main__":
    chosen_parsers = [p for p in parsers if p in sys.argv] or parsers
    chosen_modules = [m for m in test_modules if m in sys.argv] or test_modules
    tests = testall(chosen_parsers, chosen_modules)
