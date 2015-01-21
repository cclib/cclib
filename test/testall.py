# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Aggregate and run all unit tests for cclib"""

from __future__ import print_function

import os
import sys
import unittest

from cclib.parser import (ADF, GAMESS, GAMESSUK, Gaussian,
                          Jaguar, Molpro, NWChem, ORCA, Psi, QChem)


# All supported parsers.
parsers = [ "ADF", "GAMESS", "GAMESSUK", "Gaussian", "Jaguar", "Molpro", "NWChem", "ORCA", "Psi", "QChem" ]

# The modules to be included in the global test testall().
test_modules = [ "SP", "SPun", "GeoOpt", "Basis", "Core",   # Basic calculations.
                 "MP", "CC", "CI", "TD", "TDun",            # Post-SCF calculations.
                 "vib", "Scan" ]                            # Other property calculations.


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
        **kwds - currently accepts 'stream' keyword argument
    
    Outputs:
        data - the resulting data object
        logfile - the parser object used for parsing
    """

    location = os.path.join(("..", "data", get_program_dir(parser.__name__)) + location)
    stream = kwds.get('stream', sys.stdout)

    # Construct the proper full path(s). Multiple paths will be in a list only if more than one
    # data file given. Presently, location contains only one subdirectory (basic*), so this is easy
    # since there are normally 5 elements in location.
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
    lines = [line for line in lines if (line.strip() and line[0] != '#')]

    # Remove comment at end of lines (everything after a '#').
    lines = [line.split('#')[0] for line in lines]

    # Split up each line into columns.
    lines = [line.split() for line in lines]

    # Filter for lines only for the given module.
    if module:
        lines = [line for line in lines if line[0] == module]

    # Each dictionary in this list contains the information needed to run one unit test.
    testdata = []
    for line in lines:
        test = {}
        test["module"] = line[0]
        test["parser"] = line[1]
        test["class"] = line[2]
        test["location"] = line[3:]
        testdata.append(test)

    return testdata


def visualtests(stream=sys.stdout):
    """These are not formal tests -- but they should be eyeballed."""

    parsers_to_test = {
        'ADF2013.01' : getfile(ADF, "basicADF2013.01", "dvb_gopt.adfout")[0],
        'Firefly8.0' : getfile(GAMESS, "basicFirefly8.0", "dvb_gopt_a.out")[0],
        'Gaussian09' : getfile(Gaussian, "basicGaussian09", "dvb_gopt.out")[0],
        'GAMESS-US' : getfile(GAMESS, "basicGAMESS-US2012", "dvb_gopt_a.out")[0],
        'Jaguar8.0' : getfile(Jaguar, "basicJaguar8.3", "dvb_gopt_ks.out")[0],
        'Molpro2012' : getfile(Molpro, "basicMolpro2012", "dvb_gopt.log", "dvb_gopt.out")[0],
        'NWChem6.0' : getfile(NWChem, "basicNWChem6.0", "dvb_gopt_ks.out")[0],
        'ORCA3.0' : getfile(ORCA, "basicORCA3.0", "dvb_gopt.out")[0],
        'QChem4.2' : getfile(QChem, "basicQChem4.2", "dvb_gopt.out")[0],
    }
    parser_names = sorted(parsers_to_test.keys())
    output = [parsers_to_test[pn] for pn in parser_names]

    print("\n*** Visual tests ***", file=stream)
    print("MO energies of optimised dvb", file=stream)
    print("      ", "".join(["%-12s" % pn for pn in parser_names]), file=stream)
    print("HOMO", "   ".join(["%+9.4f" % out.moenergies[0][out.homos[0]] for out in output]), file=stream)
    print("LUMO", "   ".join(["%+9.4f" % out.moenergies[0][out.homos[0]+1] for out in output]), file=stream)
    print("H-L ", "   ".join(["%9.4f" % (out.moenergies[0][out.homos[0]+1]-out.moenergies[0][out.homos[0]],) for out in output]), file=stream)


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


def testall(parsers=parsers, modules=test_modules, status=False, terse=False, stream=sys.stdout):
    """Run all unittests in all modules.

    Run unit tests for all or a subset of parsers and modules. Arguments:
        parsers - list of parsers to test (all by default)
        modules - list of modules to test (all by default)
        status - exists with error status when any errors in tests found
        terse - do not print output from tests
        stream - stream used for all output
    """

    # Make sure we are in the test directory of this script, so that getfile()
    # can access the data files.
    curdir = os.path.abspath(os.curdir)
    destdir = os.path.dirname(__file__)
    if destdir:
        os.chdir(destdir)

    # We want to gather the unit tests and results in several lists/dicts,
    # in order to easily generate summaries at the end.
    errors = []
    failures = []
    alltests = []
    perpackage = {}

    stream_test = stream
    if terse:
        devnull = open(os.devnull, 'w')
        stream_test = devnull

    for module in modules:

        testdata = gettestdata(module)

        # Filter the test data by parser, which has an effect only if a list of requested parsers
        # is passed to this function othet than the default. We used to assume here that all unit tests
        # for a given test class use the same parser, but that is no longer true. We now use generic
        # test classes in many cases, in order to reduce the total number of classes needed.
        testdata = [t for t in testdata if t['parser'] in parsers]

        for test_instance in testdata:

            name = test_instance['class']
            path = '/'.join(test_instance["location"])
            program = test_instance["location"][0][5:]
            fname = test_instance["location"][-1]

            try:
                test = importName("test%s" %module, name)
            except:
                errors.append("ERROR: could not import %s from %s." %(name, module))
            else:
                description = "%s/%s: %s" %(program, fname, test.__doc__)
                print("", file=stream_test)
                print("**** %s ****" % description, file=stream)

                parser = test_instance["parser"]
                location = test_instance["location"]
                test.data, test.logfile = getfile(eval(parser), *location, stream=stream)

                myunittest = unittest.makeSuite(test)
                a = unittest.TextTestRunner(stream=stream_test, verbosity=2).run(myunittest)

                l = perpackage.setdefault(program, [0, 0, 0, 0])
                l[0] += a.testsRun
                l[1] += len(a.errors)
                l[2] += len(a.failures)

                if hasattr(a, "skipped"):
                    l[3] += len(a.skipped)
                alltests.append(test)
                errors.extend([description+"\n"+"".join(map(str,e)) for e in a.errors])
                failures.extend([description+"\n"+"".join(map(str,f)) for f in a.failures])

    if terse:
        devnull.close()

    if errors:
        print("\n********* SUMMARY OF ERRORS *********\n", file=stream)
        print("\n".join(errors), file=stream)

    if failures:
        print("\n********* SUMMARY OF FAILURES *********\n", file=stream)
        print("\n".join(failures), file=stream)

    print("\n********* SUMMARY PER PACKAGE ****************", file=stream)
    names = sorted(perpackage.keys())
    total = [0, 0, 0, 0]
    print(" "*14, "\t".join(["Total", "Passed", "Failed", "Errors", "Skipped"]), file=stream)

    fmt = "%3d\t%3d\t%3d\t%3d\t%3d"
    for name in names:
        l = perpackage[name]
        args = (l[0], l[0]-l[1]-l[2]-l[3], l[2], l[1], l[3])
        print(name.ljust(15), fmt % args, file=stream)
        for i in range(4):
            total[i] += l[i]

    print("\n********* SUMMARY OF EVERYTHING **************", file=stream)
    print("TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d\tSKIPPED: %d" \
            %(total[0], total[0]-(total[1]+total[2]+total[3]), total[2], total[1], total[3]), file=stream)

    if destdir:
        os.chdir(curdir)

    if status and len(errors) > 0:
        sys.exit(1)

    return alltests


if __name__ == "__main__":

    # These allow the parsers and modules tested to be filtered on the command line
    # with any number of arguments. No matching parsers/modules implies all of them.
    parsers = [p for p in parsers if p in sys.argv] or parsers
    modules = [m for m in test_modules if m in sys.argv] or test_modules

    # These options are used for Travis CI.
    status = "status" in sys.argv or "--status" in sys.argv
    terse = "terse" in sys.argv or "--terse" in sys.argv

    tests = testall(parsers, modules, status, terse)
    visualtests()
