# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Run parser tests for cclib."""

from __future__ import print_function

import importlib
import os
import sys
import unittest

import cclib


parser_names = [
    "ADF", "DALTON", "GAMESS", "GAMESSUK", "Gaussian",
    "Jaguar", "Molpro", "NWChem", "ORCA", "Psi", "QChem",
]
all_parsers = {name: getattr(cclib.parser, name) for name in parser_names}


module_names = [
    "SP", "SPun", "GeoOpt", "Basis", "Core",    # Basic calculations.
    "MP", "CC", "CI", "TD", "TDun",             # Post-SCF calculations.
    "vib", "Scan",                              # Other property calculations.
]
all_modules = {tn: importlib.import_module('data.test' + tn) for tn in module_names}
    

def get_program_dir(parser_name):
    """Return a directory name given a parser name.

    In at least one case (GAMESS-UK) the directory is named differently.
    """

    if parser_name == "GAMESSUK":
        return "GAMESS-UK"
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

    if isinstance(parser, str):
        parser = all_parsers[parser]

    location = os.path.join(("..", "data", get_program_dir(parser.__name__)) + location)
    stream = kwds.get('stream', sys.stdout)

    # Construct the proper full path(s). Multiple paths will be in a list only if more than one
    # data file given. Presently, location contains only one subdirectory (basic*), so this is easy
    # since there are normally 5 elements in location.
    if len(location) == 5:
        filename = os.path.join(*location)
    else:
        filename = [os.path.join(*(location[:4]+location[n:n+1])) for n in range(4, len(location))]

    logfile = parser(filename, logstream=stream)
    logfile.logger.setLevel(0)
    data = logfile.parse()

    return data, logfile


class DataSuite(object):
    """Suite containing data (logfile) tests in cclib.

    This is supposed to represent a single run of the entire data test suite in cclib or
    a subset of it. The main functions are to load data, run test cases in the data/
    subdirectory, and do some basic bookkeeping.
    """

    def __init__(self, argv=[]):

        if argv:
            self.argparse(argv)

    def argparse(self, argv):
        """Parse a list of command line arguments tfor anything relevant."""

        # These allow the parsers and modules tested to be filtered on the command line
        # with any number of arguments. No matching parsers/modules implies all of them.
        self.parsers = {p: all_parsers[p] for p in parser_names if p in argv} or all_parsers
        self.modules = {m: all_modules[m] for m in module_names if m in argv} or all_modules

        # These options modify the output and are used for Travis CI.
        self.status = "status" in argv or "--status" in argv
        self.terse = "terse" in argv or "--terse" in argv

    def gettestdata(self, module=None):
        """Return a dict of test files for a given module."""

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

    def testall(self, stream=sys.stdout):
        """Run all unittests in all modules.

        Run unit tests for all or a subset of parsers and modules. Arguments:
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
        self.errors = []
        self.failures = []
        self.alltests = []
        self.perpackage = {}

        stream_test = stream
        if self.terse:
            devnull = open(os.devnull, 'w')
            stream_test = devnull

        for module in self.modules:

            testdata = self.gettestdata(module)

            # Filter the test data by parser, which has an effect only if a list of requested parsers
            # is passed to this function othet than the default. We used to assume here that all unit tests
            # for a given test class use the same parser, but that is no longer true. We now use generic
            # test classes in many cases, in order to reduce the total number of classes needed.
            testdata = [t for t in testdata if t['parser'] in self.parsers]

            for test_instance in testdata:

                name = test_instance['class']
                path = '/'.join(test_instance["location"])
                program = test_instance["location"][0][5:]
                fname = test_instance["location"][-1]

                try:
                    test = getattr(self.modules[module], name)
                except:
                    self.errors.append("ERROR: could not import %s from %s." %(name, module))
                else:
                    description = "%s/%s: %s" %(program, fname, test.__doc__)
                    print("", file=stream_test)
                    print("**** %s ****" % description, file=stream)

                    parser = test_instance["parser"]
                    location = test_instance["location"]
                    print(parser)
                    test.data, test.logfile = getfile(self.parsers[parser], *location, stream=stream)

                    myunittest = unittest.makeSuite(test)
                    a = unittest.TextTestRunner(stream=stream_test, verbosity=2).run(myunittest)

                    l = self.perpackage.setdefault(program, [0, 0, 0, 0])
                    l[0] += a.testsRun
                    l[1] += len(a.errors)
                    l[2] += len(a.failures)

                    if hasattr(a, "skipped"):
                        l[3] += len(a.skipped)
                    self.alltests.append(test)
                    self.errors.extend([description+"\n"+"".join(map(str,e)) for e in a.errors])
                    self.failures.extend([description+"\n"+"".join(map(str,f)) for f in a.failures])

        if self.terse:
            devnull.close()

        if self.errors:
            print("\n********* SUMMARY OF ERRORS *********\n", file=stream)
            print("\n".join(self.errors), file=stream)

        if self.failures:
            print("\n********* SUMMARY OF FAILURES *********\n", file=stream)
            print("\n".join(self.failures), file=stream)

        print("\n********* SUMMARY PER PACKAGE ****************", file=stream)
        names = sorted(self.perpackage.keys())
        total = [0, 0, 0, 0]
        print(" "*14, "\t".join(["Total", "Passed", "Failed", "Errors", "Skipped"]), file=stream)

        fmt = "%3d\t%3d\t%3d\t%3d\t%3d"
        for name in names:
            l = self.perpackage[name]
            args = (l[0], l[0]-l[1]-l[2]-l[3], l[2], l[1], l[3])
            print(name.ljust(15), fmt % args, file=stream)
            for i in range(4):
                total[i] += l[i]

        print("\n********* SUMMARY OF EVERYTHING **************", file=stream)
        print("TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d\tSKIPPED: %d" \
                %(total[0], total[0]-(total[1]+total[2]+total[3]), total[2], total[1], total[3]), file=stream)

        if destdir:
            os.chdir(curdir)

        if self.status and len(self.errors) > 0:
            sys.exit(1)

        return self.alltests

    def visualtests(self, stream=sys.stdout):
        """These are not formal tests -- but they should be eyeballed."""

        # Make sure we are in the test directory of this script, so that getfile()
        # can access the data files.
        curdir = os.path.abspath(os.curdir)
        destdir = os.path.dirname(__file__)
        if destdir:
            os.chdir(destdir)

        parsers_to_test = {
            'ADF2013.01' : getfile('ADF', "basicADF2013.01", "dvb_gopt.adfout")[0],
            'Firefly8.0' : getfile('GAMESS', "basicFirefly8.0", "dvb_gopt_a.out")[0],
            'Gaussian09' : getfile('Gaussian', "basicGaussian09", "dvb_gopt.out")[0],
            'GAMESS-US' : getfile('GAMESS', "basicGAMESS-US2012", "dvb_gopt_a.out")[0],
            'Jaguar8.0' : getfile('Jaguar', "basicJaguar8.3", "dvb_gopt_ks.out")[0],
            'Molpro2012' : getfile('Molpro', "basicMolpro2012", "dvb_gopt.log", "dvb_gopt.out")[0],
            'NWChem6.0' : getfile('NWChem', "basicNWChem6.0", "dvb_gopt_ks.out")[0],
            'ORCA3.0' : getfile('ORCA', "basicORCA3.0", "dvb_gopt.out")[0],
            'QChem4.2' : getfile('QChem', "basicQChem4.2", "dvb_gopt.out")[0],
        }
        parser_names = sorted(parsers_to_test.keys())
        output = [parsers_to_test[pn] for pn in parser_names]

        print("\n*** Visual tests ***", file=stream)
        print("MO energies of optimised dvb", file=stream)
        print("      ", "".join(["%-12s" % pn for pn in parser_names]), file=stream)
        print("HOMO", "   ".join(["%+9.4f" % out.moenergies[0][out.homos[0]] for out in output]), file=stream)
        print("LUMO", "   ".join(["%+9.4f" % out.moenergies[0][out.homos[0]+1] for out in output]), file=stream)
        print("H-L ", "   ".join(["%9.4f" % (out.moenergies[0][out.homos[0]+1]-out.moenergies[0][out.homos[0]],) for out in output]), file=stream)

        if destdir:
            os.chdir(curdir)


if __name__ == "__main__":

    suite = DataSuite(sys.argv)
    suite.testall()
    suite.visualtests()