# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Run data tests for cclib."""

from __future__ import print_function

import importlib
import os
import sys
import unittest

import cclib


__filedir__ = os.path.realpath(os.path.dirname(__file__))


# We need this in Python3 for importing things from the same directory
# within the unit test files.
sys.path.append(os.path.join(__filedir__, 'data'))


parser_names = [
    "ADF", "DALTON", "GAMESS", "GAMESSUK", "Gaussian",
    "Jaguar", "Molpro", "NWChem", "ORCA", "Psi", "QChem",
]
all_parsers = {name: getattr(cclib.parser, name) for name in parser_names}


module_names = [
    "SP", "SPun", "GeoOpt", "Basis", "Core",    # Basic calculations.
    "MP", "CC", "CI", "TD", "TDun",             # Post-SCF calculations.
    "vib", "Scan", "BOMD",                      # Other property calculations.
]
all_modules = {tn: importlib.import_module('data.test' + tn) for tn in module_names}


def gettestdata():
    """Return a dict of the test file data."""

    testdatadir = os.path.dirname(os.path.realpath(__file__))
    lines = open(testdatadir + '/testdata').readlines()

    # Remove blank lines and those starting with '#'.
    lines = [line for line in lines if (line.strip() and line[0] != '#')]

    # Remove comment at end of lines (everything after a '#').
    lines = [line.split('#')[0] for line in lines]

    # Transform remaining lines into dictionaries.
    cols = [line.split() for line in lines]
    labels = ('module', 'parser', 'class', 'subdir', 'files')
    testdata = [dict(zip(labels, (c[0], c[1], c[2], c[3], c[4:]))) for c in cols]

    return testdata


def get_program_dir(parser_name):
    """Return a directory name given a parser name.

    In at least one case (GAMESS-UK) the directory is named differently.
    """

    if parser_name == "GAMESSUK":
        return "GAMESS-UK"
    return parser_name


def getdatafile(parser, subdir, *files, **kwds):
    """Returns a parsed logfile.

    Inputs:
        parser - a logfile parser class (subclass of LogFile)
        subdir - subdirectory containing data files (program version)
        *files - data filename(s)
        **kwds - currently accepts 'stream' keyword argument

    Outputs:
        data - the resulting data object
        logfile - the parser object used for parsing
    """

    # Convert any string into the parser object we will be using.
    if isinstance(parser, str):
        parser = all_parsers[parser]

    datadir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
    programdir = os.path.join(get_program_dir(parser.__name__), subdir)
    inputs = [os.path.join(datadir, programdir, fn) for fn in files]

    # We should be able to pass a list of length one here, but for some reason
    # this does not work with some parsers and we get errors.
    if len(inputs) == 1:
        inputs = inputs[0]

    stream = kwds.get('stream', sys.stdout)
    logfile = parser(inputs, logstream=stream)
    logfile.logger.setLevel(0)

    data = logfile.parse()
    return data, logfile


def ccdata_getattribute_with_coverage(self, attr):
    """A bookkeeping version of __getattribute__ for ccData objects."""
    if attr != '_attrlist' and attr in self._attrlist:
        if not hasattr(self, 'coverage'):
            self.coverage = {}
        self.coverage[attr] = self.coverage.get(attr, 0) + 1
    return object.__getattribute__(self, attr)


class DataSuite(object):
    """Suite containing data (logfile) tests in cclib.

    This is supposed to represent a single run of the entire data test suite in cclib or
    a subset of it. The main functions are to load data, run test cases in the data/
    subdirectory, and do some basic bookkeeping.
    """

    def __init__(self, argv=[], stream=sys.stdout):

        self.stream = stream
        self.argparse(argv)

        # Load the test data and filter with parsers and modules.
        self.testdata = gettestdata()
        self.testdata = [td for td in self.testdata if td['parser'] in self.parsers]
        self.testdata = [td for td in self.testdata if td['module'] in self.modules]

        # We want to gather the unit tests and results in several lists/dicts,
        # in order to easily generate summaries at the end.
        self.errors = []
        self.failures = []
        self.alltests = []
        self.perpackage = {p: [0, 0, 0, 0] for p in self.parsers}

    def argparse(self, argv):
        """Parse a list of command line arguments tfor anything relevant."""

        # These allow the parsers and modules tested to be filtered on the command line
        # with any number of arguments. No matching parsers/modules implies all of them.
        self.parsers = {p: all_parsers[p] for p in parser_names if p in argv} or all_parsers
        self.modules = {m: all_modules[m] for m in module_names if m in argv} or all_modules

        # These options modify the output and are used for Travis CI.
        self.status = "status" in argv or "--status" in argv
        self.terse = "terse" in argv or "--terse" in argv

    def testall(self):
        """Run all unittests in all modules.

        Run unit tests for all or a subset of parsers and modules. Arguments:
            stream - stream used for all output
        """

        stream_test = self.stream
        if self.terse:
            devnull = open(os.devnull, 'w')
            stream_test = devnull

        for td in self.testdata:

            module = self.modules[td['module']]
            parser = self.parsers[td['parser']]
            test = getattr(module, td['class'])

            description = "%s/%s: %s" % (td['subdir'], td['files'], test.__doc__)
            print("", file=stream_test)
            print("*** %s ***" % description, file=self.stream)

            test.data, test.logfile = getdatafile(parser, td['subdir'], *td['files'], stream=self.stream)

            # By overriding __getattribute__ temporarily with a custom method, we collect
            # coverage information for data attributes while the tests are run. This slightly
            # hacky approach is very convenient since it is self-contained and we don't
            # need to worry about it when writing the actual test cases.
            test.data.__class__.__getattribute__ = ccdata_getattribute_with_coverage

            # Here we actually run the tests for this line in testdata.
            myunittest = unittest.makeSuite(test)
            results = unittest.TextTestRunner(stream=stream_test, verbosity=2).run(myunittest)

            # We don't want to collect coverage stats beyond this point, so set __getattribute__
            # back to its original value. Note that we are setting the class method.
            test.data.__class__.__getattribute__ = object.__getattribute__

            self.perpackage[td['parser']][0] += results.testsRun
            self.perpackage[td['parser']][1] += len(results.errors)
            self.perpackage[td['parser']][2] += len(results.failures)
            self.perpackage[td['parser']][3] += len(getattr(results, 'skipped', []))

            self.alltests.append(test)
            self.errors.extend([description + "\n" + "".join(map(str, e)) for e in results.errors])
            self.failures.extend([description + "\n" + "".join(map(str, f)) for f in results.failures])

        if self.terse:
            devnull.close()

    def summary(self):
        """Prints a summary of the suite after it has been run."""

        if self.errors:
            print("\n********* SUMMARY OF ERRORS *********\n", file=self.stream)
            print("\n".join(self.errors), file=self.stream)

        if self.failures:
            print("\n********* SUMMARY OF FAILURES *********\n", file=self.stream)
            print("\n".join(self.failures), file=self.stream)

        print("\n********* SUMMARY PER PACKAGE ****************", file=self.stream)
        names = sorted(self.perpackage.keys())
        total = [0, 0, 0, 0]
        print(" "*14, "\t".join(["Total", "Passed", "Failed", "Errors", "Skipped"]), file=self.stream)

        fmt = "%3d\t%3d\t%3d\t%3d\t%3d"
        for name in names:
            l = self.perpackage[name]
            args = (l[0], l[0]-l[1]-l[2]-l[3], l[2], l[1], l[3])
            print(name.ljust(15), fmt % args, file=self.stream)
            for i in range(4):
                total[i] += l[i]

        print("\n********* SUMMARY OF EVERYTHING **************", file=self.stream)
        print("TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d\tSKIPPED: %d" \
                %(total[0], total[0]-(total[1]+total[2]+total[3]), total[2], total[1], total[3]), file=self.stream)

        if self.status and len(self.errors) > 0:
            sys.exit(1)

    def visualtests(self, stream=sys.stdout):
        """These are not formal tests -- but they should be eyeballed."""

        parsers_to_test = {
            'ADF2013.01' : getdatafile('ADF', "basicADF2013.01", "dvb_gopt.adfout")[0],
            'DALTON2015' : getdatafile('DALTON', "basicDALTON-2015", "dvb_gopt_ks.out")[0],
            'Firefly8.0' : getdatafile('GAMESS', "basicFirefly8.0", "dvb_gopt_a.out")[0],
            'Gaussian09' : getdatafile('Gaussian', "basicGaussian09", "dvb_gopt.out")[0],
            'GAMESS-US' : getdatafile('GAMESS', "basicGAMESS-US2012", "dvb_gopt_a.out")[0],
            'Jaguar8.0' : getdatafile('Jaguar', "basicJaguar8.3", "dvb_gopt_ks.out")[0],
            'Molpro2012' : getdatafile('Molpro', "basicMolpro2012", "dvb_gopt.log", "dvb_gopt.out")[0],
            'NWChem6.0' : getdatafile('NWChem', "basicNWChem6.0", "dvb_gopt_ks.out")[0],
            'ORCA3.0' : getdatafile('ORCA', "basicORCA3.0", "dvb_gopt.out")[0],
            'QChem4.2' : getdatafile('QChem', "basicQChem4.2", "dvb_gopt.out")[0],
        }
        parser_names = sorted(parsers_to_test.keys())
        output = [parsers_to_test[pn] for pn in parser_names]

        print("\n*** Visual tests ***", file=self.stream)
        print("MO energies of optimised dvb", file=self.stream)
        print("      ", "".join(["%-12s" % pn for pn in parser_names]), file=self.stream)
        print("HOMO", "   ".join(["%+9.4f" % out.moenergies[0][out.homos[0]] for out in output]), file=self.stream)
        print("LUMO", "   ".join(["%+9.4f" % out.moenergies[0][out.homos[0]+1] for out in output]), file=stream)
        print("H-L ", "   ".join(["%9.4f" % (out.moenergies[0][out.homos[0]+1]-out.moenergies[0][out.homos[0]],) for out in output]), file=self.stream)


if __name__ == "__main__":

    suite = DataSuite(sys.argv)
    suite.testall()
    suite.summary()
    suite.visualtests()
