# -*- coding: utf-8 -*-

"""Generate the coverage.rst and coverage.rst files from test
results."""

from __future__ import print_function

import os
import sys

from docs_common import check_cclib

# Import cclib and check we are using the version from a subdirectory.
import cclib
check_cclib(cclib)


def generate_coverage():
    """Generate a string containing a reStructuredTest table
    representation of which parsers support which attributes, based on
    test results.
    """
    lines = []

    # Change directory to where tests are and add it to the path. Because there are
    # separate directories for different branches/versions, and we use a symlink to
    # point to the one we want, we need to test the real path this link resolves to.
    if "cclib_prod" in os.path.realpath('cclib'):
        testpath = "_build/cclib_prod"
    else:
        assert "cclib_dev" in os.path.realpath('cclib')
        testpath = "_build/cclib_dev"

    os.chdir(testpath)

    thispath = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(1, thispath)

    from test.test_data import (all_modules, all_parsers, parser_names, DataSuite)
    import inspect
    ds_args = inspect.getargspec(DataSuite.__init__).args
    logpath = thispath + "/coverage.tests.log"
    try:
        with open(logpath, "w") as flog:
            stdout_backup = sys.stdout
            sys.stdout = flog
            alltests = {}
            for p in parser_names:
                assert 'parsers' in ds_args
                suite = DataSuite(parsers={p: all_parsers[p]}, modules=all_modules, stream=flog)
                suite.testall()
                alltests[p] = [{'data': t.data} for t in suite.alltests]
            sys.stdout = stdout_backup
    except Exception as e:
        print("Unit tests did not run correctly. Check log file for errors:")
        with open(logpath) as fh:
            print(fh.read())
        print(e)
        sys.exit(1)

    ncols = len(parser_names) + 1
    colwidth = 20
    colfmt = "%%-%is" % colwidth
    dashes = ("=" * (colwidth - 1) + " ") * ncols

    lines.append(dashes)
    lines.append(colfmt * ncols % tuple(["attributes"] + parser_names))
    lines.append(dashes)

    # Eventually we want to move this to cclib, too.
    not_applicable = {
        'ADF' : ['aonames', 'ccenergies', 'mpenergies'],
        'DALTON' : ['fonames', 'fooverlaps', 'fragnames', 'frags'],
        'GAMESS' : ['fonames', 'fooverlaps', 'fragnames', 'frags'],
        'GAMESSUK' : ['fonames', 'fooverlaps', 'fragnames', 'frags'],
        'Gaussian' : ['fonames', 'fooverlaps', 'fragnames', 'frags'],
        'Jaguar' : ['fonames', 'fooverlaps', 'fragnames', 'frags'],
        'Molpro' : ['fonames', 'fooverlaps', 'fragnames', 'frags'],
        'NWChem' : ['fonames', 'fooverlaps', 'fragnames', 'frags'],
        'ORCA' : ['fonames', 'fooverlaps', 'fragnames', 'frags'],
        'Psi' : ['fonames', 'fooverlaps', 'fragnames', 'frags'],
        'QChem' : ['fonames', 'fooverlaps', 'fragnames', 'frags'],
    }
    not_possible = {
        'Psi' : ['aooverlaps', 'vibirs'],
        'QChem' : ['aooverlaps', 'etrotats'],
    }

    # For each attribute, get a list of Boolean values for each parser that flags
    # if it has been parsed by at least one unit test. Substitute an OK sign or
    # T/D appropriately, with the exception of attributes that have been explicitely
    # designated as N/A.
    attributes = sorted(cclib.parser.data.ccData._attrlist)
    for attr in attributes:
        parsed = [any([attr in t['data'].__dict__ for t in alltests[p]]) for p in parser_names]
        for ip, p in enumerate(parsed):
            if p:
                parsed[ip] = "√"
            else:
                if attr in not_applicable.get(parser_names[ip], []):
                    parsed[ip] = "N/A"
                elif attr in not_possible.get(parser_names[ip], []):
                    parsed[ip] = "N/P"
                else:
                    parsed[ip] = "T/D"
        lines.append(colfmt*ncols % tuple(["`%s`_" % attr] + parsed))

    lines.append(dashes)
    lines.append("")

    for attr in attributes:
        lines.append(".. _`%s`: data_notes.html#%s" % (attr, attr))

    return "\n".join(lines)

if __name__ == "__main__":
    print(generate_coverage())
