#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Script for loading data from computational chemistry files."""


from __future__ import print_function

import getopt
import glob
import logging
import os.path
import sys
from functools import partial
from pprint import pprint

import numpy

from cclib.parser import ccData
from cclib.io import ccread, URL_PATTERN


MSG_USAGE = """\
Usage:  ccget <attribute> [<attribute>] <compchemlogfile> [<compchemlogfile>]
Try     ccget --help for more information\
"""

MSG_USAGE_LONG = """\
Usage:  ccget <attribute> [<attribute>] <compchemlogfile> [<compchemlogfile>]
    where <attribute> is one of the attributes to be parsed by cclib
    from each of the compchemlogfiles.
For a list of attributes available in a file, use --list (or -l):
    ccget --list <compchemlogfile>
To extract attributes available in a CJSON file, use --json (or -j):
    ccget --json  <attr> [<attr>]  <compchemlogfile>
To parse multiple files as one input stream, use --multi (or -m):
    ccget --multi <attr> [<attr>]  <cclogfile> <cclogfile> [<cclogfile>]
Additional options:
    -v or --verbose: more verbose parsing output (only errors by default)
    -u or --future: use experimental features (currently optdone_as_list)
    -f or --full: toggle full print behaviour for attributes\
"""

# These are the options ccget accepts and their one letter versions.
OPTS_LONG = ["help", "list", "json", "multi", "verbose", "future", "full"]
OPTS_SHORT = "hljmvuf"


# Set up options for pretty-printing output.
if sys.version_info < (3, 4):
    pprint = partial(pprint, width=120)
else:
    pprint = partial(pprint, width=120, compact=True)
numpy.set_printoptions(linewidth=120)


def ccget():
    """Parse files with cclib based on command line arguments."""

    # Parse the arguments and pass them to ccget, but print help information
    # and exit if it fails.
    try:
        optlist, arglist = getopt.getopt(sys.argv[1:], OPTS_SHORT, OPTS_LONG)
    except getopt.GetoptError:
        print(MSG_USAGE_LONG)
        sys.exit(1)

    future = False
    showattr = False
    cjsonfile = False
    multifile = False
    verbose = False
    full = False
    for opt, arg in optlist:
        if opt in ("-h", "--help"):
            print(MSG_USAGE_LONG)
            sys.exit()
        if opt in ("-l", "--list"):
            showattr = True
        if opt in ("-j", "--json"):
            cjsonfile = True
        if opt in ("-m", "--multi"):
            multifile = True
        if opt in ("-v", "--verbose"):
            verbose = True
        if opt in ("-u", "--future"):
            future = True
        if opt in ("-f", "--full"):
            full = True

    # Toggle full print behaviour for numpy arrays.
    if full:
        numpy.set_printoptions(threshold=numpy.nan)

    # We need at least one attribute and the filename, so two arguments, or
    # just one filename if we want to list attributes that can be extracted.
    # In multifile mode, we generally want at least two filenames, so the
    # expected number of arguments is a bit different.
    if not multifile:
        correct_number = (not showattr and len(arglist) > 1) or (showattr and len(arglist) > 0)
    else:
        correct_number = (not showattr and len(arglist) > 2) or (showattr and len(arglist) > 1)
    if not correct_number:
        print("The number of arguments does not seem to be correct.")
        print(MSG_USAGE)
        sys.exit(1)

    # Figure out which are the attribute names and which are the filenames or links.
    # Note that in Linux, the shell expands wild cards, but not so in Windows,
    # so try to do that here using glob.
    attrnames = []
    filenames = []
    for arg in arglist:
        if arg in ccData._attrlist:
            attrnames.append(arg)
        elif URL_PATTERN.match(arg) or os.path.isfile(arg):
            filenames.append(arg)
        else:
            wildcardmatches = glob.glob(arg)
            if wildcardmatches:
                filenames.extend(wildcardmatches)
            else:
                print("%s is neither a filename nor an attribute name." % arg)
                print(MSG_USAGE)
                sys.exit(1)

    # Since there is some ambiguity to the correct number of arguments, check
    # that there is at least one filename (or two in multifile mode), and also
    # at least one attribute to parse if the -l option was not passed.
    if len(filenames) == 0:
        print("No logfiles given")
        sys.exit(1)
    if multifile and len(filenames) == 1:
        print("Expecting at least two logfiles in multifile mode")
        sys.exit(1)
    if not showattr and len(attrnames) == 0:
        print("No attributes given")
        sys.exit(1)

    # This should be sufficient to correctly handle multiple files, that is to
    # run the loop below only once with all logfiles in the variable `filename`.
    # Although, perhaps it would be clearer to abstract the contents of the loop
    # into another function.
    if multifile:
        filenames = [filenames]

    # Now parse each file and print out the requested attributes.
    for filename in filenames:

        if multifile:
            name = ", ".join(filename[:-1]) + " and " + filename[-1]
        else:
            name = filename

        # The keyword dictionary are not used so much. but could be useful for
        # passing options downstream. For example, we might use --future for
        # triggering experimental or alternative behavior (as with optdone).
        kwargs = {}
        if verbose:
            kwargs['verbose'] = True
            kwargs['loglevel'] = logging.INFO
        else:
            kwargs['verbose'] = False
            kwargs['loglevel'] = logging.ERROR
        if future:
            kwargs['future'] = True
        if cjsonfile:
            kwargs['cjson'] = True

        print("Attempting to read %s" % name)
        data = ccread(filename, **kwargs)

        if data == None:
            print("Cannot figure out the format of '%s'" % name)
            print("Report this to the cclib development team if you think it is an error.")
            print("\n" + MSG_USAGE)
            sys.exit()

        if showattr:
            print("cclib can parse the following attributes from %s:" % name)
            if cjsonfile:
                for key in data:
                    print(key)
                break
            for attr in data._attrlist:
                if hasattr(data, attr):
                    print("  %s" % attr)
        else:
            invalid = False
            for attr in attrnames:
                if cjsonfile:
                    if attr in data:
                        print("%s:\n%s" % (attr, data[attr]))
                        continue
                else:
                    if hasattr(data, attr):
                        print(attr)
                        attr_val = getattr(data, attr)
                        # List of attributes to be printed with new lines
                        if attr in data._listsofarrays and full:
                            for val in attr_val:
                                pprint(val)
                        else:
                            pprint(attr_val)
                        continue

                print("Could not parse %s from this file." % attr)
                invalid = True
            if invalid:
                print(MSG_USAGE_LONG)


if __name__ == "__main__":

    ccget()
