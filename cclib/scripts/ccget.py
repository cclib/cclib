#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Script for loading data from computational chemistry files."""


import glob
import logging
import os.path
import difflib
from functools import partial
from pprint import pprint

# This is needed for testing purposes only.
import sys

import numpy

from cclib.parser import ccData
from cclib.io import ccread, URL_PATTERN


# Set up options for pretty-printing output.
pprint = partial(pprint, width=120, compact=True)
numpy.set_printoptions(linewidth=120)


def ccget():
    """Parse files with cclib based on command line arguments."""

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "attribute_or_compchemlogfile", nargs="+",
        help="one or more attributes to be parsed from one ore more logfiles",
    )

    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        "--list", "-l",
        action="store_true",
        help="print a list of attributes available in each file",
    )
    group.add_argument(
        "--json", "-j",
        action="store_true",
        help="the given logfile is in CJSON format",
    )
    group.add_argument(
        "--multi", "-m",
        action="store_true",
        help="parse multiple input files as one input stream",
    )

    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="more verbose parsing output (only errors by default)",
    )
    parser.add_argument(
        "--future", "-u",
        action="store_true",
        help="use experimental features (currently optdone_as_list)",
    )
    parser.add_argument(
        "--full", "-f",
        action="store_true",
        help="toggle full print behaviour for attributes",
    )

    args = parser.parse_args()

    arglist = args.attribute_or_compchemlogfile
    showattr = args.list
    cjsonfile = args.json
    multifile = args.multi
    verbose = args.verbose
    future = args.future
    full = args.full

    # Toggle full print behaviour for numpy arrays.
    if full:
        numpy.set_printoptions(threshold=numpy.inf)

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
        parser.print_usage()
        parser.exit(1)

    # Figure out which are the attribute names and which are the filenames or links.
    # Note that in Linux, the shell expands wild cards, but not so in Windows,
    # so try to do that here using glob.
    attrnames = []
    filenames = []
    for arg in arglist:
        if arg not in ccData._attrlist:
            fuzzy_attr = difflib.get_close_matches(arg, ccData._attrlist, n=1, cutoff=0.85)
            if len(fuzzy_attr) > 0:
                fuzzy_attr = fuzzy_attr[0]
                logging.warning(
                    f"Attribute '{arg}' not found, but attribute '{fuzzy_attr}' is close. "
                    f"Using '{fuzzy_attr}' instead."
                )
                arg = fuzzy_attr
        if arg in ccData._attrlist:
            attrnames.append(arg)
        elif URL_PATTERN.match(arg) or os.path.isfile(arg):
            filenames.append(arg)
        else:
            wildcardmatches = glob.glob(arg)
            if wildcardmatches:
                filenames.extend(wildcardmatches)
            else:
                print(f"{arg} is neither a filename nor an attribute name.")
                parser.print_usage()
                parser.exit(1)

    # Since there is some ambiguity to the correct number of arguments, check
    # that there is at least one filename (or two in multifile mode), and also
    # at least one attribute to parse if the -l option was not passed.
    if len(filenames) == 0:
        print("No logfiles given")
        parser.exit(1)
    if multifile and len(filenames) == 1:
        print("Expecting at least two logfiles in multifile mode")
        parser.exit(1)
    if not showattr and len(attrnames) == 0:
        print("No attributes given")
        parser.exit(1)

    # This should be sufficient to correctly handle multiple files, that is to
    # run the loop below only once with all logfiles in the variable `filename`.
    # Although, perhaps it would be clearer to abstract the contents of the loop
    # into another function.
    if multifile:
        filenames = [filenames]

    # Now parse each file and print out the requested attributes.
    for filename in filenames:

        if multifile:
            name = f"{', '.join(filename[:-1])} and {filename[-1]}"
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

        print(f"Attempting to read {name}")
        data = ccread(filename, **kwargs)

        if data is None:
            print(f"Cannot figure out the format of '{name}'")
            print("Report this to the cclib development team if you think it is an error.")
            print(f"\n{parser.format_usage()}")
            parser.exit(1)

        if showattr:
            print(f"cclib can parse the following attributes from {name}:")
            if cjsonfile:
                for key in data:
                    print(key)
                break
            for attr in data._attrlist:
                if hasattr(data, attr):
                    print(f"  {attr}")
        else:
            invalid = False
            for attr in attrnames:
                if cjsonfile:
                    if attr in data:
                        print(f"{attr}:\n{data[attr]}")
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

                print(f"Could not parse {attr} from this file.")
                invalid = True
            if invalid:
                parser.print_help()


if __name__ == "__main__":

    ccget()
