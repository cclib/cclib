#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import argparse
import logging
import os.path
import sys

from cclib.parser import ccData
from cclib.io import ccopen
from cclib.io import ccwrite


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('outputtype',
                        choices=('json', 'cjson', 'cml', 'xyz', 'molden', 'wfx'),
                        help='the output format to write (json/cjson are identical)')
    parser.add_argument('compchemlogfile',
                        nargs='+',
                        help='one or more computational chemistry output files to parse and convert')

    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='more verbose parsing output (only errors by default)')

    parser.add_argument('-g', '--ghost',
                        type=str,
                        default=None,
                        help='Symbol to use for ghost atoms')

    parser.add_argument('-t', '--terse',
                        action='store_true',
                        help='CJSON by default is not indented for readability, saves space (indented for readability\'s sake)')

    parser.add_argument('-u', '--future',
                        action='store_true',
                        help='use experimental features (currently optdone_as_list)')

    parser.add_argument('-i', '--index',
                        type=int,
                        default=None,
                        help='optional zero-based index for which structure to extract')

    args = parser.parse_args()

    outputtype = args.outputtype
    filenames = args.compchemlogfile
    verbose = args.verbose
    terse = args.terse
    future = args.future
    index = args.index
    ghost = args.ghost

    for filename in filenames:

        # We might want to use this option in the near future.
        ccopen_kwargs = dict()
        if future:
            ccopen_kwargs['future'] = True

        print(f"Attempting to parse {filename}")
        log = ccopen(filename, **ccopen_kwargs)

        if not log:
            print(
                f"Cannot figure out what type of computational chemistry output file '{filename}' is."
            )
            print("Report this to the cclib development team if you think this is an error.")
            sys.exit()

        if verbose:
            log.logger.setLevel(logging.INFO)
        else:
            log.logger.setLevel(logging.ERROR)
        data = log.parse()

        print(f"cclib can parse the following attributes from {filename}:")
        hasattrs = [f"  {attr}" for attr in ccData._attrlist if hasattr(data, attr)]
        print("\n".join(hasattrs))

        # Write out to disk.
        outputdest = '.'.join([os.path.splitext(os.path.basename(filename))[0], outputtype])
        ccwrite_kwargs = dict()
        if future:
            ccwrite_kwargs['future'] = True
        if ghost:
            ccwrite_kwargs['ghost'] = ghost
        # For XYZ files, write the last geometry unless otherwise
        # specified.
        if not index:
            index = -1
        ccwrite_kwargs['jobfilename'] = filename

        # The argument terse presently is only applicable to
        # CJSON/JSON formats
        ccwrite(data, outputtype, outputdest, indices=index, terse=terse,
                **ccwrite_kwargs)


if __name__ == "__main__":
    main()
