#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
"""Script for writing data tables from computational chemistry files."""

import argparse
import os.path
import sys

from cclib.io import ccopen
from cclib.io import ccframe
from cclib.parser.utils import find_package

_has_pandas = find_package("pandas")
if _has_pandas:
    import pandas as pd


def process_logfiles(filenames, output, identifier):
    df = ccframe([ccopen(path) for path in filenames])
    if output is not None:
        outputtype = os.path.splitext(os.path.basename(output))[1][1:]
        if not outputtype:
            raise RuntimeWarning(
                "The output type could not be determined from the given path, "
                "not writing DataFrame to disk"
            )

        if outputtype in {'csv'}:
            df.to_csv(output, mode='w')
        elif outputtype in {'h5', 'hdf', 'hdf5'}:
            df.to_hdf(output, mode='w', key=identifier)
        elif outputtype in {'json'}:
            df.to_json(output)
        elif outputtype in {'pickle', 'pkl'}:
            df.to_pickle(output)
        elif outputtype in {'xlsx'}:
            writer = pd.ExcelWriter(output, mode='w')
            # This overwrites previous sheets
            # (see https://stackoverflow.com/a/42375263/4039050)
            df.to_excel(writer, sheet_name=identifier)
            writer.save()
    else:
        print(df)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-O', '--output',
                        help=('the output document to write, including an '
                              'extension supported by pandas '
                              '(csv, h5/hdf/hdf5, json, pickle/pkl, xlsx)'))
    parser.add_argument('compchemlogfiles', metavar='compchemlogfile',
                        nargs='+',
                        help=('one or more computational chemistry output '
                              'files to parse and convert'))
    parser.add_argument('--identifier',
                        default='logfiles',
                        help=('name of sheet which will contain DataFrame, if '
                              'writing to an Excel file, or identifier for '
                              'the group in HDFStore, if writing a HDF file'))
    parser.add_argument('-f', '--force', action='store_true',
                        help=('overwrite output file in case it already exists'))
    args = parser.parse_args()
    if args.output is not None and not args.force and os.path.exists(args.output):
        parser.exit(
            1, f'failure: exiting to avoid overwriting existing file "{args.output}"\n'
        )

    process_logfiles(args.compchemlogfiles, args.output, args.identifier)


if __name__ == "__main__":
    main()
