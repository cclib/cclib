#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Script for writing data tables from computational chemistry files."""

import sys
import os.path
import argparse

try:
    import pandas as pd
    _has_pandas = True
except ImportError:
    _has_pandas = False

from cclib.io import ccopen
from cclib.io import ccframe


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-O', '--outputdest',
                        help=('the output document to write, including an '
                              'extension supported by pandas '
                              '(csv, hdf/hdf5, json/cjson, pickle/pkl, xlsx)'))
    parser.add_argument('compchemlogfile',
                        nargs='+',
                        help=('one or more computational chemistry output '
                              'files to parse and convert'))
    parser.add_argument('--identifier',
                        default='logfiles',
                        help=('name of sheet which will contain DataFrame, if '
                              'writing to an Excel file, or identifier for '
                              'the group in HDFStore, if writing a HDF file'))
    args = parser.parse_args()

    outputdest = args.outputdest
    identifier = args.identifier
    filenames = args.compchemlogfile

    df = ccframe([ccopen(path) for path in filenames])

    if outputdest is not None:
        outputtype = os.path.splitext(os.path.basename(outputdest))[1][1:]

        if outputtype in {'csv'}:
            df.to_csv(outputdest)
        elif outputtype in {'hdf', 'hdf5'}:
            df.to_hdf(outputdest, key=identifier)
        elif outputtype in {'json', 'cjson'}:
            df.to_json(outputdest)
        elif outputtype in {'pickle', 'pkl'}:
            df.to_pickle(outputdest)
        elif outputtype in {'xlsx'}:
            writer = pd.ExcelWriter(outputdest)
            # This overwrites previous sheets
            # (see https://stackoverflow.com/a/42375263/4039050)
            df.to_excel(writer, sheet_name=identifier)
            writer.save()
    else:
        print(df)


if __name__ == "__main__":
    main()
