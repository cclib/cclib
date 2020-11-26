#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import logging
from argparse import ArgumentParser

from cclib.io import ccread
from cclib.method import CDA


def main():
    parser = ArgumentParser()
    parser.add_argument("file1", help="logfile containing the supermolecule")
    parser.add_argument("file2", help="logfile containing the first fragment")
    parser.add_argument("file3", help="logfile containing the second fragment")
    args = parser.parse_args()

    loglevel = logging.ERROR

    data1 = ccread(args.file1, loglevel=loglevel)
    data2 = ccread(args.file2, loglevel=loglevel)
    data3 = ccread(args.file3, loglevel=loglevel)

    fa = CDA(data1, None, loglevel)
    retval = fa.calculate([data2, data3])

    if retval:

        print("Charge decomposition analysis of {}\n".format(args.file1))

        if len(data1.homos) == 2:
            print("ALPHA SPIN:")
            print("===========")

        print(" MO#      d       b       r       s")
        print("-------------------------------------")

        for spin in range(len(data1.homos)):

            if spin == 1:
                print("\nBETA SPIN:")
                print("==========")

            for i in range(len(fa.donations[spin])):

                print("%4i: %7.3f %7.3f %7.3f %7.3f" % \
                        (i + 1, fa.donations[spin][i],
                                fa.bdonations[spin][i],
                                fa.repulsions[spin][i],
                                fa.residuals[spin][i]))

                if i == data1.homos[spin]:
                    print("------ HOMO - LUMO gap ------")
                    

            print("-------------------------------------")
            print(" T:   %7.3f %7.3f %7.3f %7.3f" % \
                    (fa.donations[spin].sum(),
                        fa.bdonations[spin].sum(),
                        fa.repulsions[spin].sum(),
                        fa.residuals[spin].sum()))


if __name__ == '__main__':
    main()
