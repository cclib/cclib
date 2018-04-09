#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

from __future__ import print_function
import os
import sys
import glob
import getopt
import logging

import numpy

from cclib.io import ccopen
from cclib.method import CDA


def main():
    parser1 = ccopen(sys.argv[1], logging.ERROR)
    parser2 = ccopen(sys.argv[2], logging.ERROR)
    parser3 = ccopen(sys.argv[3], logging.ERROR)

    data1 = parser1.parse(); data2 = parser2.parse(); data3 = parser3.parse()
    fa = CDA(data1, None, logging.ERROR)
    retval = fa.calculate([data2, data3])

    if retval:

        print("Charge decomposition analysis of %s\n"%(sys.argv[1]))

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
