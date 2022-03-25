# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import sys


class TextProgress:

    def __init__(self):

        self.nstep = 0
        self.text = None
        self.oldprogress = 0
        self.progress = 0
        self.calls = 0

    def initialize(self, nstep, text=None):

        self.nstep = float(nstep)
        self.text = text

        #sys.stdout.write("\n")

    def update(self, step, text=None):

        self.progress = int(step * 100 / self.nstep)

        if self.progress/2 >= self.oldprogress/2 + 1 or self.text != text:
        # just went through at least an interval of ten, ie. from 39 to 41,
        # so update

            mystr = "\r["
            prog = int(self.progress / 10)
            mystr += prog * "=" + (10-prog) * "-"
            mystr += f"] {int(self.progress):3}%"

            if text:
                mystr += f"    {text}"

            sys.stdout.write(f"\r{70 * ' '}")
            sys.stdout.flush()
            sys.stdout.write(mystr)
            sys.stdout.flush()
            self.oldprogress = self.progress

            if self.progress >= 100 and text == "Done":
                print(" ")

        return
