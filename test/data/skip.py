# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Tools for skipping data tests in cclib."""


def skipForParser(parser, msg):
    """Return a decorator that skips the test for specified parser."""
    def testdecorator(testfunc):
        def testwrapper(self, *args, **kwargs):
            if self.logfile.logname == parser:
                self.skipTest(msg)
            else:
                testfunc(self, *args, **kwargs)
        return testwrapper
    return testdecorator


def skipForLogfile(fragment, msg):
    """Return a decorator that skips the test for logfiles containing fragment."""
    def testdecorator(testfunc):
        def testwrapper(self, *args, **kwargs):
            if fragment in self.logfile.filename:
                self.skipTest(msg)
            else:
                testfunc(self, *args, **kwargs)
        return testwrapper
    return testdecorator
