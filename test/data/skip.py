# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import sys

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
            # self.logfile.filename may be a string or list of strings.
            if fragment in self.logfile.filename or any(fragment in filename for filename in self.logfile.filename):
                self.skipTest(msg)
            else:
                testfunc(self, *args, **kwargs)
        return testwrapper
    return testdecorator

def skipForVersion(version, msg):
    """Return a decorator that skips the test for specified python version."""
    def testdecorator(testfunc):
        def testwrapper(self, *args, **kwargs):
            curr_ver = sys.version_info
            skip_test = True
            if "major" in version:
                if version["major"] != curr_ver.major:
                    skip_test = False
            if "minor" in version:
                if version["minor"] != curr_ver.minor:
                    skip_test = False
            if "micro" in version:
                if version["micro"] != curr_ver.micro:
                    skip_test = False
            if skip_test:
                self.skipTest(msg)
            else:
                testfunc(self, *args, **kwargs)
        return testwrapper
    return testdecorator


