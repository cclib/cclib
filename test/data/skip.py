# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

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