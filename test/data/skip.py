# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Tools for skipping data tests in cclib."""

from typing import Callable

from cclib.parser import ccData


def skipForParser(parser: str, msg: str):
    """Return a decorator that skips the test for specified parser."""

    def tstdecorator(testfunc: Callable):
        # def tstdecorator(*args, **kwargs):
        print(f"parser: {parser} msg: {msg}")

        # def tstwrapper(self, *args, **kwargs):
        # breakpoint()
        def tstwrapper(data: ccData) -> None:
            print(data)
            breakpoint()
            # if self.logfile.logname == parser:
            #     self.skipTest(msg)
            # else:
            #     testfunc(self, *args, **kwargs)
            pass

        return tstwrapper

    return tstdecorator


def skipForLogfile(fragment: str, msg: str):
    """Return a decorator that skips the test for logfiles containing fragment."""

    def tstdecorator(testfunc: Callable):
        # def tstdecorator(*args, **kwargs):
        # def tstwrapper(self, *args, **kwargs):
        print(f"fragment: {fragment} msg: {msg}")

        # breakpoint()
        def tstwrapper(data: ccData) -> None:
            print(data)
            breakpoint()
            # self.logfile.filename may be a string or list of strings.
            # if fragment in self.logfile.filename or any(fragment in filename for filename in self.logfile.filename):
            #     self.skipTest(msg)
            # else:
            #     testfunc(self, *args, **kwargs)
            pass

        return tstwrapper

    return tstdecorator
