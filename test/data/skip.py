# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Tools for skipping data tests in cclib."""

from typing import Callable

from cclib.parser import ccData

import pytest


def skipForParser(parser: str, msg: str):
    """Return a decorator that skips the test for specified parser."""

    def tstdecorator(testfunc: Callable) -> Callable[[], None]:
        # breakpoint()
        def tstwrapper(self, data: ccData) -> None:
            # breakpoint()
            if data.parsername == parser:
                pytest.skip(reason=f"{parser}: {msg}")
            else:
                testfunc(self, data)

        return tstwrapper

    return tstdecorator


def skipForLogfile(fragment: str, msg: str):
    """Return a decorator that skips the test for logfiles containing fragment."""

    def tstdecorator(testfunc: Callable) -> Callable[[], None]:
        # breakpoint()
        def tstwrapper(self, data: ccData) -> None:
            # breakpoint()
            if any(fragment in filename for filename in data.filenames):
                pytest.skip(reason=f"{fragment}: {msg}")
            else:
                testfunc(self, data)

        return tstwrapper

    return tstdecorator
