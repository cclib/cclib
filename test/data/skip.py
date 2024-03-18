# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Tools for skipping data tests in cclib."""

from inspect import signature
from typing import Callable

import pytest


def skipForParser(parser: str, msg: str):
    """Return a decorator that skips the test for specified parser."""

    def tstdecorator(testfunc: Callable) -> Callable[[], None]:
        func_args = list(signature(testfunc).parameters.keys())
        if "numvib" in func_args:

            def tstwrapper(self, data, numvib) -> None:
                if data.parsername == parser:
                    pytest.skip(reason=f"{parser}: {msg}")
                else:
                    testfunc(self, data, numvib)
        elif "extra" in func_args:

            def tstwrapper(self, data, extra) -> None:
                if data.parsername == parser:
                    pytest.skip(reason=f"{parser}: {msg}")
                else:
                    testfunc(self, data, extra)
        else:

            def tstwrapper(self, data) -> None:
                if data.parsername == parser:
                    pytest.skip(reason=f"{parser}: {msg}")
                else:
                    testfunc(self, data)

        return tstwrapper

    return tstdecorator


def skipForLogfile(fragment: str, msg: str):
    """Return a decorator that skips the test for logfiles containing fragment."""

    def tstdecorator(testfunc: Callable) -> Callable[[], None]:
        func_args = list(signature(testfunc).parameters.keys())
        if "numvib" in func_args:

            def tstwrapper(self, data, numvib) -> None:
                if any(fragment in filename for filename in data.filenames):
                    pytest.skip(reason=f"{fragment}: {msg}")
                else:
                    testfunc(self, data, numvib)
        elif "extra" in func_args:

            def tstwrapper(self, data, extra) -> None:
                if any(fragment in filename for filename in data.filenames):
                    pytest.skip(reason=f"{fragment}: {msg}")
                else:
                    testfunc(self, data, extra)
        else:

            def tstwrapper(self, data) -> None:
                if any(fragment in filename for filename in data.filenames):
                    pytest.skip(reason=f"{fragment}: {msg}")
                else:
                    testfunc(self, data)

        return tstwrapper

    return tstdecorator
