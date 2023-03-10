"""test/conftest.py: dynamic testing configuration for pytest

See the pytest documentation for more details:
https://docs.pytest.org/en/latest/contents.html
"""

import logging
import os
import sys
from pathlib import Path
from typing import Dict, Mapping

import pytest

from test.test_data import (
    all_modules,
    all_parsers,
    module_names,
    parser_names,
)
from cclib.io import ccopen
from cclib.parser.logfileparser import Logfile

version_major = sys.version_info.major

# Paths that should be ignored for all Python versions.
paths_ignore_allver = [
    'cclib/progress/qt4progress.py',
]


def match_path(path, partial_paths):
    """Does the given path contain any of the stubs in partial_paths?"""
    return any(partial_path in str(path)
               for partial_path in partial_paths)


def pytest_ignore_collect(path, config):
    """pytest automatically runs this on every discovered path. If this
    returns True for a given path, pytest will ignore it.
    """
    if match_path(path, paths_ignore_allver):
        return True
    return False


def pytest_addoption(parser):
    parser.addoption("--terse", action="store_true")
    parser.addoption("--silent", action="store_true")


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == "test_all":
        metafunc.parametrize("parsers", [{p: all_parsers[p] for p in parser_names}])
        metafunc.parametrize("modules", [{p: all_modules[p] for p in module_names}])
        metafunc.parametrize("terse", [metafunc.config.getoption("--terse")])
        metafunc.parametrize("silent", [metafunc.config.getoption("--silent")])
        metafunc.parametrize("loglevel",
                             [logging.DEBUG if metafunc.config.getoption("--debug")
                              else logging.ERROR])
        metafunc.parametrize("summary", [True])
        metafunc.parametrize("visual_tests", [True])


def normalisefilename(filename: str) -> str:
    """Replace all non-alphanumeric symbols by underscores.

    >>> from . import regression
    >>> for x in [ "Gaussian/Gaussian03/Mo4OSibdt2-opt.log" ]:
    ...     print(regression.normalisefilename(x))
    ...
    Gaussian_Gaussian03_Mo4OSibdt2_opt_log
    """
    ans = []
    for y in filename:
        x = y.lower()
        if (x >= 'a' and x <= 'z') or (x >= '0' and x <= '9'):
            ans.append(y)
        else:
            ans.append("_")
    return "".join(ans)


@pytest.fixture(scope="session")
def filenames() -> Dict[str, Path]:
    __filedir__ = Path(__file__).resolve().parent
    __regression_dir__ = (__filedir__ / ".." / "data" / "regression").resolve()
    regfile = __regression_dir__ / "regressionfiles.txt"
    regfilenames = [
        os.sep.join(x.strip().split("/"))
        for x in regfile.read_text(encoding="utf-8").splitlines()
    ]
    return {
        normalisefilename(filename): __regression_dir__ / filename
        for filename in regfilenames
    }


@pytest.fixture
def filename(request, filenames: Mapping[str, Path]) -> Path:
    prefix = "testnoparse"
    assert request.node.name[:len(prefix)] == prefix
    normalized_name = request.node.name[len(prefix):]
    if normalized_name in filenames:
        return filenames[normalized_name]
    # Allow explicitly skipped tests through.
    if "__unittest_skip__" in request.node.keywords:
        return None  # type: ignore
    raise RuntimeError


@pytest.fixture
def logfile(request, filenames: Mapping[str, Path]) -> Logfile:
    prefix = "test"
    assert request.node.name[:len(prefix)] == prefix
    normalized_name = request.node.name[len(prefix):]
    if normalized_name in filenames:
        fn = filenames[normalized_name]
        if fn.is_dir():
            # FIXME List[Path] not allowed yet
            fn = [str(x) for x in sorted(fn.iterdir())]
        lfile = ccopen(fn)
        lfile.data = lfile.parse()
        return lfile
    if normalized_name.endswith("__") and normalized_name[:-2] in filenames:
        normalized_name = normalized_name[:-2]
        fn = filenames[normalized_name]
        if fn.is_dir():
            # FIXME List[Path] not allowed yet
            fn = [str(x) for x in sorted(fn.iterdir())]
        lfile = ccopen(fn)
        lfile.data = lfile.parse()
        return lfile
    # Allow explicitly skipped tests through.
    if "__unittest_skip__" in request.node.keywords:
        return None  # type: ignore
    raise RuntimeError
