"""test/conftest.py: dynamic testing configuration for pytest

See the pytest documentation for more details:
https://docs.pytest.org/en/latest/contents.html
"""

import logging
import sys

from test.test_data import (
    all_modules,
    all_parsers,
    module_names,
    parser_names,
)

version_major = sys.version_info.major

# Paths that should be ignored for all Python versions.
paths_ignore_allver = [
    "cclib/progress/qt4progress.py",
]

# Paths that should run only for Python 2.7.
paths_ignore_only_2_7 = [
    "cclib/bridge/cclib2pyquante.py",
]


def match_path(path, partial_paths):
    """Does the given path contain any of the stubs in partial_paths?"""
    return any(partial_path in str(path) for partial_path in partial_paths)


def pytest_ignore_collect(path, config):
    """pytest automatically runs this on every discovered path. If this
    returns True for a given path, pytest will ignore it.
    """
    if match_path(path, paths_ignore_allver):
        return True
    if version_major != 2:
        if match_path(path, paths_ignore_only_2_7):
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
        metafunc.parametrize(
            "loglevel", [logging.DEBUG if metafunc.config.getoption("--debug") else logging.ERROR]
        )
        metafunc.parametrize("summary", [True])
        metafunc.parametrize("visual_tests", [True])
