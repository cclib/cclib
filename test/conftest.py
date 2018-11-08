"""test/conftest.py: dynamic testing configuration for pytest

See the pytest documentation for more details:
https://docs.pytest.org/en/latest/contents.html
"""

import sys


version_major = sys.version_info.major

# Paths that should be ignored for all Python versions.
paths_ignore_allver = [
    'cclib/progress/qt4progress.py',
]

# Paths that should run only for Python 2.7.
paths_ignore_only_2_7 = [
    'cclib/bridge/cclib2pyquante.py',
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
    if version_major != 2:
        if match_path(path, paths_ignore_only_2_7):
            return True
    return False
