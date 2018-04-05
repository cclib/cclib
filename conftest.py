import sys


version_major = sys.version_info.major
version_minor = sys.version_info.minor

# Paths that should be ignored for all Python versions.
paths_allver = [
    'src/cclib/progress/qt4progress.py',
]

# Paths that should run only for Python 2.7.
paths_only_2_7 = [
    'src/cclib/bridge/cclib2pyquante.py',
]

# Paths that should run everywhere except Python 3.2.
paths_not_3_2 = [
    'src/cclib/bridge/cclib2biopython.py',
]


def match_path(path, partial_paths):
    """Does the given path contain any of the stubs in partial_paths?"""
    return any(partial_path in str(path)
               for partial_path in partial_paths)


def pytest_ignore_collect(path, config):
    """If this returns True for a given path, pytest will ignore it."""
    if match_path(path, paths_allver):
        return True
    if version_major != 2:
        if match_path(path, paths_only_2_7):
            return True
    if version_major == 3 and version_minor == 2:
        if match_path(path, paths_not_3_2):
            return True
    return False
