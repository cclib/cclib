import sys


version_major = sys.version_info.major
version_minor = sys.version_info.minor


paths_allver = [
    'src/cclib/progress/qt4progress.py',
]

paths_2_7_only = [
    'src/cclib/bridge/cclib2pyquante.py',
]


def match_path(path, partial_paths):
    return any(partial_path in str(path)
               for partial_path in partial_paths)


def pytest_ignore_collect(path, config):
    if match_path(path, paths_allver):
        return True
    if version_major != 2:
        if match_path(path, paths_2_7_only):
            return True
    return False
