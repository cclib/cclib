"""test/conftest.py: dynamic testing configuration for pytest

See the pytest documentation for more details:
https://docs.pytest.org/en/latest/contents.html
"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Mapping, Union

from cclib.io import ccopen
from cclib.parser.data import ccData
from cclib.parser.logfileparser import Logfile

import pytest

version_major = sys.version_info.major


def normalisefilename(filename: str) -> str:
    """Replace all non-alphanumeric symbols by underscores.

    >>> from . import regression
    >>> for x in [ "Gaussian/Gaussian03/Mo4OSibdt2-opt.log" ]:
    ...     print(regression.normalisefilename(x))
    ...
    Gaussian_Gaussian03_Mo4OSibdt2_opt_log
    """
    # TODO This is duplicated from test/regression.py, where it can be removed
    # after all (dynamically generated) regression tests are converted to run
    # via pytest.
    ans = []
    for y in filename:
        x = y.lower()
        if (x >= "a" and x <= "z") or (x >= "0" and x <= "9"):
            ans.append(y)
        else:
            ans.append("_")
    return "".join(ans)


@pytest.fixture(scope="session")
def regression_filenames() -> Dict[str, Path]:
    """Map normalized filenames suitable for test function names to their
    absolute location on the filesystem.
    """
    __filedir__ = Path(__file__).resolve().parent
    __regression_dir__ = (__filedir__ / ".." / "data" / "regression").resolve()
    regfile = __regression_dir__ / "regressionfiles.txt"
    regfilenames = [
        os.sep.join(x.strip().split("/")) for x in regfile.read_text(encoding="utf-8").splitlines()
    ]
    return {normalisefilename(filename): __regression_dir__ / filename for filename in regfilenames}


@pytest.fixture
def filename(request, regression_filenames: Mapping[str, Path]) -> Path:
    """For a test function whose name corresponds to a normalized filename,
    get the absolution location on the filesystem of the corresponding test
    data.

    The only tests that can use this fixture are those marked as 'noparse',
    which typically instantiate the logfile object manually for manipulation.
    Most tests require a parse and should use the logfile fixture.
    """
    prefix = "testnoparse"
    assert request.node.name[: len(prefix)] == prefix
    normalized_name = request.node.name[len(prefix) :]
    if normalized_name in regression_filenames:
        return regression_filenames[normalized_name]
    # Allow explicitly skipped tests through.
    if "__unittest_skip__" in request.node.keywords:
        return None  # type: ignore
    raise RuntimeError


def get_parsed_logfile(regression_filenames: Mapping[str, Path], normalized_name: str) -> Logfile:
    """For a normalized filename suitable for a test function name and a
    mapping of these names to the absolute locations on the filesystem of
    their test files, parse the test file and return its data on the logfile
    instance.
    """
    fn = regression_filenames[normalized_name]
    if fn.is_dir():
        # TODO List[Path] not allowed yet by ccopen?
        fn = [str(x) for x in sorted(fn.iterdir())]
    lfile = ccopen(fn)
    lfile.data = lfile.parse()
    return lfile


@pytest.fixture
def logfile(request, regression_filenames: Mapping[str, Path]) -> Logfile:
    """For a test function whose name corresponds to a normalized filename,
    parse the corresponding data and return the logfile with data attached.
    """
    prefix = "test"
    assert request.node.name[: len(prefix)] == prefix
    normalized_name = request.node.name[len(prefix) :]
    if normalized_name in regression_filenames:
        return get_parsed_logfile(regression_filenames, normalized_name)
    # Workaround (?) for locations that are full directories (e.g. Turbomole)
    if normalized_name.endswith("__") and normalized_name[:-2] in regression_filenames:
        return get_parsed_logfile(regression_filenames, normalized_name[:-2])
    # Allow explicitly skipped tests through.
    if "__unittest_skip__" in request.node.keywords:
        return None  # type: ignore
    raise RuntimeError


def gettestdata() -> List[Dict[str, Union[str, List[str]]]]:
    """Return a dict of the test file data."""

    testdatadir = Path(__file__).resolve().parent
    with open(testdatadir / "testdata", encoding="utf-8") as testdatafile:
        lines = testdatafile.readlines()

    # Remove blank lines and those starting with '#'.
    lines = [line for line in lines if (line.strip() and line[0] != "#")]

    # Remove comment at end of lines (everything after a '#').
    lines = [line.split("#")[0] for line in lines]

    # Transform remaining lines into dictionaries.
    cols = [line.split() for line in lines]
    labels = ("module", "parser", "class", "subdir", "files")
    testdata = [dict(zip(labels, (c[0], c[1], c[2], c[3], c[4:]))) for c in cols]

    return testdata


_TESTDATA = gettestdata()


def get_program_dir(parser_name: str) -> str:
    """Return a directory name given a parser name.

    In at least one case (GAMESS-UK) the directory is named differently.
    """

    if parser_name == "GAMESSUK":
        return "GAMESS-UK"
    return parser_name


_CACHE = {}


@pytest.fixture(scope="session")
def data(request) -> ccData:
    files = request.param
    first = files[0]
    if first not in _CACHE:
        logfile = ccopen([str(f) for f in files])
        data = logfile.parse()
        data.filenames = logfile.filename
        data.parsername = logfile.logname
        _CACHE[first] = data
    return _CACHE[first]


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    module_components = metafunc.module.__name__.split(".")
    if module_components[:2] == ["test", "data"]:
        target_class = metafunc.cls.__name__
        matches = [entry for entry in _TESTDATA if entry["class"] == target_class]
        filegroups = []
        ids = []
        basedir = (Path(__file__) / ".." / "..").resolve()
        datadir = basedir / "data"
        for m in matches:
            subdir = datadir / get_program_dir(m["parser"]) / m["subdir"]
            files = [subdir / fn for fn in m["files"]]
            filegroups.append(files)
            ids.append(str(files[0].relative_to(datadir)))
        metafunc.parametrize(argnames="data", argvalues=filegroups, ids=ids, indirect=True)
