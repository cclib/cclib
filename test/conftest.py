"""test/conftest.py: dynamic testing configuration for pytest

See the pytest documentation for more details:
https://docs.pytest.org/en/latest/contents.html
"""

import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Union

from cclib.io import ccopen
from cclib.parser.data import ccData
from cclib.parser.logfileparser import Logfile

import pytest


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
        if (x >= "a" and x <= "z") or (x >= "0" and x <= "9"):
            ans.append(y)
        else:
            ans.append("_")
    return "".join(ans)


@dataclass(frozen=True)
class Regression:
    filename: Path
    normalisedfilename: str
    tests: Optional[List[str]]


def read_regressionfiles_txt(regression_dir: Path) -> List[Regression]:
    """Create a Regression for every entry in regressionfiles.txt."""
    entries = list()
    regfile = regression_dir / "regressionfiles.txt"
    if regfile.is_file():
        contents = [line.split() for line in regfile.read_text(encoding="utf-8").splitlines()]
        for tokens in contents:
            assert len(tokens) >= 2
            filename = os.sep.join(tokens[0].split("/"))
            normed = normalisefilename(filename)
            # TODO error checking once format is more formalized
            if tokens[1] == "None":
                tests = None
            else:
                tests = tokens[1:]
            entries.append(
                Regression(
                    filename=regression_dir / filename, normalisedfilename=normed, tests=tests
                )
            )
    return entries


def make_regression_entries() -> List[Regression]:
    """Create a Regression for every entry in regressionfiles.txt."""
    __filedir__ = Path(__file__).resolve().parent
    __regression_dir__ = (__filedir__ / ".." / "data" / "regression").resolve()
    return read_regressionfiles_txt(__regression_dir__)


@pytest.fixture(scope="session")
def regression_filenames() -> Dict[str, Path]:
    """Map normalized filenames suitable for test function names to their
    absolute location on the filesystem.
    """
    return {entry.normalisedfilename: entry.filename for entry in make_regression_entries()}


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
    lfile = ccopen(fn, future=True)
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
_REGRESSIONDATA = make_regression_entries()
_REGRESSION_CLS_ENTRIES = [entry for entry in _REGRESSIONDATA if entry.tests is not None]


def get_program_dir(parser_name: str) -> str:
    """Return a directory name given a parser name.

    In at least one case (GAMESS-UK) the directory is named differently.
    """

    if parser_name == "GAMESSUK":
        return "GAMESS-UK"
    return parser_name


_CACHE: Dict[str, ccData] = {}


@pytest.fixture(scope="session")
def data(request) -> ccData:
    files = request.param
    first = files[0]
    if first not in _CACHE:
        logfile = ccopen([str(f) for f in files], future=True)
        data = logfile.parse()
        filenames = logfile.filename
        if not isinstance(filenames, list):
            filenames = [filenames]
        data.filenames = filenames
        data.parsername = logfile.logname
        _CACHE[first] = data
    return _CACHE[first]


def pytest_collection_modifyitems(session, config, items) -> None:
    """
    Official docstring:  Called after collection has been performed. May filter or re-order the items in-place.

    For us, that means:

      - All tests get marked with whatever their containing subdirectory is
        (is_parser, is_data, is_method, ...).  This is because `-k` does not
        operate on filenames, only function/class/method names, and many of
        these names do not contain the general category of feature we want to
        group by during testing.
    """
    rootdir = config.rootpath
    test_basedir = rootdir / "test"
    for item in items:
        try:
            test_subpath = item.path.relative_to(test_basedir)
            mark_name = f"is_{str(test_subpath.parent)}"
            item.add_marker(mark_name)
        except ValueError:
            # Not being a subpath is ok, that just means it isn't a test file.
            pass


def pytest_configure(config) -> None:
    """Automatically add a marker named for each test subdirectory, prefixed with `is_`."""
    rootdir = config.rootpath
    test_basedir = rootdir / "test"
    for path in test_basedir.glob("*"):
        if path.is_dir() and "__pycache__" not in str(path):
            mark_name = f"is_{path.stem}"
            config.addinivalue_line("markers", mark_name)


def pytest_generate_tests(metafunc: "pytest.Metafunc") -> None:
    module_components = metafunc.module.__name__.split(".")
    if module_components[:2] == ["test", "data"]:
        target_class = metafunc.cls.__name__
        matches = [entry for entry in _TESTDATA if entry["class"] == target_class]
        filegroups = []
        ids = []
        basedir = (Path(__file__) / ".." / "..").resolve()
        datadir = basedir / "data"
        for mtch in matches:
            subdir = datadir / get_program_dir(mtch["parser"]) / mtch["subdir"]
            files = [subdir / fn for fn in mtch["files"]]
            filegroups.append(files)
            ids.append(str(files[0].relative_to(datadir)))
        metafunc.parametrize(argnames="data", argvalues=filegroups, ids=ids, indirect=True)
    elif module_components[:2] == ["test", "regression"] and metafunc.cls is not None:
        target_class = metafunc.cls.__name__
        matches = [entry for entry in _REGRESSION_CLS_ENTRIES if target_class in entry.tests]
        filegroups = []
        ids = []
        for mtch in matches:
            if not mtch.filename.exists():
                raise RuntimeError(f"{mtch.filename} doesn't exist")
            if mtch.filename.is_dir():
                files = sorted(mtch.filename.glob("*"))
            else:
                files = [mtch.filename]
            filegroups.append(files)
            # TODO factor out of loop
            regressiondir = (mtch.filename / ".." / ".." / "..").resolve()
            ids.append(str(mtch.filename.relative_to(regressiondir)))
        metafunc.parametrize(argnames="data", argvalues=filegroups, ids=ids, indirect=True)
