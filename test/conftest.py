"""test/conftest.py: dynamic testing configuration for pytest

See the pytest documentation for more details:
https://docs.pytest.org/en/latest/contents.html
"""

import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Mapping, Optional, Tuple, Union

from cclib.io import ccopen
from cclib.parser.data import ccData
from cclib.parser.logfileparser import Logfile

import pytest
import yaml


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
    """Representation of a regression."""

    # The fully-resolved location of the output file or directory of files to parse.
    loc_entry: Path
    all_files: Tuple[Path, ...]
    # The name of the outputfile transformed by `normalisefilename` for use as
    # part of a function (not class) in regression.py for testing, if one exists.
    normalisedfilename: str
    # Names of test classes (not functions) present in `regression.py`, either
    # defined there or imported (TODO), that should be parameterized with a
    # parsed `ccData` object, just like in the main unit tests.  Each test
    # class must exist: no automatic generation is done.
    tests: Optional[Tuple[str, ...]]
    # Should this regression entry be parsed?
    parse: bool
    # Testing via implicit discovery of functions and explicit listing of
    # classes is independent: the requested test classes given above may be
    # present but a special test function with the matching normalized name
    # may not be, and vice-versa.
    #
    # In the event that neither are present, TODO


class RegressionItem(pytest.Item):
    def __init__(self, *, regression: Regression, **kwargs) -> None:
        super().__init__(**kwargs)
        self.regression = regression

    def runtest(self) -> None:
        lfile = ccopen([str(x) for x in self.regression.all_files], future=True)
        if self.regression.parse:
            data = lfile.parse()
            assert self.regression not in _REGCACHE
            _REGCACHE[self.regression] = data


class RegressionFile(pytest.File):
    def collect(self) -> Iterator[RegressionItem]:
        rootdir = self.config.rootpath
        regression_dir = rootdir / "data" / "regression"
        regressions = read_regressionfiles_txt(regression_dir)
        for regression in regressions:
            yield RegressionItem.from_parent(
                self, name=f"parse{regression.normalisedfilename}", regression=regression
            )


def read_regressionfiles_txt(regression_dir: Path) -> List[Regression]:
    """Create a Regression for every entry in regressionfiles.txt."""
    entries = list()
    regfile = regression_dir / "regressionfiles_withtests.txt"
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
                tests = tuple(tokens[1:])
            loc_full = regression_dir / filename
            assert loc_full.exists()
            if loc_full.is_dir():
                all_files = tuple(sorted(loc_full.iterdir()))
            else:
                all_files = (loc_full,)
            entries.append(
                Regression(
                    loc_entry=Path(filename),
                    all_files=all_files,
                    normalisedfilename=normed,
                    tests=tests,
                    parse=True,
                )
            )
    return entries


def make_regression_entries() -> List[Regression]:
    """Create a Regression for every entry in regressionfiles.yaml."""
    __filedir__ = Path(__file__).resolve().parent
    __regression_dir__ = (__filedir__ / ".." / "data" / "regression").resolve()
    return read_regressionfiles_txt(__regression_dir__)


@pytest.fixture(scope="session")
def regression_entries() -> Dict[str, Regression]:
    return {entry.normalisedfilename: entry for entry in make_regression_entries()}


@pytest.fixture
def filename(request, regression_entries: Mapping[str, Regression]) -> Path:
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
    if normalized_name in regression_entries:
        return regression_entries[normalized_name].loc_entry
    # Allow explicitly skipped tests through.
    if "__unittest_skip__" in request.node.keywords:
        return None  # type: ignore
    raise RuntimeError(f"file not found for {normalized_name}")


def get_parsed_logfile(
    regression_entries: Mapping[str, Regression], normalized_name: str
) -> Logfile:
    """For a normalized filename suitable for a test function name and a
    mapping of these names to the absolute locations on the filesystem of
    their test files, parse the test file and return its data on the logfile
    instance.
    """
    # TODO List[Path] not allowed yet by ccopen?
    lfile = ccopen([str(fn) for fn in regression_entries[normalized_name].all_files], future=True)
    # TODO switch to cache
    lfile.data = lfile.parse()
    return lfile


@pytest.fixture
def logfile(request, regression_entries: Mapping[str, Regression]) -> Logfile:
    """For a test function whose name corresponds to a normalized filename,
    parse the corresponding data and return the logfile with data attached.
    """
    prefix = "test"
    assert request.node.name[: len(prefix)] == prefix
    normalized_name = request.node.name[len(prefix) :]
    if normalized_name in regression_entries:
        return get_parsed_logfile(regression_entries, normalized_name)
    # Workaround (?) for locations that are full directories (e.g. Turbomole)
    if normalized_name.endswith("__") and normalized_name[:-2] in regression_entries:
        return get_parsed_logfile(regression_entries, normalized_name[:-2])
    # Allow explicitly skipped tests through.
    if "__unittest_skip__" in request.node.keywords:
        return None  # type: ignore
    raise RuntimeError(f"file not found for {normalized_name}")


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
_REGCACHE: Dict[Regression, ccData] = {}


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


def pytest_collect_file(file_path: Path, parent) -> Optional[pytest.Collector]:
    # Ensure that each file specified in regressionfiles.txt is parsed, even
    # if a function- or class-based test doesn't exist for it.
    if file_path.name.startswith("regression"):
        # Only collect (and therefore run) if explicitly asked for
        # regression*.py.
        config = parent.config
        params = config.invocation_params
        if any(arg for arg in params.args if arg in str(file_path)):
            return RegressionFile.from_parent(parent, path=file_path)


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
        data_basedir = metafunc.config.rootpath / "data"
        for mtch in matches:
            subdir = data_basedir / get_program_dir(mtch["parser"]) / mtch["subdir"]
            files = [subdir / fn for fn in mtch["files"]]
            filegroups.append(files)
            ids.append(str(files[0].relative_to(data_basedir)))
        metafunc.parametrize(argnames="data", argvalues=filegroups, ids=ids, indirect=True)
    elif module_components[:2] == ["test", "regression"] and metafunc.cls is not None:
        # TODO rewrite to use Regressions
        target_class = metafunc.cls.__name__
        matches = [entry for entry in _REGRESSION_CLS_ENTRIES if target_class in entry.tests]
        filegroups = []
        ids = []
        regression_basedir = metafunc.config.rootpath / "data" / "regression"
        for mtch in matches:
            loc_entry = regression_basedir / mtch.loc_entry
            if not loc_entry.exists():
                raise RuntimeError(f"{loc_entry} doesn't exist")
            if loc_entry.is_dir():
                files = sorted(loc_entry.glob("*"))
            else:
                files = [loc_entry]
            filegroups.append(files)
            ids.append(str(mtch.loc_entry))
        metafunc.parametrize(argnames="data", argvalues=filegroups, ids=ids, indirect=True)
