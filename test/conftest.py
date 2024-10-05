# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""test/conftest.py: dynamic testing configuration for pytest

See the pytest documentation for more details:
https://docs.pytest.org/en/latest/contents.html
"""

import os
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Mapping, Optional, Tuple, Union

from cclib.attribute_parsers.data import ccData
from cclib.file_handler import FileHandler
from cclib.io import ccopen

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
    """Representation of a regression.

    A regression, at a minimum, consists of one or more outputs
    from a run of a computational chemistry program.
    These outputs are either
    - old unit tests or
    - one-offs that failed to parse or did parse but presented incorrect results.

    There are three kinds of testing that are currently performed:
    - Parsing the file without erroring
    - Testing the output with test classes, used for old unit tests
    - Testing the output with a standalone function, used for one-offs

    All regressions must have an entry in `regressionfiles.yaml`,
    found in the cclib-data repository, where each dictionary under
    the top-level `regressions` key must contain `loc_entry` and
    optionally contain `tests` and `parse`, all described below.

    All the kinds of testing are independent of each other;
    in particular, using the test classes or functions does not require parsing.
    In order to parse a regression, `parse` must be true.
    In order to test using classes, the name of each class must be added to `tests`.
    In order to test using a function, a function must exist whose name starts
    with either 'test' or 'testnoparse' followed by the contents of `normalisedfilename`.
    """

    # The non-fully-resolved location of the output file or directory of files to parse.
    loc_entry: Path
    # The fully-resolved location of all files to parse.
    # For example, Turbomole output is spread over multiple files.
    all_files: Tuple[Path, ...]
    # The name of the output file transformed by `normalisefilename` for use as
    # part of a function (not class) in regression.py for testing,
    # if such a function exists.
    normalisedfilename: str
    # Names of test classes (not functions) present in `regression.py`,
    # either defined there or imported from the unit test suite,
    # that should be parameterized with a parsed `ccData` object,
    # just like in the main unit tests.
    # Each test class must exist: no automatic generation is done.
    tests: Optional[Tuple[str, ...]]
    # Should this regression entry be parsed?
    # - If true (default), the regression will be parsed even if
    #   no test classes were requested and no test function was defined.
    # - If false, the regression files are still collected into a Logfile
    #   but parsing is not performed.
    parse: bool


class RegressionItem(pytest.Item):
    """A pytest item, representing a single pytest test, that parses a regression."""

    def __init__(self, *, regression: Regression, **kwargs) -> None:
        super().__init__(**kwargs)
        self.regression = regression

    def runtest(self) -> None:
        if self.regression not in _REGCACHE:
            _REGCACHE[self.regression] = parse(self.regression)


class RegressionFile(pytest.File):
    """A pytest collector that ensures all requested regressions can be parsed
    when regression.py is collected, even if no tests are defined for that regression.
    """

    def collect(self) -> Iterator[RegressionItem]:
        rootdir = self.config.rootpath
        regression_dir = rootdir / "data" / "regression"
        regressions = read_regressionfiles_yaml(regression_dir)
        for regression in regressions:
            # It's not sufficient to check for this inside
            # RegressionItem.runtest(), because then the test will show as
            # passing even though it was never parsed.
            if regression.parse:
                yield RegressionItem.from_parent(
                    self, name=f"parse{regression.normalisedfilename}", regression=regression
                )


def read_regressionfiles_yaml(regression_dir: Path) -> List[Regression]:
    """Create a Regression for every entry in regressionfiles.yaml."""
    regressions = list()
    regfile = regression_dir / "regressionfiles.yaml"
    if regfile.is_file():
        entries = yaml.safe_load(regfile.read_text(encoding="utf-8"))
        assert set(entries.keys()) == {"regressions"}
        for entry in entries["regressions"]:
            loc_entry = os.sep.join(entry["loc_entry"].split("/"))
            tests = entry.get("tests", None)
            if tests is not None:
                tests = tuple(tests)
            loc_full = regression_dir / loc_entry
            assert loc_full.exists()
            if loc_full.is_dir():
                all_files = tuple(sorted(loc_full.iterdir()))
            else:
                all_files = (loc_full,)
            regressions.append(
                Regression(
                    loc_entry=Path(loc_entry),
                    all_files=all_files,
                    normalisedfilename=normalisefilename(loc_entry),
                    tests=tests,
                    parse=entry.get("parse", True),
                )
            )
    return regressions


def make_regression_entries() -> List[Regression]:
    """Create a Regression for every entry in regressionfiles.yaml."""
    __filedir__ = Path(__file__).resolve().parent
    __regression_dir__ = (__filedir__ / ".." / "data" / "regression").resolve()
    return read_regressionfiles_yaml(__regression_dir__)


@pytest.fixture(scope="session")
def _regression_entries() -> Dict[str, Regression]:
    """Not meant to be used in test code outside of conftest.py."""
    return {entry.normalisedfilename: entry for entry in make_regression_entries()}


@pytest.fixture
def filename(request, _regression_entries: Mapping[str, Regression]) -> Path:
    """For a test function whose name corresponds to a normalized filename,
    starting with 'testnoparse', get the absolute location on the filesystem
    of the corresponding test data.

    The only tests that can use this fixture are those marked as 'noparse',
    which typically instantiate the logfile object manually for manipulation.
    Most tests require a parse and should use the `logfile` fixture.
    """
    prefix = "testnoparse"
    assert request.node.name[: len(prefix)] == prefix
    normalized_name = request.node.name[len(prefix) :]
    if normalized_name in _regression_entries:
        return (
            request.config.rootpath
            / "data"
            / "regression"
            / _regression_entries[normalized_name].loc_entry
        )
    # Allow explicitly skipped tests through.
    if "__unittest_skip__" in request.node.keywords:
        return None  # type: ignore
    raise RuntimeError(f"file not found for {normalized_name}")


def parse(regression: Regression) -> FileHandler:
    """Collect the contents of a Regression into a Logfile,
    optionally parsing the logfile.
    """
    lfile = ccopen([str(x) for x in regression.all_files], future=True)
    if regression.parse:
        data = lfile.process_combinator()
        lfile.data = data
    return lfile


def get_parsed_logfile(
    regression_entries: Mapping[str, Regression], normalized_name: str
) -> FileHandler:
    """For a normalized filename suitable for a test function name and a
    mapping of these names to the absolute locations on the filesystem of
    their test files, parse the test file and return its data on the logfile
    instance.
    """
    assert normalized_name in regression_entries
    regression = regression_entries[normalized_name]
    if regression not in _REGCACHE:
        _REGCACHE[regression] = parse(regression)
    return _REGCACHE[regression]


@pytest.fixture
def logfile(request, _regression_entries: Mapping[str, Regression]) -> FileHandler:
    """For a test function whose name corresponds to a normalized filename,
    starting with 'test', parse the corresponding data and return the logfile
    with data attached.
    """
    prefix = "test"
    assert request.node.name[: len(prefix)] == prefix
    normalized_name = request.node.name[len(prefix) :]
    if normalized_name in _regression_entries:
        return get_parsed_logfile(_regression_entries, normalized_name)
    # Workaround (?) for locations that are full directories (e.g. Turbomole)
    if normalized_name.endswith("__") and normalized_name[:-2] in _regression_entries:
        return get_parsed_logfile(_regression_entries, normalized_name[:-2])
    # Allow explicitly skipped tests through.
    if "__unittest_skip__" in request.node.keywords:
        return None  # type: ignore
    raise RuntimeError(f"file not found for {normalized_name}")


def gettestdata() -> List[Dict[str, Union[str, List[str]]]]:
    """Return a dict of the unit test file data."""

    lines = (Path(__file__).resolve().parent / "testdata").read_text(encoding="utf-8").splitlines()

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
    # This is duplicated from test/test_data.py.
    if parser_name == "GAMESSUK":
        return "GAMESS-UK"
    return parser_name


# A naive caching system so that parsed unit and regression tests are not
# re-parsed.  This is necessary since not all code in conftest.py that
# requires parsed results can take advantage of session-scoped fixtures.
#
# TODO see if this can be replaced with
# https://docs.pytest.org/en/latest/how-to/writing_hook_functions.html#storing-data-on-items-across-hook-functions
_CACHE: Dict[str, ccData] = {}
# Each logfile, if Regression.parse == True, will have a `.data` member.
_REGCACHE: Dict[Regression, FileHandler] = {}


@pytest.fixture(scope="session")
def data(request) -> ccData:
    """Parse a unit test, placing the filenames and parser name
    on the data instance.
    """
    files = request.param
    first = files[0]
    if first not in _CACHE:
        ccdriver_instance = ccopen([str(f) for f in files], future=True)
        ccdriver_instance.process_combinator()
        filenames = ccdriver_instance._fileHandler.filenames
        if not isinstance(filenames, list):
            filenames = [filenames]
        ccdriver_instance._fileHandler.filenames = filenames
        ccdriver_instance.parsername = ccdriver_instance.identified_program
        # cc.parserclassname = str(ccdriver_instance.__class__).split(".")[-1][:-2]
        _CACHE[first] = ccdriver_instance
    return _CACHE[first]


def pytest_collect_file(file_path: Path, parent) -> Optional[pytest.Collector]:
    """Ensure that each regression specified in regressionfiles.yaml is parsed,
    even if a function- or class-based test doesn't exist for it,
    as long as the regression isn't explicitly marked to not parse.
    """
    if file_path.name == "regression.py":
        # Only collect (and therefore run) if explicitly asked for
        # regression.py.
        #
        # If regression_io.py or future files ever want to separate parsing
        # from testing, this condition will need to be changed.
        config = parent.config
        params = config.invocation_params
        if any(arg for arg in params.args if arg in str(file_path)):
            return RegressionFile.from_parent(parent, path=file_path)


def pytest_collection_modifyitems(config, items) -> None:
    """
    Official docstring:  Called after collection has been performed.
    May filter or re-order the items in-place.

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
            parent = str(test_subpath.parent)
            if parent != ".":
                mark_name = f"is_{parent}"
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
    """For class-based tests that accept the data fixture, ensure they are parameterized
    over all requested unit and regression test entries.
    """
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


_EXCLUDE = {"filenames", "parserclassname", "parsername"}


def pytest_sessionfinish(session: "pytest.Session") -> None:
    """Write out coverage information used for building docs.

    The coverage data is a dictionary that maps each parser name to
    all the attribute names created across that parser's unit tests.
    Strictly speaking, it is the parser's class name, not the 'logname'.

    We place coverage collection here rather than in a pytest reporting hook since
    this seems to be the only relevant global hook that runs at the end of a test session.
    """
    if _CACHE:
        coverage_accumulate = defaultdict(set)  # noqa: F841
        # for data in _CACHE.values():
        #    coverage_accumulate[data.__dict__["parserclassname"]].update(data.__dict__.keys())
        # for parserclassname in coverage_accumulate:
        #   coverage_accumulate[parserclassname] -= _EXCLUDE
        # coverage = {
        #    parserclassname: sorted(attributenames)
        #    for parserclassname, attributenames in coverage_accumulate.items()
        # }
        # coverage_dir = pytest.Cache.cache_dir_from_config(session.config)
        # (coverage_dir / "coverage_unit.json").write_text(json.dumps(coverage), encoding="utf-8")
