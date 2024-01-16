# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Generate the coverage.rst and coverage.rst files from test results."""

import argparse
import itertools
import json
import os
import sys
from pathlib import Path

# Import cclib and check we are using the version from a subdirectory.
import cclib


def generate_coverage(cclib_base_dir: Path) -> str:
    """Generate a string containing a reStructuredTest table
    representation of which parsers support which attributes, based on
    test results.
    """
    lines = []

    os.chdir(cclib_base_dir)

    thispath = Path(__file__).parent
    # TODO this isn't correct but doesn't seem to hurt things right now.
    sys.path.insert(1, os.path.join(thispath, cclib_base_dir))

    logpath = thispath / "coverage.tests.log"

    from test.test_data import all_parsers, parser_names

    try:
        # unittest-based testing
        import inspect
        from test.test_data import DataSuite, all_modules

        ds_args = inspect.getfullargspec(DataSuite.__init__).args
        try:
            with open(logpath, "w") as flog:
                stdout_backup = sys.stdout
                sys.stdout = flog
                alltests = {}
                for p in parser_names:
                    assert "parsers" in ds_args
                    suite = DataSuite(parsers={p: all_parsers[p]}, modules=all_modules, stream=flog)
                    suite.testall()
                    alltests[p] = [{"data": t.data} for t in suite.alltests]
                sys.stdout = stdout_backup
        except Exception as e:
            print("Unit tests did not run correctly. Check log file for errors:")
            with open(logpath) as fh:
                print(fh.read())
            print(e)
            sys.exit(1)
        coverage = {
            parser_name: set(
                itertools.chain.from_iterable(
                    [datadict["data"].__dict__.keys() for datadict in alltests[parser_name]]
                )
            )
            for parser_name in alltests
        }
    except ImportError:
        # pytest-based testing
        import pytest

        class CaptureCoverageDir:
            def pytest_sessionfinish(self, session: "pytest.Session") -> None:
                # See conftest.py.
                self.coverage_dir = pytest.Cache.cache_dir_from_config(session.config)

        capture_coverage_dir = CaptureCoverageDir()
        with open(logpath, "w") as flog:
            stdout_backup = sys.stdout
            sys.stdout = flog
            # Ignore one or more checked-out cclib source trees under doc/sphinx/,
            # since there will be conftest.py multiple detection issues.
            # TODO specify in pytest config?
            retcode = pytest.main(["--ignore=doc", "-m", "is_data"], plugins=[capture_coverage_dir])
            sys.stdout = stdout_backup
        if retcode != 0:
            print("Unit tests did not run correctly. Check log file for errors:")
            with open(logpath) as fh:
                print(fh.read())
            sys.exit(1)
        coverage = json.loads(
            (capture_coverage_dir.coverage_dir / "coverage_unit.json").read_text(encoding="utf-8")
        )

    attributes = sorted(cclib.parser.data.ccData._attrlist)

    ncols = len(parser_names) + 1
    colwidth = 4 + max(len(attribute) for attribute in attributes)
    colfmt = f"%-{int(colwidth)}s"
    dashes = f"{'=' * (colwidth - 1)} " * ncols

    lines.append(dashes)
    lines.append(colfmt * ncols % tuple(["attributes"] + parser_names))
    lines.append(dashes)

    # Eventually we want to move this to cclib, too.
    not_applicable = {
        "ADF": ["aonames", "ccenergies", "mpenergies"],
        "DALTON": ["fonames", "fooverlaps", "fragnames", "frags"],
        "GAMESS": ["fonames", "fooverlaps", "fragnames", "frags"],
        "GAMESSUK": ["fonames", "fooverlaps", "fragnames", "frags"],
        "Gaussian": ["fonames", "fooverlaps", "fragnames", "frags"],
        "Jaguar": ["fonames", "fooverlaps", "fragnames", "frags"],
        "Molpro": ["fonames", "fooverlaps", "fragnames", "frags"],
        "NWChem": ["fonames", "fooverlaps", "fragnames", "frags"],
        "ORCA": ["fonames", "fooverlaps", "fragnames", "frags"],
        "Psi4": ["fonames", "fooverlaps", "fragnames", "frags"],
        "QChem": ["fonames", "fooverlaps", "fragnames", "frags"],
    }
    not_possible = {
        "NWChem": ["vibfconsts", "vibrmasses"],
        "Psi4": ["aooverlaps", "vibirs"],
        "QChem": ["aooverlaps", "etrotats"],
    }

    # For each attribute, get a list of Boolean values for each parser that flags
    # if it has been parsed by at least one unit test. Substitute an OK sign or
    # T/D appropriately, with the exception of attributes that have been explicitely
    # designated as N/A.
    for attr in attributes:
        parsed = [attr in coverage[p] for p in parser_names]
        for ip, p in enumerate(parsed):
            if p:
                parsed[ip] = "√"
            else:
                if attr in not_applicable.get(parser_names[ip], []):
                    parsed[ip] = "N/A"
                elif attr in not_possible.get(parser_names[ip], []):
                    parsed[ip] = "N/P"
                else:
                    parsed[ip] = "T/D"
        lines.append(colfmt * ncols % tuple([f"`{attr}`_"] + parsed))

    lines.append(dashes)
    lines.append("")

    for attr in attributes:
        lines.append(f".. _`{attr}`: data_notes.html#{attr}")

    return "\n".join(lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("cclib_base_dir", type=Path)
    args = parser.parse_args()
    print(generate_coverage(args.cclib_base_dir))
