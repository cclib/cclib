# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
import logging
from pathlib import Path

from cclib.parser import Gaussian
from cclib.parser.gaussianparser import parse_route_line, parse_route_lines

import numpy


def test_parse_route_line() -> None:
    original = " 3/6=3,11=9,16=1,25=1,30=1/1,2,3;"
    ref = (3, {6: 3, 11: 9, 16: 1, 25: 1, 30: 1})
    assert parse_route_line(original) == ref


def test_parse_route_lines() -> None:
    original = [
        " 1/38=1/1;",
        " 2/12=2,17=6,18=5,40=1/2;",
        " 3/6=3,11=9,16=1,25=1,30=1/1,2,3;",
        " 4//1;",
        " 5/5=2,38=5/2;",
        " 6/7=2,8=2,9=2,10=2,28=1,79=1/1;",
        " 99/5=1,9=1/99;",
    ]
    ref = {
        1: {38: 1},
        2: {12: 2, 17: 6, 18: 5, 40: 1},
        3: {6: 3, 11: 9, 16: 1, 25: 1, 30: 1},
        4: {},
        5: {5: 2, 38: 5},
        6: {7: 2, 8: 2, 9: 2, 10: 2, 28: 1, 79: 1},
        99: {5: 1, 9: 1},
    }
    assert parse_route_lines(original) == ref


def test_rotconsts_asterisk(tmp_path: Path) -> None:
    """Rotational constants with asterisks (linear molecules) parse as NaN.

    Gaussian prints asterisks when a rotational constant overflows, e.g. for
    linear molecules (see #1958, also #1270). `utils.float` maps such values
    to NaN instead of raising ValueError.
    """
    logfile = tmp_path / "gaussian_asterisk_rotconsts.log"
    logfile.write_text(
        " Entering Gaussian System, Link 0=g16\n"
        " Rotational constants (GHZ): ************ 12.73690 12.73690\n"
        " Normal termination of Gaussian 16.\n"
    )
    data = Gaussian(str(logfile), loglevel=logging.ERROR).parse()
    assert data.rotconsts.shape == (1, 3)
    assert numpy.isnan(data.rotconsts[0][0])
    numpy.testing.assert_allclose(data.rotconsts[0][1:], [12.73690, 12.73690])
