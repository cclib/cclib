# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from cclib.parser.gaussianparser import parse_route_line, parse_route_lines


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
