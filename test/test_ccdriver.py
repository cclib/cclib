# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import TYPE_CHECKING

from cclib.driver import ccDriver

if TYPE_CHECKING:
    from pathlib import Path


def test_ccdriver_identified_program(tmp_path: "Path") -> None:
    """Functionality for storing current program"""
    empty_source = tmp_path / "empty_source.txt"
    empty_source.touch()
    driver = ccDriver(empty_source)

    # basic usage
    assert driver.identified_program is None
    driver.identified_program = "Psi4"
    assert driver.identified_program == "Psi4"
    driver.identified_program = None
    assert driver.identified_program is None

    # two programs
    driver.identified_program = "hello"
    driver.identified_program = "world"
    assert driver.identified_program == "world"
    driver.identified_program = None
    assert driver.identified_program == "hello"
    driver.identified_program = None
    assert driver.identified_program is None
