# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

from pathlib import Path

from cclib import ccDriver
from cclib.tree import Tree


def test_code_in_code() -> None:
    my_tree = Tree()
    my_tree.add_root()
    my_tree.add_child(0)

    # TODO replace with importlib
    __dir__ = Path(__file__).resolve().parent
    driver = ccDriver(
        str(__dir__ / ".." / "data" / "NBO" / "basicNBO7.0" / "basicORCA5.0" / "dvb_sp.out"),
        tree=my_tree,
    )
    data = driver.process_combinator()

    print("ORCA atomcharges")
    for i, j in data.parsed_data[0].atomcharges.items():
        print("\t", i, j)
    print("NBO atomcharges")
    for i, j in data.parsed_data[1].atomcharges.items():
        print("\t", i, j)
