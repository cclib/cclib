# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

from pathlib import Path

from cclib import ccDriver
from cclib.tree import Tree


def test_tree_psi4_bsse() -> None:
    my_tree = Tree()
    my_tree.add_root()  # supramolecular calculation
    my_tree.add_child(0)  # monomer1 calculation
    my_tree.add_child(0)  # monomer2 calculation

    # TODO replace with importlib
    __dir__ = Path(__file__).resolve().parent
    driver = ccDriver(
        str(__dir__ / ".." / "data" / "Psi4" / "basicPsi4-1.9" / "mp2_water_dimer_bsse.out"),
        tree=my_tree,
    )
    _data = driver.process_combinator()

    print(f"Supramolecular energy {driver._ccCollection._parsed_data[0].scfenergies[-1]}")
    for i in range(1, 3):
        print(f"Monomer {i:} energy {driver._ccCollection._parsed_data[i].scfenergies[-1]}")
