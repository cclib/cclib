# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for parser utils module."""

import cclib


class convertorTest:
    def test_basic(self) -> None:
        """Are some basic conversions correct?"""
        convertor = cclib.parser.utils.convertor
        assert round(abs(convertor(1.89, "bohr", "Angstrom") - 1.0), 3) == 0
        assert round(abs(convertor(0.529, "Angstrom", "bohr") - 1.0), 3) == 0
        assert round(abs(convertor(627.5, "kcal/mol", "hartree") - 1.0), 3) == 0

    def test_wavenumber_to_frequency(self) -> None:
        """Is a wavenumber converted to a frequency in Hz, not GHz?"""

        convertor = cclib.parser.utils.convertor

        # nu [Hz] = c * nu_tilde, with c = 299792458 m/s exactly (SI) and
        # 1 cm^-1 = 100 m^-1.
        assert convertor(1.0, "wavenumber", "Hz") == 299792458.0 * 100

        # Chaining the table's own hartree entry has to land on the CODATA 2010
        # hartree-hertz relationship, 6.579683920729e15 Hz.
        hartree_in_hz = convertor(convertor(1.0, "hartree", "wavenumber"), "wavenumber", "Hz")
        assert abs(hartree_in_hz / 6.579683920729e15 - 1.0) < 1e-10

    def test_pairs(self) -> None:
        """Do flipped conversions correspond to each other?"""

        convertor = cclib.parser.utils.convertor

        pairs_proportional = (
            ("Angstrom", "bohr"),
            ("wavenumber", "eV"),
            ("wavenumber", "kcal/mol"),
            ("eV", "kJ/mol"),
            ("coulomb", "e"),
        )
        pairs_inverse = (("nm", "wavenumber"),)

        for unit1, unit2 in pairs_proportional:
            conv1 = convertor(1.0, unit1, unit2)
            conv2 = 1.0 / convertor(1.0, unit2, unit1)
            assert round(abs((conv1 - conv2) / conv1), 7) == 0
        for unit1, unit2 in pairs_inverse:
            conv1 = convertor(1.0, unit1, unit2)
            conv2 = convertor(1.0, unit2, unit1)
            assert round(abs((conv1 - conv2) / conv1), 7) == 0


class PeriodicTableTest:
    def setup_method(self) -> None:
        self.pt = cclib.parser.utils.PeriodicTable()

    def test_elements(self) -> None:
        """Does the periodic table give correct elements?"""
        assert self.pt.element[6] == "C"
        assert self.pt.element[44] == "Ru"
        assert self.pt.element[0] is None

    def test_numbers(self) -> None:
        """Does the periodic table give correct atom numbers?"""
        assert self.pt.number["C"] == 6
        assert self.pt.number["Au"] == 79
