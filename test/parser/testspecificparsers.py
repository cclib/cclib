# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for specific parser behaviors, such as overriden methods."""

import io


class NormalisesymTest:
    # Not needed: DALTON, MOPAC, NWChem, ORCA, QChem

    def test_normalisesym_adf(self) -> None:
        from cclib.parser.adfparser import ADF

        sym = ADF(io.StringIO("dummyfile")).normalisesym
        labels = [
            "A",
            "s",
            "A1",
            "A1.g",
            "Sigma",
            "Pi",
            "Delta",
            "Phi",
            "Sigma.g",
            "A.g",
            "AA",
            "AAA",
            "EE1",
            "EEE1",
        ]
        ref = [
            "A",
            "s",
            "A1",
            "A1g",
            "sigma",
            "pi",
            "delta",
            "phi",
            "sigma.g",
            "Ag",
            "A'",
            'A"',
            "E1'",
            'E1"',
        ]
        assert list(map(sym, labels)) == ref

    def test_normalisesym_cfour(self) -> None:
        from cclib.parser.cfourparser import CFOUR

        sym = CFOUR(io.StringIO("dummyfile")).normalisesym
        labels = [
            "A",
            "A1",
            "Ag",
            "A1g",
            "A'",
            "A''",
            "SG",
            "SG+",
            "SGg",
            "SGg+",
            "PI",
            "PIg",
            "DE",
            "DEu",
            "PH",
            "PHu",
            "E",
            "1g",
        ]
        ref = [
            "A",
            "A1",
            "Ag",
            "A1g",
            "A'",
            'A"',
            "sigma",
            "sigma",
            "sigma.g",
            "sigma.g",
            "pi",
            "pi.g",
            "delta",
            "delta.u",
            "phi",
            "phi.u",
            "E",
            "E1g",
        ]
        assert list(map(sym, labels)) == ref

    def test_normalisesym_gamess(self) -> None:
        from cclib.parser.gamessparser import GAMESS

        sym = GAMESS(io.StringIO("dummyfile")).normalisesym
        labels = ["A", "A1", "A1G", "A'", "A''", "AG"]
        ref = ["A", "A1", "A1g", "A'", 'A"', "Ag"]
        assert list(map(sym, labels)) == ref

    def test_normalisesym_gamessuk(self) -> None:
        from cclib.parser.gamessukparser import GAMESSUK

        sym = GAMESSUK(io.StringIO("dummyfile")).normalisesym
        labels = ["a", "a1", "ag", "a'", 'a"', "a''", "a1''", 'a1"', "e1+", "e1-"]
        ref = ["A", "A1", "Ag", "A'", 'A"', 'A"', 'A1"', 'A1"', "E1", "E1"]
        assert list(map(sym, labels)) == ref

    def test_normalisesym_gaussian(self) -> None:
        from cclib.parser.gaussianparser import Gaussian

        sym = Gaussian(io.StringIO("dummyfile")).normalisesym
        labels = ["A1", "AG", "A1G", "SG", "PI", "PHI", "DLTA", "DLTU", "SGG"]
        ref = ["A1", "Ag", "A1g", "sigma", "pi", "phi", "delta", "delta.u", "sigma.g"]
        assert list(map(sym, labels)) == ref

    def test_normalisesym_jaguar(self) -> None:
        from cclib.parser.jaguarparser import Jaguar

        sym = Jaguar(io.StringIO("dummyfile")).normalisesym
        labels = ["A", "A1", "Ag", "Ap", "App", "A1p", "A1pp", "E1pp/Ap"]
        ref = ["A", "A1", "Ag", "A'", 'A"', "A1'", 'A1"', 'E1"']
        assert list(map(sym, labels)) == ref

    def test_normalisesym_molcas(self) -> None:
        from cclib.parser.molcasparser import Molcas

        sym = Molcas(io.StringIO("dummyfile")).normalisesym
        labels = ["a", "a1", "ag"]
        ref = ["A", "A1", "Ag"]
        assert list(map(sym, labels)) == ref

    def test_normalisesym_molpro(self) -> None:
        from cclib.parser.molproparser import Molpro

        sym = Molpro(io.StringIO("dummyfile")).normalisesym
        labels = ["A`", "A``"]
        ref = ["A'", "A''"]
        assert list(map(sym, labels)) == ref

    def test_normalisesym_psi4(self) -> None:
        from cclib.parser.psi4parser import Psi4

        sym = Psi4(io.StringIO("dummyfile")).normalisesym
        labels = ["Ap", "App"]
        ref = ["A'", 'A"']
        assert list(map(sym, labels)) == ref

    def test_normalisesym_turbomole(self) -> None:
        from cclib.parser.turbomoleparser import Turbomole

        sym = Turbomole(io.StringIO("dummyfile")).normalisesym
        labels = ["a", "a1", "ag"]
        ref = ["A", "A1", "Ag"]
        assert list(map(sym, labels)) == ref
