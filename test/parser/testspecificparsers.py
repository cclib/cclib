# Copyright (c) 2025-2026, the cclib development team
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


def _nwchem_n2_output(symbol: str) -> str:
    """A minimal NWChem logfile for N2, with atoms named as typed in the input.

    Modeled on the output attached to https://github.com/cclib/cclib/issues/1801;
    NWChem echoes atom names with the capitalization the user typed (e.g. "n"),
    and D4h is a non-abelian group whose order (16) exceeds its irrep count (10).
    """
    return f"""\
                             Geometry "geometry" -> ""
                             -------------------------

 Output coordinates in angstroms (scale by  1.889725989 to convert to a.u.)

  No.       Tag          Charge          X              Y              Z
 ---- ---------------- ---------- -------------- -------------- --------------
    1 {symbol}                    7.0000     0.00000000     0.00000000    -0.54000000
    2 {symbol}                    7.0000     0.00000000     0.00000000     0.54000000

      Atomic Mass
      -----------

      {symbol}                 14.003070


      Symmetry information
      --------------------

 Group name             D4h
 Group number             28
 Group order              16
 No. of unique centers     1

 Summary of "ao basis" -> "ao basis" (cartesian)
 ------------------------------------------------------------------------------
       Tag                 Description            Shells   Functions and Types
 ---------------- ------------------------------  ------  ---------------------
 {symbol}                          cc-pvdz                  6       15   3s2p1d


      Symmetry analysis of basis
      --------------------------

        a1g         7
        a1u         0
        a2g         0
        a2u         7
        b1g         1
        b1u         1
        b2g         1
        b2u         1
        eg          6
        eu          6

 Forming initial guess at       0.0s
"""


class NWChemTest:
    def test_atom_names_as_typed_in_input(self) -> None:
        """NWChem echoes atom names as typed (e.g. lowercase "n"), while masses,
        atombasis and aonames are resolved via standard PeriodicTable symbols (#1801)."""
        from cclib.parser.nwchemparser import NWChem

        data = NWChem(io.StringIO(_nwchem_n2_output("n"))).parse()
        assert data.atommasses.tolist() == [14.003070, 14.003070]
        assert data.atombasis == [list(range(15)), list(range(15, 30))]
        assert data.aonames[:4] == ["N1_1S", "N1_2S", "N1_3S", "N1_2PX"]
        assert len(data.aonames) == 30

    def test_symlabels_nonabelian_group(self) -> None:
        """The symmetry analysis block has one line per irrep, which is fewer
        than the group order for non-abelian groups (D4h: 10 irreps, order 16)."""
        from cclib.parser.nwchemparser import NWChem

        data = NWChem(io.StringIO(_nwchem_n2_output("N"))).parse()
        assert data.atommasses.tolist() == [14.003070, 14.003070]
