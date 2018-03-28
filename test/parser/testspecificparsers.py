# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for specific parser behaviors, such as overriden methods."""

import unittest


class NormalisesymTest(unittest.TestCase):

    # Not needed: DALTON, MOPAC, NWChem, ORCA, Psi, QChem

    def test_normalisesym_adf(self):
        from cclib.parser.adfparser import ADF
        sym = ADF("dummyfile").normalisesym
        labels = ['A', 's', 'A1', 'A1.g', 'Sigma', 'Pi', 'Delta', 'Phi', 'Sigma.g', 'A.g', 'AA', 'AAA', 'EE1', 'EEE1']
        ref = ['A', 's', 'A1', 'A1g', 'sigma', 'pi', 'delta', 'phi', 'sigma.g', 'Ag', "A'", 'A"', "E1'", 'E1"']
        self.assertEqual(list(map(sym, labels)), ref)

    def test_normalisesym_gamess(self):
        from cclib.parser.gamessparser import GAMESS
        sym = GAMESS("dummyfile").normalisesym
        labels = ['A', 'A1', 'A1G', "A'", "A''", "AG"]
        ref = ['A', 'A1', 'A1g', "A'", 'A"', 'Ag']
        self.assertEqual(list(map(sym, labels)), ref)

    def test_normalisesym_gamessuk(self):
        from cclib.parser.gamessukparser import GAMESSUK
        sym = GAMESSUK("dummyfile.txt").normalisesym
        labels = ['a', 'a1', 'ag', "a'", 'a"', "a''", "a1''", 'a1"', "e1+", "e1-"]
        ref = ['A', 'A1', 'Ag', "A'", 'A"', 'A"', 'A1"', 'A1"', 'E1', 'E1']
        self.assertEqual(list(map(sym, labels)), ref)

    def test_normalisesym_gaussian(self):
        from cclib.parser.gaussianparser import Gaussian
        sym = Gaussian("dummyfile").normalisesym
        labels = ['A1', 'AG', 'A1G', "SG", "PI", "PHI", "DLTA", 'DLTU', 'SGG']
        ref = ['A1', 'Ag', 'A1g', 'sigma', 'pi', 'phi', 'delta', 'delta.u', 'sigma.g']
        self.assertEqual(list(map(sym, labels)), ref)

    def test_normalisesym_jaguar(self):
        from cclib.parser.jaguarparser import Jaguar
        sym = Jaguar("dummyfile").normalisesym
        labels = ['A', 'A1', 'Ag', 'Ap', 'App', "A1p", "A1pp", "E1pp/Ap"]
        ref = ['A', 'A1', 'Ag', "A'", 'A"', "A1'", 'A1"', 'E1"']
        self.assertEqual(list(map(sym, labels)), ref)

    def test_normalisesym_molpro(self):
        from cclib.parser.molproparser import Molpro
        sym = Molpro("dummyfile").normalisesym
        labels = ["A`", "A``"]
        ref = ["A'", "A''"]
        self.assertEqual(list(map(sym, labels)), ref)


if __name__ == "__main__":
    unittest.main()
