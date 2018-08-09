# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test metadata in cclib."""

import unittest

from skip import skipForParser


class MetadataTest(unittest.TestCase):

    @skipForParser('ADF', 'Does not fully support metadata yet')
    @skipForParser('DALTON', 'Does not fully support metadata yet')
    @skipForParser('GAMESS', 'Does not fully support metadata yet')
    @skipForParser('GAMESSUK', 'Does not fully support metadata yet')
    @skipForParser('Gaussian', 'Does not fully support metadata yet')
    @skipForParser('Jaguar', 'Does not fully support metadata yet')
    @skipForParser('Molcas', 'Does not full support metadata yet')
    @skipForParser('Molpro', 'Does not full support metadata yet')
    @skipForParser('MOPAC', 'Does not full support metadata yet')
    @skipForParser('NWChem', 'Does not full support metadata yet')
    @skipForParser('ORCA', 'Does not full support metadata yet')
    @skipForParser('Psi3', 'Does not fully support metadata yet')
    @skipForParser('Psi4', 'Does not fully support metadata yet')
    @skipForParser('QChem', 'Does not fully support metadata yet')
    @skipForParser('Turbomole', 'Does not fully support metadata yet')
    def testkeys(self):
        """Does metadata have expected keys and values?"""
        self.assertTrue(hasattr(self.data, "metadata"))
        if self.logfile.logname not in ['ORCA', 'Psi']:
            self.assertIn("basis_set", self.data.metadata)
        if self.logfile.logname == 'ORCA':
            self.assertIn("input_file_name", self.data.metadata)
            self.assertIn("input_file_contents", self.data.metadata)
        self.assertIn("methods", self.data.metadata)
        self.assertIn("package", self.data.metadata)
        self.assertIn("package_version", self.data.metadata)
        self.assertIn("success", self.data.metadata)

    @skipForParser('ADF', 'Does not fully support metadata yet')
    @skipForParser('DALTON', 'Does not fully support metadata yet')
    @skipForParser('GAMESS', 'Does not fully support metadata yet')
    @skipForParser('GAMESSUK', 'Does not fully support metadata yet')
    @skipForParser('Gaussian', 'Does not fully support metadata yet')
    @skipForParser('Jaguar', 'Does not fully support metadata yet')
    @skipForParser('Molcas', 'Does not full support metadata yet')
    @skipForParser('Molpro', 'Does not full support metadata yet')
    @skipForParser('MOPAC', 'Does not full support metadata yet')
    @skipForParser('NWChem', 'Does not full support metadata yet')
    @skipForParser('ORCA', 'Does not full support metadata yet')
    @skipForParser('Psi3', 'Does not fully support metadata yet')
    @skipForParser('Psi4', 'Does not fully support metadata yet')
    @skipForParser('QChem', 'Does not fully support metadata yet')
    @skipForParser('Turbomole', 'Does not fully support metadata yet')
    def testtypes(self):
        self.assertIsInstance(self.data.metadata, dict)
        if self.logfile.logname not in ['ORCA', 'Psi']:
            self.assertIsInstance(self.data.metadata["basis_set"], str)
        if self.logfile.logname == 'ORCA':
            self.assertIsInstance(self.data.metadata["input_file_name"], str)
            self.assertIsInstance(self.data.metadata["input_file_contents"], str)
        self.assertIsInstance(self.data.metadata["methods"], list)
        for method in self.data.metadata["methods"]:
            self.assertIsInstance(method, str)
        self.assertIsInstance(self.data.metadata["package"], str)
        self.assertIsInstance(self.data.metadata["package_version"], str)
        self.assertIsInstance(self.data.metadata["success"], bool)
