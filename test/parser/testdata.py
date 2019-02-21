# -*- coding: utf-8 -*-
#
# Copyright (c) 2019, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for parser data module."""

import unittest

import numpy

import cclib

from six import add_move, MovedModule
add_move(MovedModule('mock', 'mock', 'unittest.mock'))
from six.moves import mock


class StringContains(str):
    def __eq__(self, other):
        return self in other


class ccDataTest(unittest.TestCase):
    """Unit tests for the ccData class."""

    check_array = ['atomcoords', 'scfenergies']
    check_arrlist = ['mocoeffs', 'scfvalues']
    check_arrdict = ['atomcharges', 'atomspins']

    def setUp(self):
        self.data = cclib.parser.logfileparser.ccData()
        self.mock_logger = mock.Mock()

    def _set_attributes(self, names, val):
        for n in names:
            setattr(self.data, n, val)

    def test_valuecheck_negative_etenergies(self):
        self.data.etenergies = [1, -1]
        self.data.check_values(logger=self.mock_logger)
        self.mock_logger.error.assert_called_once_with(
            StringContains("At least one excitation energy is negative"))

    def test_listify_ndarray(self):
        """Does the method convert ndarrays as expected?"""
        self._set_attributes(self.check_array, numpy.array([1,2,3]))
        self.data.listify()
        for attr in self.check_array:
            self.assertIsInstance(getattr(self.data, attr), list)

    def test_listify_arraylist(self):
        """Does the method convert lists of arrays as expected?"""
        self._set_attributes(self.check_arrlist, [numpy.array([1,2,3]), numpy.array([4,5,6])])
        self.data.listify()
        for attr in self.check_arrlist:
            for a in getattr(self.data, attr):
                self.assertIsInstance(a, list)

    def test_listify_arraydict(self):
        """Does the method convert dicts of arrays as expected?"""
        self._set_attributes(self.check_arrdict, {1: numpy.array([1,2,3]), 2: numpy.array([4,5,6])})
        self.data.listify()
        for attr in self.check_arrdict:
            for a in getattr(self.data, attr).values():
                self.assertIsInstance(a, list)

    def test_arrayify_ndarray(self):
        """Does the method convert lists as expected?"""
        self._set_attributes(self.check_array, [1,2,3])
        self.data.arrayify()
        for attr in self.check_array:
            self.assertIsInstance(getattr(self.data, attr), numpy.ndarray)

    def test_arrayify_arraylist(self):
        """Does the method convert lists of lists as expected?"""
        self._set_attributes(self.check_arrlist, [[1,2,3], [4,5,6]])
        self.data.arrayify()
        for attr in self.check_arrlist:
            for a in getattr(self.data, attr):
                self.assertIsInstance(a, numpy.ndarray)

    def test_arrayify_arraydict(self):
        """Does the method convert dicts of lists as expected?"""
        self._set_attributes(self.check_arrdict, {1: [1,2,3], 2: [4,5,6]})
        self.data.arrayify()
        for attr in self.check_arrdict:
            for a in getattr(self.data, attr).values():
                self.assertIsInstance(a, numpy.ndarray)


if __name__ == "__main__":
    unittest.main()
