# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the Moments method in cclib"""

from __future__ import print_function

import sys
import unittest
from unittest import mock

import numpy as np
from numpy.testing import assert_equal, assert_almost_equal

from cclib.method import Moments


class TestIdealizedInputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.linear_dipole_attrs = {
            'atomcoords': np.array([[[-1, 0, 0], [ 1,  0, 0]]]),
            'atomcharges': {'mulliken': [-1, 1]},
            'atomnos': [1, 1],
            'charge': 0,
        }
    
    @mock.patch('cclib.parser.ccData', spec=True)
    def test_dipole_moment(self, mock):
        mock.configure_mock(**self.linear_dipole_attrs)

        x = Moments(mock).calculate()[1]
        assert_almost_equal(x / 4.803204252680386, [2, 0, 0])

    @mock.patch('cclib.parser.ccData', spec=True)
    def test_quarupole_moment(self, mock):
        mock.charge = 0
        mock.atomcoords = np.array([[
            [-0.75,  1, 0],
            [ 1,  1, 0],
            [ 0.75, -1, 0],
            [-1, -1, 0]]])
        mock.atomcharges = {'mulliken': [-1, 1, -1, 1]}
        mock.atomnos = np.ones(mock.atomcoords.shape[1])

        x = Moments(mock).calculate()[2]
        assert np.isclose(x[0] + x[3] + x[5], 0)
        assert np.isclose(x[0], -x[3] * 2)
        assert np.isclose(x[0], -x[5] * 2)

    @mock.patch('cclib.parser.ccData', spec=True)
    def test_invariant_to_origin_dislacement(self, mock):
        mock.configure_mock(**self.linear_dipole_attrs)

        x = Moments(mock).calculate(origin=[0, 0, 0])[1]
        y = Moments(mock).calculate(origin=[1, 1, 1])[1]
        assert_equal(x, y)
    
    @mock.patch('cclib.parser.ccData', spec=True)
    def test_variant_to_origin_dislacement(self, mock):
        attrs = dict(self.linear_dipole_attrs, **{
            'charge': 1,
            'atomcharges': {'mulliken': [-1, 2]}
        })
        mock.configure_mock(**attrs)

        x = Moments(mock).calculate(origin=[0, 0, 0])[1]
        y = Moments(mock).calculate(origin=[1, 1, 1])[1]
        assert not np.array_equal(x, y)

    @mock.patch('cclib.parser.ccData', spec=True)
    def test_origin_at_center_of_nuclear_charge(self, mock):
        mock.configure_mock(**self.linear_dipole_attrs)

        x = Moments(mock).calculate()[0]
        assert_equal(x, [0, 0, 0])
        
    @mock.patch('cclib.parser.ccData', spec=True)
    def test_origin_at_center_of_mass(self, mock):
        mock.configure_mock(**self.linear_dipole_attrs)
        mock.atommasses = np.ones(mock.atomcoords.shape[1])

        x = Moments(mock).calculate(origin='mass')[0]
        assert_equal(x, [0, 0, 0])

    @mock.patch('cclib.parser.ccData', spec=True)
    def test_user_provided_origin(self, mock):
        mock.configure_mock(**self.linear_dipole_attrs)

        x = Moments(mock).calculate(origin=[1, 1, 1])
        assert_almost_equal(x[0], [1, 1, 1])

    @mock.patch('cclib.parser.ccData', spec=True)
    def test_user_provided_masses(self, mock):
        mock.configure_mock(**self.linear_dipole_attrs)

        x = Moments(mock).calculate(masses=[1, 1], origin='mass')
        assert_almost_equal(x[0], [0, 0, 0])

    @mock.patch('cclib.parser.ccData', spec=True)
    def test_not_providing_masses(self, mock):
        mock.configure_mock(**self.linear_dipole_attrs)
        with self.assertRaisesRegex(ValueError, 'masses'):
            Moments(mock).calculate(origin='mass')

    @mock.patch('cclib.parser.ccData', spec=True)
    def test_results_storing(self, mock):
        mock.configure_mock(**self.linear_dipole_attrs)
        mock.atomcharges.update({'lowdin': [-0.5, 0.5]})
        
        m = Moments(mock)
        m.calculate(population='mulliken')
        m.calculate(population='lowdin')

        a, b = m.results['mulliken'][1], m.results['lowdin'][1]
        assert not np.array_equal(a, b)

        
if __name__ == '__main__':
    suite = unittest.makeSuite(TestIdealizedInputs)
    unittest.TextTestRunner(verbosity=2).run(suite)
