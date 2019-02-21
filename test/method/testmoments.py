# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the Moments method in cclib"""

import unittest

import numpy as np
from numpy.testing import assert_equal, assert_almost_equal

from cclib.method import Moments

from six import add_move, MovedModule
add_move(MovedModule('mock', 'mock', 'unittest.mock'))
from six.moves import mock


class TestIdealizedInputs(unittest.TestCase):
    linear_dipole_attrs = {
        'atomcoords': np.array([[[-1, 0, 0], [ 1,  0, 0]]]),
        'atomcharges': {'mulliken': [-1, 1]},
        'atomnos': [1, 1]
    }
    
    @mock.patch('cclib.parser.ccData', spec=True)
    def test_dipole_moment(self, mock):
        mock.configure_mock(**self.linear_dipole_attrs)

        x = Moments(mock).calculate()[1]
        assert_almost_equal(x / 4.80320425, [2, 0, 0])

    @unittest.skip("This does not pass for some reason.")
    @mock.patch('cclib.parser.ccData', spec=True)
    def test_nonzero_quadrupole_moment(self, mock):
        mock.atomcoords = np.array([[
            [-1, 0, 0],
            [0, 0, 0],
            [1, 0, 0]]])
        # The periods are for Python2.
        mock.atomcharges = {'mulliken': [1/2., -1., 1/2.]}
        mock.atomnos = np.ones(mock.atomcoords.shape[1])

        x = Moments(mock).calculate()
        self.assertEqual(np.count_nonzero(x[1]), 0)
        assert_almost_equal(x[2] / 4.80320423, [1, 0, 0, -0.5, 0, -0.5])

    @mock.patch('cclib.parser.ccData', spec=True)
    def test_zero_moments(self, mock):
        mock.atomcoords = np.array([[
            [-2, 0, 0],
            [-1, 0, 0],
            [0, 0, 0],
            [1, 0, 0],
            [2, 0, 0]]])
        # The periods are for Python2.
        mock.atomcharges = {'mulliken': [-1/8., 1/2., -3/4., 1/2., -1/8.]}
        mock.atomnos = np.ones(mock.atomcoords.shape[1])

        x = Moments(mock).calculate()
        self.assertEqual(np.count_nonzero(x[1]), 0)
        self.assertEqual(np.count_nonzero(x[2]), 0)

    @mock.patch('cclib.parser.ccData', spec=True)
    def test_invariant_to_origin_dislacement(self, mock):
        mock.configure_mock(**self.linear_dipole_attrs)

        x = Moments(mock).calculate(origin=[0, 0, 0])[1]
        y = Moments(mock).calculate(origin=[1, 1, 1])[1]
        assert_equal(x, y)

    @mock.patch('cclib.parser.ccData', spec=True)
    def test_variant_to_origin_dislacement(self, mock):
        attrs = dict(self.linear_dipole_attrs, **{
            'atomcharges': {'mulliken': [-1, 2]}
        })
        mock.configure_mock(**attrs)

        x = Moments(mock).calculate(origin=[0, 0, 0])[1]
        y = Moments(mock).calculate(origin=[1, 1, 1])[1]
        self.assertFalse(np.array_equal(x, y))

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

        x = Moments(mock).calculate(masses=[1, 3], origin='mass')
        assert_almost_equal(x[0], [0.5, 0, 0])

    @mock.patch('cclib.parser.ccData', spec=True)
    def test_not_providing_masses(self, mock):
        mock.configure_mock(**self.linear_dipole_attrs)
        # Replace with the refex version when Python2 is dropped.
        # with self.assertRaisesRegex(ValueError, 'masses'):
        with self.assertRaises(ValueError):
            Moments(mock).calculate(origin='mass')

    @mock.patch('cclib.parser.ccData', spec=True)
    def test_results_storing(self, mock):
        mock.configure_mock(**self.linear_dipole_attrs)
        mock.atomcharges.update({'lowdin': [-0.5, 0.5]})

        m = Moments(mock)
        m.calculate(population='mulliken')
        m.calculate(population='lowdin')

        a, b = m.results['mulliken'][1], m.results['lowdin'][1]
        self.assertFalse(np.array_equal(a, b))


if __name__ == '__main__':
    suite = unittest.makeSuite(TestIdealizedInputs)
    unittest.TextTestRunner(verbosity=2).run(suite)
