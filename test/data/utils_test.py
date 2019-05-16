# -*- coding: utf-8 -*-
#
# Copyright (c) 2019, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Tests for data test utilities."""

import mock
import os
import unittest

import cclib

import utils


def load_and_run(test_case_class):
    """Loads and runs test methods from a TestCase class."""
    loader = unittest.loader.TestLoader()
    test_suite = loader.loadTestsFromTestCase(test_case_class)
    results = unittest.TestResult()
    test_suite.run(results)
    return test_suite, results


def get_fake_ccread(filename_to_attribute_dict):
    """Generates a function that fakes the behavior or ccread.
    
    This is for use by tests that want to simulate the behavior of ccread
    and return ccData object with specific attributes set.
    """
    def fake_ccread(path):
        data = cclib.parser.ccData()
        for filename, attribute_dict in filename_to_attribute_dict.items():
            if os.path.basename(path) == filename:
                for key, value in attribute_dict.items():
                    setattr(data, key, value)
                return data
        return data
    return fake_ccread


class DataTestCaseTest(unittest.TestCase):        
    """Test case class used by cclib for testing logfile parsing for multiple programs."""

    def setUp(self):
      self.addCleanup(mock.patch.stopall)

      # We don't actually want to parse anything in this test.
      self.mock_cclib = mock.patch.object(utils, 'cclib').start()
      
      # This mock will let us set unit test path data for specific job
      # types and programs in test methods. Note that overriding the
      # unit_test_paths dictionary will break tests, because it needs to
      # be the same object that the mock references. The default dict
      # below is useful for test that don't realy on the data.
      self.mock_unit_test_paths = mock.patch.object(utils, 'UNIT_TEST_PATHS').start()
      self.unit_test_paths = {
            'some_job_type': {
                'some_program': {
                    'some_version': ('some_file',)
                }
            }
        }
      self.mock_unit_test_paths.__getitem__.side_effect = self.unit_test_paths.__getitem__

    def test_for_job_type(self):
        self.unit_test_paths['some_job_type'] = {
            'some_program': {
                'some_version1': ('some_file1',),
                'some_version2': ('some_file2',),
            },
        }
        class FakeTestCase(utils.DataTestCase):
            @utils.for_job_type('some_job_type')
            def test_pass(self, data): pass
        self.assertIn('test_pass_some_job_type_some_version1', dir(FakeTestCase))
        self.assertIn('test_pass_some_job_type_some_version2', dir(FakeTestCase))

        # Only the job-specific methods get loaded, because the original
        # test_pass method is set to None since the for_job_type decorator
        # doesn't return a callable object.
        test_suite, results = load_and_run(FakeTestCase)
        self.assertEqual(test_suite.countTestCases(), 2)
        self.assertEqual(results.testsRun, 2)
        self.assertEqual(len(results.errors), 0)
        self.assertEqual(len(results.failures), 0)

    def test_for_job_type_with_parameter(self):
        self.unit_test_paths['some_job_type'] = {
            'some_program': {
                'some_version1': ('some_file1',),
                'some_version2': ('some_file2',),
            },
        }
        class FakeTestCase(utils.DataTestCase):
            @utils.for_job_type('some_job_type', some_param=5)
            def test_param(self, data, some_param):
                self.assertEqual(some_param, 5)

        test_suite, results = load_and_run(FakeTestCase)
        self.assertEqual(results.testsRun, 2)
        self.assertEqual(len(results.errors), 0)
        self.assertEqual(len(results.failures), 0)

    def test_for_job_type_with_parameter_and_program_override(self):
        self.unit_test_paths['some_job_type'] = {
            'some_program': {
                'some_version': ('some_file',),
            },
            'other_program': {
                'other_version': ('other_file',),
            },
        }
        self.mock_cclib.io.ccread = get_fake_ccread({
            'some_file': {'some_param': 5},
            'other_file': {'some_param': 6},
        })
        class FakeTestCase(utils.DataTestCase):
            @utils.for_job_type('some_job_type', some_param=5)
            @utils.for_program('other_program', some_param=6)
            def test_param(self, data, some_param):
                self.assertEqual(some_param, data.some_param)

        test_suite, results = load_and_run(FakeTestCase)
        self.assertEqual(results.testsRun, 2)
        self.assertEqual(len(results.errors), 0)
        self.assertEqual(len(results.failures), 0)

    def test_for_job_type_with_multiple_job_types(self):
        self.unit_test_paths['some_job_type'] = {
            'some_program': {
                'some_version1': ('some_file1',),
                'some_version2': ('some_file2',),
            },
        }
        self.unit_test_paths['other_job_type'] = {
            'some_program': {
                'some_version3': ('some_file3',),
                'some_version4': ('some_file4',),
            },
        }
        class FakeTestCase(utils.DataTestCase):
            @utils.for_job_type('some_job_type', 'other_job_type')
            def test_pass(self, data): pass

        test_suite, results = load_and_run(FakeTestCase)
        self.assertEqual(results.testsRun, 4)
        self.assertEqual(len(results.errors), 0)
        self.assertEqual(len(results.failures), 0)

    def test_for_job_type_with_multiple_test_methods(self):
        self.unit_test_paths['some_job_type'] = {
            'some_program': {
                'some_version1': ('some_file1',),
                'some_version2': ('some_file2',),
            },
        }
        self.unit_test_paths['other_job_type'] = {
            'some_program': {
                'some_version2': ('some_file2',),
                'some_version3': ('some_file3',),
                'some_version4': ('some_file4',),
            },
        }
        class FakeTestCase(utils.DataTestCase):
            @utils.for_job_type('some_job_type')
            def test_pass1(self, data): pass
            @utils.for_job_type('some_job_type', 'other_job_type')
            def test_pass2(self, data): pass

        # We need to start "parsing" from scratch in order to make
        # meaningful assertions about it below.
        FakeTestCase._logfile_data_cache.clear()

        test_suite, results = load_and_run(FakeTestCase)
        self.assertEqual(results.testsRun, 7)
        self.assertEqual(len(results.errors), 0)
        self.assertEqual(len(results.failures), 0)

        # There are five job-program-filename combinations in the paths
        # dictionary, but the cache should eliminate one call to ccread,
        # because only program-filename combos have different data.
        self.assertEqual(self.mock_cclib.io.ccread.call_count, 4)

    def test_for_program_without_for_job_type(self):
        class FakeTestCase(utils.DataTestCase):
            def runTest(self): pass
            @utils.for_program('some_program')
            def test_pass(self): pass
        self.assertIn('_test_pass', dir(FakeTestCase))

        # The test_pass method is wrapped and set, but it will error because
        # there are no program-specific test methods and the for_program decorator
        # only makes sense in that context.
        test_suite, results = load_and_run(FakeTestCase)
        self.assertEqual(test_suite.countTestCases(), 1)
        self.assertEqual(results.testsRun, 1)
        self.assertEqual(len(results.errors), 1)
        self.assertIn('Use for_job_type decorator before for_program',
                      results.errors[0][1])

    def test_for_job_type_before_for_program(self):
        with self.assertRaises(TypeError):
            class FakeTestCase(utils.DataTestCase):
                def runTest(self): pass
                @utils.for_program('some_program')
                @utils.for_job_type('some_job_type')
                def test_pass(self): pass

    def test_for_job_type_with_different_class(self):
        class FakeTestCase(unittest.TestCase):
            @utils.for_job_type('some_job_type')
            def test_pass(self): pass

        test_suite, results = load_and_run(FakeTestCase)
        self.assertEqual(results.testsRun, 1)
        self.assertEqual(len(results.errors), 1)
        self.assertIn('method assumes it will be bound to a cclib DataTestCase',
                      results.errors[0][1])


if __name__ == '__main__':
    unittest.main()
