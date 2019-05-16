# -*- coding: utf-8 -*-
#
# Copyright (c) 2019, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Helpers used for data tests."""

import inspect
import os
import unittest

import cclib

from paths import UNIT_TEST_PATHS


class DataTestCase(unittest.TestCase):
    """A test case class the support parameterizing with job types.
    
    The purpose of this class is to make it easier to write compact tests
    for cclib data unit tests. There are two basic features on top of the
    base TestCase class:
      1. The ability to parameterize test methods with job types such as
         SP, GOpt, which generates a test method for each program version
         that has a logfile in the UNIT_TEST_PATHS dictionary. Additional
         parameters may be passed to the test method as well.
      2. The ability to override a test parameter for specific programs.
    
    Below is an example test case with two methods. The first test method
    applies to two different job types, and all programs for which we have
    log files for those job types. The second one uses a different expected
    value for one program.
    
    class NAtomsTestCase(DataTestCase):

        @for_job_type("SP", "GOpt", expected_natoms=5)
        test_value(self, data, expected_natoms):
            self.assertEqual(data.natoms, expected_natoms)

        @for_job_type("vib", expected_natoms=5)
        @for_program("Psi", expected_natoms=6, msg="Some mix up in the input file")
        test_value(self, data, expected_natoms):
            self.assertEqual(data.natoms, expected_natoms)

    Note that the for_job_type decorator needs to always be on opt (last).
    """

    # This is a class member, which means any unit tests inheriting from
    # it will access and write to the same cache. And this is desirable,
    # since we only want to parse each logfile once regardless of the
    # scope of the test run.
    _logfile_data_cache = {}

    def setUp(self):
        self.current_program = None

    def _get_data_for_logfile_path(self, path):
        if path not in self._logfile_data_cache:
            data = cclib.io.ccread(path)
            if not data:
                raise IOError("Could not parse: %s" % path)
            self._logfile_data_cache[path] = data
        return self._logfile_data_cache[path]


def for_job_type(*job_types, **kwargs):
    """Returns a decorator for parameterizing a DataTestCase method with one or more job types."""
    def wrapper(wrapped_func):
        new_methods = {}
        for jtype, program, version, file_name in _flattened_unit_test_paths(job_types):
            path = os.path.join('data', program, version, file_name)
            original_name = wrapped_func.__name__.lstrip('_')
            new_func = _standalone_job_type_func(path, wrapped_func, program, original_name, **kwargs)
            new_func.__test__ = True
            new_name = '_'.join((original_name, jtype, version))
            new_methods[new_name] = new_func

        stack = inspect.stack()
        frame = stack[1]
        frame_locals = frame[0].f_locals
        for new_name, new_func in new_methods.items():
            frame_locals[new_name] = new_func

        wrapped_func.__test__ = False

    return wrapper


def for_program(program_type, msg=None, **kwargs):
    """Returns a decorator for modifying a DataTestCase method with one or more program-specific changes."""
    def decorator(wrapped_func):
        if not wrapped_func:
            raise TypeError('The for_job_type decorator needs to be applied last.')

        new_func = _standalone_program_func(wrapped_func, program_type, **kwargs)
        new_func.__test__ = True
        new_name = '_' + wrapped_func.__name__
        new_func.__name__ = new_name

        stack = inspect.stack()
        frame = stack[1]
        frame_locals = frame[0].f_locals
        frame_locals[new_name] = new_func

        wrapped_func.__test__ = False
        return new_func
    return decorator


def _standalone_job_type_func(path, wrapped_func, program, original_name, **kwargs):
  """Returns a new function that wraps a job-specific test method."""
  def standalone_func(self):
    if not isinstance(self, DataTestCase):
        raise TypeError("This method assumes it will be bound to a cclib DataTestCase.")
    data = self._get_data_for_logfile_path(path)
    self.current_program = program
    wrapped_func(self, data, **kwargs)
  return standalone_func


def _standalone_program_func(wrapped_func, program_type, **wrapped_kwargs):
  """Returns a new function that wraps a program-specific test method."""
  def standalone_func(self, *args, **kwargs):
    if not isinstance(self, DataTestCase):
        raise TypeError('This method assumes it will be bound to a cclib DataTestCase.')
    new_kwargs = kwargs
    if not self.current_program:
        raise ValueError('No program set. Use for_job_type decorator before for_program.')
    if self.current_program == program_type:
        new_kwargs = _merge_two_dicts(kwargs, wrapped_kwargs)
    wrapped_func(self, *args, **new_kwargs)
  return standalone_func


def _flattened_unit_test_paths(job_types):
  """Generates flattened tuples from the unit test path dictionary of selected job types."""
  for jtype in job_types:
      for program, program_files in UNIT_TEST_PATHS[jtype].items():
          for version, version_files in program_files.items():
              for file_name in version_files:
                  yield jtype, program, version, file_name


def _merge_two_dicts(x, y):
    z = x.copy()
    z.update(y)
    return z
