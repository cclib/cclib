# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Run all unit tests for cclib."""

import unittest

from . import test_data

from .test_bridge import *
from .test_io import *
from .test_method import *
from .test_parser import *
from .test_utils import *


if __name__ == "__main__":
    print("Running unit tests for data...")
    test_data.test_all(silent=True, summary=False, visual_tests=False)
    print("Running all other unit tests...")
    unittest.main()
