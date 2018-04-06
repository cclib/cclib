# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Run parser unit tests for cclib."""

import sys
import unittest

sys.path.insert(1, 'parser')

from .parser.testdata import *
from .parser.testlogfileparser import *
from .parser.testspecificparsers import *
from .parser.testutils import *


if __name__ == "__main__":
    unittest.main()
