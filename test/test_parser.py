# -*- coding: utf-8 -*-
#
# Copyright (c) 2016, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Run parser unit tests for cclib."""

import sys
import unittest

sys.path.append('parser')
from testdata import *
from testlogfileparser import *
from testutils import *


if __name__ == "__main__":
    unittest.main()
