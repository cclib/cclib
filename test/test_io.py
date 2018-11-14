# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Run writer unit tests for cclib."""

import sys
import unittest

sys.path.insert(1, 'io')

from .io.testccio import *
from .io.testcjsonreader import *
from .io.testcjsonwriter import *
from .io.testfilewriter import *
from .io.testmoldenwriter import *
from .io.testscripts import *
from .io.testwfxwriter import *
from .io.testxyzreader import *
from .io.testxyzwriter import *

if __name__ == "__main__":
    unittest.main()
