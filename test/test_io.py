# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Run writer unit tests for cclib."""

import sys
import unittest

sys.path.append('io')
from testccio import *
from testfilewriter import *
from testxyzwriter import *
from testcjsonreader import *
from testcjsonwriter import *
from testmoldenwriter import *

from testwfxwriter import *

if __name__ == "__main__":
    unittest.main()
