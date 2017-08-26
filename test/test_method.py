# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Run method unit tests for cclib."""

import sys
import unittest

sys.path.append("method")
from testcda import *
from testmbo import *
from testnuclear import *
from testorbitals import *
from testpopulation import *


if __name__ == "__main__":
    unittest.main()
