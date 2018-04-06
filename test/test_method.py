# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Run method unit tests for cclib."""

import sys
import unittest

sys.path.insert(1, "method")

from .method.testcda import *
from .method.testmbo import *
from .method.testnuclear import *
from .method.testorbitals import *
from .method.testpopulation import *
if sys.version_info[0] == 2:
    from .method.testvolume import *


if __name__ == "__main__":
    unittest.main()
