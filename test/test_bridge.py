# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Run bridge unit tests for cclib."""

import sys
import unittest

sys.path.insert(1, "bridge")

if sys.version_info[0] == 3 and sys.version_info[1] >= 4:
    from .bridge.testbiopython import *
from .bridge.testopenbabel import *
if sys.version_info[0] == 2:
    from .bridge.testpyquante import *


if __name__ == "__main__":
    unittest.main()
