# -*- coding: utf-8 -*-
#
# Copyright (c) 2016, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Run method unit tests for cclib."""

import sys
import unittest


if __name__ == "__main__":

    sys.path.append("method")
    from testcda import *
    from testmbo import *
    from testnuclear import *
    from testpopulation import *
    unittest.main()
