# -*- coding: utf-8 -*-
#
# Copyright (c) 2019, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Run data tests for cclib."""

import unittest

from data.test_atomcoords import *
from data.test_scfenergies import *
from data.test_time import *


if __name__ == '__main__':
    unittest.main()
