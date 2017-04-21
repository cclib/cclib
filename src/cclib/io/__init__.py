# -*- coding: utf-8 -*-
#
# Copyright (c) 2016, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Contains all writers for standard chemical representations"""

from .cjsonwriter import CJSON as CJSONWriter
from .cmlwriter import CML
from .xyzwriter import XYZ
from .cjsonreader import CJSON as CJSONReader

# This allows users to type:
#   from cclib.io import ccopen
#   from cclib.io import ccread
#   from cclib.io import ccwrite
#   from cclib.io import URL_PATTERN
from .ccio import ccopen
from .ccio import ccread
from .ccio import ccwrite
from .ccio import URL_PATTERN
