# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Contains all writers for standard chemical representations."""

from cclib.io.cjsonreader import CJSON as CJSONReader
from cclib.io.cjsonwriter import CJSON as CJSONWriter
from cclib.io.cmlwriter import CML
from cclib.io.moldenwriter import MOLDEN
from cclib.io.wfxwriter import WFXWriter
from cclib.io.xyzreader import XYZ as XYZReader
from cclib.io.xyzwriter import XYZ as XYZWriter

# This allows users to type:
#   from cclib.io import ccframe
#   from cclib.io import ccopen
#   from cclib.io import ccread
#   from cclib.io import ccwrite
#   from cclib.io import URL_PATTERN
from cclib.io.ccio import ccframe
from cclib.io.ccio import ccopen
from cclib.io.ccio import ccread
from cclib.io.ccio import ccwrite
from cclib.io.ccio import URL_PATTERN
