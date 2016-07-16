# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2014-2016, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Contains all writers for standard chemical representations"""

from .cjsonwriter import CJSON
from .cmlwriter import CML
from .xyzwriter import XYZ

# This allows users to type:
#   from cclib.io import ccopen
#   from cclib.io import ccread
#   from cclib.io import ccwrite
from .ccio import ccopen
from .ccio import ccread
from .ccio import ccwrite

