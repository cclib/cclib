# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2013 the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

__version__ = "1.2b"

from . import parser
from . import progress
from . import method
from . import bridge

# The test module can be imported if it was installed with cclib.
try:
    from . import test
except ImportError:
    pass
