# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.


__revision__ = "$Revision$"

# These import statements are added for the convenience of users...

# Rather than having to type:
#         from cclib.parser.gaussianparser import Gaussian
# they can use:
#         from cclib.parser import Gaussian

from adfparser import ADF
from gamessparser import GAMESS
from gamessukparser import GAMESSUK
from gaussianparser import Gaussian
from jaguarparser import Jaguar
from molproparser import Molpro
from orcaparser import ORCA

# This allow users to type:
#         from cclib.parser import ccopen

from ccopen import ccopen
