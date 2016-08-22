# -*- coding: utf-8 -*-
#
# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2016, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Contains parsers for all supported programs"""


# These import statements are added for the convenience of users...
# Rather than having to type:
#         from cclib.parser.gaussianparser import Gaussian
# they can use:
#         from cclib.parser import Gaussian

from .adfparser import ADF
from .daltonparser import DALTON
from .gamessparser import GAMESS
from .gamessukparser import GAMESSUK
from .gaussianparser import Gaussian
from .jaguarparser import Jaguar
from .molproparser import Molpro
from .mopacparser import MOPAC
from .nwchemparser import NWChem
from .orcaparser import ORCA
from .psiparser import Psi
from .qchemparser import QChem

from .data import ccData