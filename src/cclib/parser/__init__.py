# -*- coding: utf-8 -*-
#
# Copyright (c) 2016, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

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

# This allows users to type:
#         from cclib.parser import ccopen
from ..io.ccio import ccopen
