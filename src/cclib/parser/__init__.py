"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision$"

# These import statements are added for the convenience of users...

# Rather than having to type:
#         from cclib.parser.gaussianparser import Gaussian
# they can use:
#         from cclib.parser import Gaussian

from gaussianparser import Gaussian
from gamessparser import GAMESS
from adfparser import ADF
from jaguarparser import Jaguar
from gamessukparser import GAMESSUK
from molproparser import Molpro

# This allow users to type:
#         from cclib.parser import ccopen

from utils import ccopen
