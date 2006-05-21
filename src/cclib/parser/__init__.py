# These import statements are added for the convenience of users...

# Rather than having to type:
#         from cclib.parser.gaussianparser import Gaussian
# they can use:
#         from cclib.parser import Gaussian

from gaussianparser import Gaussian
from gamessparser import GAMESS
from adfparser import ADF
from jaguarparser import Jaguar

# This allow users to type:
#         from cclib.parser import guesstype

from utils import guesstype
