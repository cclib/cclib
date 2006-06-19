__revision__ = "$Revision$"

from textprogress import TextProgress
import sys

if 'qt' in sys.modules.keys():
    from qtprogress import QtProgress
