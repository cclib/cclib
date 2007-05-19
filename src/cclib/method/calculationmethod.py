"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision$"

import logging
import sys

class Method(object):
    """Abstract class for logfile objects.

    Subclasses:
        Density
    
    Attributes:
    """
    def __init__(self, parser, progress=None,
                 loglevel=logging.INFO, logname="Log"):
        """Initialise the Logfile object.

        Typically called by subclasses in their own __init__ methods.
        """
        self.parser = parser
        self.progress = progress
        self.loglevel = loglevel
        self.logname = logname

        # Set up the logger
        self.logger = logging.getLogger('%s %s' % (self.logname, self.parser))
        self.logger.setLevel(self.loglevel)
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter(
                             "[%(name)s %(levelname)s] %(message)s"))
        self.logger.addHandler(handler)


if __name__ == "__main__":
    import doctest, calculationmethod
    doctest.testmod(calculationmethod, verbose=False)
