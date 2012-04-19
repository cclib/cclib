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

import logging
import sys


class Method(object):
    """Abstract class for logfile objects.

    Subclasses defined by cclib:
        Density, Fragments, OPA, Population
    
    Attributes:
        data - ccData source data object
    """
    def __init__(self, data, progress=None,
                 loglevel=logging.INFO, logname="Log"):
        """Initialise the Logfile object.

        Typically called by subclasses in their own __init__ methods.
        """

        self.data = data
        self.progress = progress
        self.loglevel = loglevel
        self.logname = logname

        # Set up the logger.
        self.logger = logging.getLogger('%s %s' % (self.logname, self.data))
        self.logger.setLevel(self.loglevel)
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter(
                             "[%(name)s %(levelname)s] %(message)s"))
        self.logger.addHandler(handler)


if __name__ == "__main__":
    import doctest, calculationmethod
    doctest.testmod(calculationmethod, verbose=False)
