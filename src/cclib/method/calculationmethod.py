"""
cclib is a parser for computational chemistry log files.

See http://cclib.sf.net for more information.

Copyright (C) 2006 Noel O'Boyle and Adam Tenderholt

 This program is free software; you can redistribute and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY, without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

Contributions (monetary as well as code :-) are encouraged.
"""
import logging, sys
import Numeric

class Method(object):
    """Abstract class for logfile objects.

    Subclasses:
        Density
    
    Attributes:
    """
    def __init__(self,parser,progress=None,
                 loglevel=logging.INFO,logname="Log"):
        """Initialise the Logfile object.

        Typically called by subclasses in their own __init__ methods.
        """
        self.parser = parser
        self.progress = progress
        self.loglevel = loglevel
        self.logname  = logname

        # Set up the logger
        self.logger = logging.getLogger('%s %s' % (self.logname,self.parser))
        self.logger.setLevel(self.loglevel)
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter("[%(name)s %(levelname)s] %(message)s"))
        self.logger.addHandler(handler)


if __name__=="__main__":
    import doctest,calculationmethod
    doctest.testmod(calculationmethod,verbose=False)
