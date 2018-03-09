# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Abstract based class for cclib methods."""

import logging
import sys


class Method(object):
    """Abstract base class for all cclib method classes.

    Subclasses defined by cclib:
        CDA - charde decomposition analysis
        CSPA - C-squared population analysis
        Density - density matrix calculation
        FragmentAnalysis - fragment analysis for ADF output
        LPA - Löwdin population analysis
        MBO - Mayer's bond orders
        MPA - Mulliken population analysis
        Nuclear - properties of atomic nuclei
        OPA - overlap population analysis
        Population - base class for population analyses
        Volume - volume/grid calculations

    All the modules containing methods should be importable:
    >>> import cda, cspa, density, fragments, lpa, mbo, mpa, nuclear, opa, population, volume
    """

    def __init__(self, data, requiredAttr, progress=None, loglevel=logging.INFO, logname="Log"):
        """Initialise the Logfile object.

        This constructor is typically called by the constructor of a subclass.
        """

        self.data = data
        self.requiredAttr = requiredAttr
        self.progress = progress
        self.loglevel = loglevel
        self.logname = logname
        self.initFlag = 0
        self.attrErrorMessage = self.init_error()
        
        self.logger = logging.getLogger('%s %s' % (self.logname, self.data))
        self.logger.setLevel(self.loglevel)
        self.logformat = "[%(name)s %(levelname)s] %(message)s"
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter(self.logformat))
        self.logger.addHandler(handler)

    def init_error(self):

        get = []
        miss = []

        for attr in self.requiredAttr:
            if attr in self.data.__dict__:
                get.append(attr)
            else:
                miss.append(attr)

        if len(miss) == 0:
            pass
            return ''
        else:
            self.initFlag =1
            missingAttr = ','.join(i for i in miss)
            return missingAttr

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=False)
