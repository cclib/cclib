# -*- coding: utf-8 -*-
#
# Copyright (c) 2016, the cclib development team
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

    def __init__(self, data, progress=None, loglevel=logging.INFO, logname="Log"):
        """Initialise the Logfile object.

        This constructor is typically called by the constructor of a subclass.
        """

        self.data = data
        self.progress = progress
        self.loglevel = loglevel
        self.logname = logname

        self.logger = logging.getLogger('%s %s' % (self.logname, self.data))
        self.logger.setLevel(self.loglevel)
        self.logformat = "[%(name)s %(levelname)s] %(message)s"
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter(self.logformat))
        self.logger.addHandler(handler)


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=False)
