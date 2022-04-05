# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Abstract based class for cclib methods."""

import logging
import sys

class MissingAttributeError(Exception):
    pass

class Method:
    """Abstract base class for all cclib method classes.

    Subclasses defined by cclib:
        CDA - charde decomposition analysis
        CSPA - C-squared population analysis
        Density - density matrix calculation
        FragmentAnalysis - fragment analysis for ADF output
        LPA - Löwdin population analysis
        MBO - Mayer's bond orders
        Moments - multipole moments calculations
        MPA - Mulliken population analysis
        Nuclear - properties of atomic nuclei
        OPA - overlap population analysis
        Population - base class for population analyses
        Volume - volume/grid calculations

    All the modules containing methods should be importable.
    """
    required_attrs = ()
    def __init__(self, data, progress=None, loglevel=logging.INFO, logname="Log"):
        """Initialise the Logfile object.

        This constructor is typically called by the constructor of a subclass.
        """

        self.data = data
        self.progress = progress
        self.loglevel = loglevel
        self.logname = logname
        self._check_required_attributes()
        self.logger = logging.getLogger(f"{self.logname} {self.data}")
        self.logger.setLevel(self.loglevel)
        self.logformat = "[%(name)s %(levelname)s] %(message)s"
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter(self.logformat))
        self.logger.addHandler(handler)

    def _check_required_attributes(self):
        """Check if required attributes are present in data."""
        missing = [x for x in self.required_attrs
                    if not hasattr(self.data, x)]
        if missing:
            missing = ' '.join(missing)
            raise MissingAttributeError(
                f"Could not parse required attributes to use method: {missing}"
            )
