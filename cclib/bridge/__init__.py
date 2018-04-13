# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Facilities for moving parsed data to other cheminformatic libraries."""

try:
    import openbabel
except ImportError:
    pass
else:
    from cclib.bridge.cclib2openbabel import makeopenbabel

try:
    import PyQuante
except ImportError:
    pass
else:
    from cclib.bridge.cclib2pyquante import makepyquante

try:
    import Bio
except ImportError:
    pass
else:
    from cclib.bridge.cclib2biopython import makebiopython
