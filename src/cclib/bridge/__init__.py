# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

try:
    import openbabel
except Exception:
    pass
else:
    from .cclib2openbabel import makeopenbabel

try:
    import PyQuante
except ImportError:
    pass
else:
    from .cclib2pyquante import makepyquante

try:
    from .cclib2biopython import makebiopython
except ImportError:
    pass    
