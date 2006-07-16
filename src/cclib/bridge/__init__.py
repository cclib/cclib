"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision$"

try:
    import openbabel
    import pyopenbabel
except ImportError:
    pass
else:
    from cclib2openbabel import makeopenbabel

try:
    import PyQuante
except ImportError:
    pass
else:
    from cclib2pyquante import makepyquante

try:
    import Bio
except ImportError:
    pass
else:
    from cclib2biopython import makebiopython
