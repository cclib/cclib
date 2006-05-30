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
