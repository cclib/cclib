# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in horton (http://theochem.github.io/horton)."""

import numpy
from cclib.parser.data import ccData
from cclib.parser.utils import find_package

# First check horton version
_old_horton = False
_found_horton = find_package('horton')
_found_iodata = find_package('iodata')

if _found_horton: # Detect whether horton 2 is present or not
    try: # Older horton packages do not have __version__, causing exceptions.
        from horton import __version__
    except:
        _old_horton = True
    else:
        if __version__[0] == "2":
            from horton.io.iodata import IOData

if _found_iodata: # Detect whether iodata (part of horton 3) is present or not; Horton 3 is divided into smaller (sub)packages that each take different functionalities.
    from iodata import IOData
    from iodata.orbitals import MolecularOrbitals
    
def check_horton():
    if _old_horton:
        raise ImportError("You must have at least version 2 of `horton` to use this function.")
    elif not _found_horton and not _found_iodata:
        raise ImportError("You must install `horton` to use this function.")
    if (_found_iodata):
        return 3
    elif (_found_horton):
        return 2

def makehorton(ccdat):
    """ Create horton IOData object from ccData object """
    
    hortonver = check_horton()
    attributes = {}
    
    if (hortonver == 2):
        pass # Populate with bridge for horton 2 in later PR
    elif (hortonver == 3):
        pass # Populate with bridge for horton 3 in later PR
    
    return IOData(**attributes) # Pass collected attributes into IOData constructor
    
def makecclib(iodat):
    """ Create cclib ccData object from horton IOData object """
    
    hortonver = check_horton()
    attributes = {}
    
    if(hortonver == 2):
        pass
    elif (hortonver == 3):
        pass
    
    return ccData(attributes)
    