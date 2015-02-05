# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2009-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Tools for identifying and working with files and streams for any supported program"""


from __future__ import print_function

import os
import sys

from . import data

from . import logfileparser

from . import adfparser
from . import gamessparser
from . import gamessukparser
from . import gaussianparser
from . import jaguarparser
from . import molproparser
from . import nwchemparser
from . import orcaparser
from . import psiparser
from . import qchemparser

try:
    from ..bridge import cclib2openbabel
except ImportError:
    print("Could not import openbabel, fallback mechanism might not work.")


def ccread(source, *args, **kargs):
    """Attempt to open and read computational chemistry data from a file.

    If the file is not appropriate for cclib parsers, a fallback mechanism
    will try to recognize some common chemistry formats and read those using
    the appropriate bridge such as OpenBabel.    

    Inputs:
        source - a single logfile, a list of logfiles, or an input stream
    Returns:
        a ccData object containing cclib data attributes
    """
    log = ccopen(source, *args, **kargs)
    if log:
        if kargs['verbose']:
            print('Identified logfile to be in %s format' % log.logname)
        return log.parse()
    else:
        if kargs['verbose']:
            print('Attempting to use fallback mechanism to read file')
        return fallback(source)

def ccopen(source, *args, **kargs):
    """Guess the identity of a particular log file and return an instance of it.
    
    Inputs:
      source - a single logfile, a list of logfiles, or an input stream

    Returns:
      one of ADF, GAMESS, GAMESS UK, Gaussian, Jaguar, Molpro, NWChem, ORCA,
        Psi, QChem, or None (if it cannot figure it out or the file does not
        exist).
    """

    # Try to open the logfile(s), using openlogfile.
    if isinstance(source, str) or \
       isinstance(source, list) and all([isinstance(s, str) for s in source]):
        try:
            inputfile = logfileparser.openlogfile(source)
        except IOError as error:
            (errno, strerror) = error.args
            print("I/O error %s (%s): %s" % (errno, source, strerror))
            return None
        isstream = False
    elif hasattr(source, "read"):
        inputfile = source
        isstream = True
    else:
        raise ValueError

    # Read through the logfile(s) and search for a clue.
    filetype = None
    for line in inputfile:

        if line.find("Amsterdam Density Functional") >= 0:
            filetype = adfparser.ADF
            break

        # Don't break in this case as it may be a GAMESS-UK file.
        elif line.find("GAMESS") >= 0:
            filetype = gamessparser.GAMESS

        # This can break, since it is non-GAMESS-UK specific.
        elif line.find("GAMESS VERSION") >= 0:
            filetype = gamessparser.GAMESS
            break

        elif line.find("G A M E S S - U K") >= 0:
            filetype = gamessukparser.GAMESSUK
            break

        elif line.find("Gaussian, Inc.") >= 0:
            filetype = gaussianparser.Gaussian
            break

        elif line.find("Jaguar") >= 0:
            filetype = jaguarparser.Jaguar
            break

        elif line.find("PROGRAM SYSTEM MOLPRO") >= 0:
            filetype = molproparser.Molpro
            break

        # Molpro log files don't have the line above. Set this only if
        #   nothing else is detected, and notice it can be overwritten,
        #   since it does not break the loop.
        elif line[0:8] == "1PROGRAM" and not filetype:
            filetype = molproparser.Molpro

        elif line.find("Northwest Computational Chemistry Package") >= 0:
            filetype = nwchemparser.NWChem
            break

        elif line.find("O   R   C   A") >= 0:
            filetype = orcaparser.ORCA
            break

        elif line.find("PSI") >= 0 and line.find("Ab Initio Electronic Structure") >= 0:
            filetype = psiparser.Psi
            break

        elif line.find("A Quantum Leap Into The Future Of Chemistry") >= 0:
            filetype = qchemparser.QChem
            break

    # Need to close file before creating a instance.
    if not isstream:
        inputfile.close()
    
    # Return an instance of the logfile if one was chosen.
    if filetype:
        return filetype(source, *args, **kargs)

def fallback(source):
    """Attempt to read standard molecular formats using other libraries.

    Currently this will read XYZ files with OpenBabel, but this can easily
    be extended to other formats and libraries, too.
    """
    if isinstance(source, str):
        ext = os.path.splitext(source)[1][1:].lower()
        if 'cclib.bridge.cclib2openbabel' in sys.modules and ext in ('xyz'):
            return data.ccData(cclib2openbabel.readfile(source, ext))
