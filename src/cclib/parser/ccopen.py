# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2009-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

from __future__ import print_function

from . import logfileparser

from . import adfparser
from . import gamessparser
from . import gamessukparser
from . import gaussianparser
from . import jaguarparser
from . import molproparser
from . import orcaparser


def ccopen(source, *args, **kargs):
    """Guess the identity of a particular log file and return an instance of it.
    
    Inputs:
      source - a single logfile, a list of logfiles, or an input stream

    Returns:
      one of ADF, GAMESS, GAMESS UK, Gaussian, Jaguar, Molpro, ORCA, or
        None (if it cannot figure it out or the file does not exist).
    """

    filetype = None

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

        elif line.find("O   R   C   A") >= 0:
            filetype = orcaparser.ORCA
            break

    # Need to close file before creating a instance.
    if not isstream:
        inputfile.close()
    
    # Return an instance of the chosen class.
    try:
        return filetype(source, *args, **kargs)
    except TypeError:
        print("Log file type not identified.")
