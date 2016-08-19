# -*- coding: utf-8 -*-
#
# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2009-2016, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Tools for identifying, reading and writing files and streams."""


from __future__ import print_function

import io
import os
import sys

# Python 2->3 changes the default file object hierarchy.
if sys.version_info[0] == 2:
    fileclass = file
else:
    import io
    fileclass = io.IOBase

from ..parser import logfileparser
from ..parser import data

from ..parser.adfparser import ADF
from ..parser.daltonparser import DALTON
from ..parser.gamessparser import GAMESS
from ..parser.gamessukparser import GAMESSUK
from ..parser.gaussianparser import Gaussian
from ..parser.jaguarparser import Jaguar
from ..parser.molproparser import Molpro
from ..parser.nwchemparser import NWChem
from ..parser.orcaparser import ORCA
from ..parser.psiparser import Psi
from ..parser.qchemparser import QChem

from . import cjsonreader
from . import cjsonwriter
from . import cmlwriter
from . import xyzwriter

try:
    from ..bridge import cclib2openbabel
    _has_cclib2openbabel = True
except ImportError:
    _has_cclib2openbabel = False


# Parser choice is triggered by certain phrases occuring the logfile. Where these
# strings are unique, we can set the parser and break. In other cases, the situation
# is a little but more complicated. Here are the exceptions:
#   1. The GAMESS trigger also works for GAMESS-UK files, so we can't break
#      after finding GAMESS in case the more specific phrase is found.
#   2. Molpro log files don't have the program header, but always contain
#      the generic string 1PROGRAM, so don't break here either to be cautious.
#   3. The Psi header has two different strings with some variation
#
# The triggers are defined by the tuples in the list below like so:
#   (parser, phrases, flag whether we should break)
triggers = [

    (ADF,       ["Amsterdam Density Functional"],                   True),
    (DALTON,    ["Dalton - An Electronic Structure Program"],       True),
    (GAMESS,    ["GAMESS"],                                         False),
    (GAMESS,    ["GAMESS VERSION"],                                 True),
    (GAMESSUK,  ["G A M E S S - U K"],                              True),
    (Gaussian,  ["Gaussian, Inc."],                                 True),
    (Jaguar,    ["Jaguar"],                                         True),
    (Molpro,    ["PROGRAM SYSTEM MOLPRO"],                          True),
    (Molpro,    ["1PROGRAM"],                                       False),
    (NWChem,    ["Northwest Computational Chemistry Package"],      True),
    (ORCA,      ["O   R   C   A"],                                  True),
    (Psi,       ["PSI", "Ab Initio Electronic Structure"],          True),
    (QChem,     ["A Quantum Leap Into The Future Of Chemistry"],    True),

]


def guess_filetype(inputfile):
    """Try to guess the filetype by searching for trigger strings."""
    if not inputfile:
        return None

    filetype = None
    for line in inputfile:
        for parser, phrases, do_break in triggers:
            if all([line.lower().find(p.lower()) >= 0 for p in phrases]):
                filetype = parser
                if do_break:
                    return filetype
    return filetype


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
        if kargs.get('verbose', None):
            print('Identified logfile to be in %s format' % log.logname)
        # If the input file is a CJSON file and not a standard compchemlog file
        cjson_as_input = kargs.get("cjson", False)
        if cjson_as_input:
            return log.read_cjson()
        else:
            return log.parse()
    else:
        if kargs.get('verbose', None):
            print('Attempting to use fallback mechanism to read file')
        return fallback(source)


def ccopen(source, *args, **kargs):
    """Guess the identity of a particular log file and return an instance of it.

    Inputs:
      source - a single logfile, a list of logfiles, or an input stream

    Returns:
      one of ADF, DALTON, GAMESS, GAMESS UK, Gaussian, Jaguar, Molpro, NWChem, ORCA,
        Psi, QChem, CJSON or None (if it cannot figure it out or the file does not
        exist).
    """

    inputfile = None
    is_stream = False

    # Try to open the logfile(s), using openlogfile, if the source is a string (filename)
    # or list of filenames. If it can be read, assume it is an open file object/stream.
    is_string = isinstance(source, str)
    is_listofstrings = isinstance(source, list) and all([isinstance(s, str) for s in source])
    if is_string or is_listofstrings:
        try:
            inputfile = logfileparser.openlogfile(source)
        except IOError as error:
            if not kargs.get('quiet', False):
                (errno, strerror) = error.args
            return None
    elif hasattr(source, "read"):
        inputfile = source
        is_stream = True

    # Streams are tricky since they don't have seek methods or seek won't work
    # by design even if it is present. We solve this now by reading in the
    # entire stream and using a StringIO buffer for parsing. This might be
    # problematic for very large streams. Slow streams might also be an issue if
    # the parsing is not instantaneous, but we'll deal with such edge cases
    # as they arise. Ideally, in the future we'll create a class dedicated to
    # dealing with these issues, supporting both files and streams.
    if is_stream:
        try:
            inputfile.seek(0, 0)
        except (AttributeError, IOError):
            contents = inputfile.read()
            try:
              inputfile = io.StringIO(contents)
            except:
              inputfile = io.StringIO(unicode(contents))
            inputfile.seek(0, 0)

    # Proceed to return an instance of the logfile parser only if the filetype
    # could be guessed. Need to make sure the input file is closed before creating
    # an instance, because parsers will handle opening/closing on their own.
    # If the input file is a CJSON file and not a standard compchemlog file, don't
    # guess the file.
    if kargs.get("cjson", False):
        filetype = cjsonreader.CJSON
    else:
        filetype = guess_filetype(inputfile)

    # Proceed to return an instance of the logfile parser only if the filetype
    # could be guessed. Need to make sure the input file is closed before creating
    # an instance, because parsers will handle opening/closing on their own.
    if filetype:
        inputfile.seek(0, 0)
        if not is_stream:
            inputfile.close()
            return filetype(source, *args, **kargs)
        return filetype(inputfile, *args, **kargs)


def fallback(source):
    """Attempt to read standard molecular formats using other libraries.

    Currently this will read XYZ files with OpenBabel, but this can easily
    be extended to other formats and libraries, too.
    """

    if isinstance(source, str):
        ext = os.path.splitext(source)[1][1:].lower()
        if _has_cclib2openbabel:
            if ext in ('xyz', ):
                return cclib2openbabel.readfile(source, ext)
        else:
            print("Could not import openbabel, fallback mechanism might not work.")


def ccwrite(ccobj, outputtype=None, outputdest=None, terse=False , returnstr=False,
            *args, **kwargs):
    """Write the parsed data from an outputfile to a standard chemical
    representation.

    Inputs:
        ccobj - Either a job (from ccopen) or a data (from job.parse()) object
        outputtype - The output format (should be one of 'cjson', 'cml', 'xyz')
        outputdest - A filename or file object for writing
        terse -  This option is currently limited to the cjson/json format. Whether to indent the cjson/json or not
        returnstr - Whether or not to return a string representation.

    The different writers may take additional arguments, which are
    documented in their respective docstrings.

    Returns:
        the string representation of the chemical datatype
          requested, or None.
    """

    # Determine the correct output format.
    outputclass = _determine_output_format(outputtype, outputdest)

    # Is ccobj an job object (unparsed), or is it a ccdata object (parsed)?
    if isinstance(ccobj, logfileparser.Logfile):
        jobfilename = ccobj.filename
        ccdata = ccobj.parse()
    elif isinstance(ccobj, data.ccData):
        jobfilename = None
        ccdata = ccobj
    else:
        raise ValueError

    # If the logfile name has been passed in through kwargs (such as
    # in the ccwrite script), make sure it has precedence.
    if 'jobfilename' in kwargs.keys():
        jobfilename = kwargs['jobfilename']
        # Avoid passing multiple times into the main call.
        del kwargs['jobfilename']

    outputobj = outputclass(ccdata, jobfilename=jobfilename, terse=terse, *args, **kwargs)
    output = outputobj.generate_repr()

    # If outputdest isn't None, write the output to disk.
    if outputdest is not None:
        if isinstance(outputdest, str):
            with open(outputdest, 'w') as outputobj:
                outputobj.write(output)
        elif isinstance(outputdest, fileclass):
            outputdest.write(output)
        else:
            raise ValueError
    # If outputdest is None, return a string representation of the output.
    else:
        return output

    if returnstr:
        return output


def _determine_output_format(outputtype, outputdest):
    """
    Determine the correct output format.

    Inputs:
      outputtype - a string corresponding to the file type
        (one of cjson/json, cml, xyz)
      outputdest - a filename string or file handle
    Returns:
      outputclass - the class corresponding to the correct output format
    """

    # Priority for determining the correct output format:
    #  1. outputtype
    #  2. outputdest

    # First check outputtype.
    if isinstance(outputtype, str):
        if outputtype.lower() in ('cjson', 'json'):
            outputclass = cjsonwriter.CJSON
        elif outputtype.lower() == 'cml':
            outputclass = cmlwriter.CML
        elif outputtype.lower() == 'xyz':
            outputclass = xyzwriter.XYZ
    else:
        # Then checkout outputdest.
        if isinstance(outputdest, str):
            extension = os.path.splitext(outputdest)[1]
        elif isinstance(outputdest, fileclass):
            extension = os.path.splitext(outputdest.name)[1]
        else:
            raise ValueError
        if extension.lower() in ('.cjson', '.json'):
            outputclass = cjsonwriter.CJSON
        elif extension.lower() == '.cml':
            outputclass = cmlwriter.CML
        elif extension.lower() == '.xyz':
            outputclass = xyzwriter.XYZ
        else:
            raise ValueError

    return outputclass
