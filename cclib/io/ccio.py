# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
"""Tools for identifying, reading and writing files and streams."""

import atexit
import io
import os
import sys
import re
import pathlib
from tempfile import NamedTemporaryFile
from urllib.request import urlopen
from urllib.error import URLError

from cclib.parser import data
from cclib.parser import logfileparser
from cclib.parser.utils import find_package

from cclib.parser.adfparser import ADF
from cclib.parser.daltonparser import DALTON
from cclib.parser.fchkparser import FChk
from cclib.parser.gamessparser import GAMESS
from cclib.parser.gamessukparser import GAMESSUK
from cclib.parser.gaussianparser import Gaussian
from cclib.parser.jaguarparser import Jaguar
from cclib.parser.molcasparser import Molcas
from cclib.parser.molproparser import Molpro
from cclib.parser.mopacparser import MOPAC
from cclib.parser.nwchemparser import NWChem
from cclib.parser.orcaparser import ORCA
from cclib.parser.psi3parser import Psi3
from cclib.parser.psi4parser import Psi4
from cclib.parser.qchemparser import QChem
from cclib.parser.turbomoleparser import Turbomole

from cclib.io import cjsonreader
from cclib.io import cjsonwriter
from cclib.io import cmlwriter
from cclib.io import moldenwriter
from cclib.io import wfxwriter
from cclib.io import xyzreader
from cclib.io import xyzwriter

_has_cclib2openbabel = find_package("openbabel")
if _has_cclib2openbabel:
    from cclib.bridge import cclib2openbabel

_has_pandas = find_package("pandas")
if _has_pandas:
    import pandas as pd

# Regular expression for validating URLs
URL_PATTERN = re.compile(

    r'^(?:http|ftp)s?://'  # http:// or https://
    r'(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|'  # domain...
    r'localhost|'  # localhost...
    r'\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})'  # ...or ip
    r'(?::\d+)?'  # optional port
    r'(?:/?|[/?]\S+)$', re.IGNORECASE

)

# Parser choice is triggered by certain phrases occurring the logfile. Where these
# strings are unique, we can set the parser and break. In other cases, the situation
# is a little but more complicated. Here are the exceptions:
#   1. The GAMESS trigger also works for GAMESS-UK files, so we can't break
#      after finding GAMESS in case the more specific phrase is found.
#   2. Molpro log files don't have the program header, but always contain
#      the generic string 1PROGRAM, so don't break here either to be cautious.
#   3. "MOPAC" is used in some packages like GAMESS, so match MOPAC20##
#
# The triggers are defined by the tuples in the list below like so:
#   (parser, phrases, flag whether we should break)
triggers = [

    (ADF,       ["Amsterdam Density Functional"],                   True),
    (DALTON,    ["Dalton - An Electronic Structure Program"],       True),
    (FChk,      ["Number of atoms", "I"],                           True),
    (GAMESS,    ["GAMESS"],                                         False),
    (GAMESS,    ["GAMESS VERSION"],                                 True),
    (GAMESSUK,  ["G A M E S S - U K"],                              True),
    (Gaussian,  ["Gaussian, Inc."],                                 True),
    (Jaguar,    ["Jaguar"],                                         True),
    (Molcas,    ["MOLCAS"],                                         True),
    (Molpro,    ["PROGRAM SYSTEM MOLPRO"],                          True),
    (Molpro,    ["1PROGRAM"],                                       False),
    (MOPAC,     ["MOPAC20"],                                        True),
    (NWChem,    ["Northwest Computational Chemistry Package"],      True),
    (ORCA,      ["O   R   C   A"],                                  True),
    (Psi3,      ["PSI3: An Open-Source Ab Initio Electronic Structure Package"],          True),
    (Psi4,      ["Psi4: An Open-Source Ab Initio Electronic Structure Package"],          True),
    (QChem,     ["A Quantum Leap Into The Future Of Chemistry"],    True),
    (Turbomole, ["TURBOMOLE"],                                      True),

]

readerclasses = {
    'cjson': cjsonreader.CJSON,
    'json': cjsonreader.CJSON,
    'xyz': xyzreader.XYZ,
}

writerclasses = {
    'cjson': cjsonwriter.CJSON,
    'json': cjsonwriter.CJSON,
    'cml': cmlwriter.CML,
    'molden': moldenwriter.MOLDEN,
    'wfx': wfxwriter.WFXWriter,
    'xyz': xyzwriter.XYZ,
}


class UnknownOutputFormatError(Exception):
    """Raised when an unknown output format is encountered."""


def guess_filetype(inputfile):
    """Try to guess the filetype by searching for trigger strings."""
    if not inputfile:
        return None

    filetype = None
    if isinstance(inputfile, str):
        for line in inputfile:
            for parser, phrases, do_break in triggers:
                if all([line.lower().find(p.lower()) >= 0 for p in phrases]):
                    filetype = parser
                    if do_break:
                        return filetype
    else:
        for fname in inputfile:
            for line in inputfile:
                for parser, phrases, do_break in triggers:
                    if all([line.lower().find(p.lower()) >= 0 for p in phrases]):
                        filetype = parser
                        if do_break:
                            return filetype
    return filetype


def ccread(source, *args, **kwargs):
    """Attempt to open and read computational chemistry data from a file.

    If the file is not appropriate for cclib parsers, a fallback mechanism
    will try to recognize some common chemistry formats and read those using
    the appropriate bridge such as Open Babel.

    Inputs:
        source - a single logfile, a list of logfiles (for a single job),
                 an input stream, or an URL pointing to a log file.
        *args, **kwargs - arguments and keyword arguments passed to ccopen
    Returns:
        a ccData object containing cclib data attributes
    """

    log = ccopen(source, *args, **kwargs)
    if log:
        if kwargs.get("verbose", None):
            print(f"Identified logfile to be in {log.logname} format")
        # If the input file is a CJSON file and not a standard compchemlog file
        cjson_as_input = kwargs.get("cjson", False)
        if cjson_as_input:
            return log.read_cjson()
        else:
            return log.parse()
    else:
        if kwargs.get('verbose', None):
            print('Attempting to use fallback mechanism to read file')
        return fallback(source)


def ccopen(source, *args, **kwargs):
    """Guess the identity of a particular log file and return an instance of it.

    Inputs:
        source - a single logfile, a list of logfiles (for a single job),
                 an input stream, or an URL pointing to a log file.
        *args, **kwargs - arguments and keyword arguments passed to filetype

    Returns:
      one of ADF, DALTON, GAMESS, GAMESS UK, Gaussian, Jaguar,
      Molpro, MOPAC, NWChem, ORCA, Psi3, Psi/Psi4, QChem, CJSON or None
      (if it cannot figure it out or the file does not exist).
    """
    inputfile = None
    is_stream = False

    # Check if source is a link or contains links. Retrieve their content.
    # Try to open the logfile(s), using openlogfile, if the source is a string (filename)
    # or list of filenames. If it can be read, assume it is an open file object/stream.
    if isinstance(source, pathlib.PurePath):
        source = str(source)
    if isinstance(source, pathlib.PurePath)\
            and all([isinstance(s, pathlib.PurePath) for s in source]):
        source = [str(item) for item in source]
    is_string = isinstance(source, str)
    is_url = True if is_string and URL_PATTERN.match(source) else False
    is_listofstrings = isinstance(source, list) and all([isinstance(s, str) for s in source])
    if is_string or is_listofstrings:
        # Process links from list (download contents into temporary location)
        if is_listofstrings:
            filelist = []
            for filename in source:
                if not URL_PATTERN.match(filename):
                    filelist.append(filename)
                else:
                    try:
                        response = urlopen(filename)
                        tfile = NamedTemporaryFile(delete=False)
                        tfile.write(response.read())
                        # Close the file because Windows won't let open it second time
                        tfile.close()
                        filelist.append(tfile.name)
                        # Delete temporary file when the program finishes
                        atexit.register(os.remove, tfile.name)
                    except (ValueError, URLError) as error:
                        if not kwargs.get('quiet', False):
                            (errno, strerror) = error.args
                        return None
            source = filelist

        if not is_url:
            try:
                inputfile = logfileparser.openlogfile(source)
            except IOError as error:
                if not kwargs.get('quiet', False):
                    (errno, strerror) = error.args
                return None
        else:
            try:
                response = urlopen(source)
                is_stream = True

                # Retrieve filename from URL if possible
                filename = re.findall(r"\w+\.\w+", source.split('/')[-1])
                filename = filename[0] if filename else ""

                inputfile = logfileparser.openlogfile(filename, object=response.read())
            except (ValueError, URLError) as error:
                if not kwargs.get('quiet', False):
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
    filetype = guess_filetype(inputfile)

    # If the input file isn't a standard compchem log file, try one of
    # the readers, falling back to Open Babel.
    if not filetype:
        if kwargs.get("cjson"):
            filetype = readerclasses['cjson']
        elif source and not is_stream:
            ext = os.path.splitext(source)[1][1:].lower()
            for extension in readerclasses:
                if ext == extension:
                    filetype = readerclasses[extension]

    # Proceed to return an instance of the logfile parser only if the filetype
    # could be guessed. Need to make sure the input file is closed before creating
    # an instance, because parsers will handle opening/closing on their own.
    if filetype:
        # We're going to close and reopen below anyway, so this is just to avoid
        # the missing seek method for fileinput.FileInput. In the long run
        # we need to refactor to support for various input types in a more
        # centralized fashion.
        if is_listofstrings:
            pass
        else:
            inputfile.seek(0, 0)
        if not is_stream:
            if is_listofstrings:
                if filetype == Turbomole:
                    source = sort_turbomole_outputs(source)
            inputfile.close()
            return filetype(source, *args, **kwargs)
        return filetype(inputfile, *args, **kwargs)


def fallback(source):
    """Attempt to read standard molecular formats using other libraries.

    Currently this will read XYZ files with OpenBabel, but this can easily
    be extended to other formats and libraries, too.
    """

    if isinstance(source, str):
        ext = os.path.splitext(source)[1][1:].lower()
        if _has_cclib2openbabel:
            # From OB 3.0 onward, Pybel is contained inside the OB module.
            try:
                import openbabel.pybel as pb
            except:
                import pybel as pb
            if ext in pb.informats:
                return cclib2openbabel.readfile(source, ext)
        else:
            print("Could not import `openbabel`, fallback mechanism might not work.")


def ccwrite(ccobj, outputtype=None, outputdest=None,
            indices=None, terse=False, returnstr=False,
            *args, **kwargs):
    """Write the parsed data from an outputfile to a standard chemical
    representation.

    Inputs:
        ccobj - Either a job (from ccopen) or a data (from job.parse()) object
        outputtype - The output format (should be a string)
        outputdest - A filename or file object for writing
        indices - One or more indices for extracting specific geometries/etc. (zero-based)
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
    if 'jobfilename' in kwargs:
        jobfilename = kwargs['jobfilename']
        # Avoid passing multiple times into the main call.
        del kwargs['jobfilename']

    outputobj = outputclass(ccdata, jobfilename=jobfilename,
                            indices=indices, terse=terse,
                            *args, **kwargs)
    output = outputobj.generate_repr()

    # If outputdest isn't None, write the output to disk.
    if outputdest is not None:
        if isinstance(outputdest, str):
            with open(outputdest, 'w') as outputobj:
                outputobj.write(output)
        elif isinstance(outputdest, io.IOBase):
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
      outputdest - a filename string or file handle
    Returns:
      outputclass - the class corresponding to the correct output format
    Raises:
      UnknownOutputFormatError for unsupported file writer extensions
    """

    # Priority for determining the correct output format:
    #  1. outputtype
    #  2. outputdest

    outputclass = None
    # First check outputtype.
    if isinstance(outputtype, str):
        extension = outputtype.lower()
        if extension in writerclasses:
            outputclass = writerclasses[extension]
        else:
            raise UnknownOutputFormatError(extension)
    else:
        # Then checkout outputdest.
        if isinstance(outputdest, str):
            extension = os.path.splitext(outputdest)[1].lower()
        elif isinstance(outputdest, io.IOBase):
            extension = os.path.splitext(outputdest.name)[1].lower()
        else:
            raise UnknownOutputFormatError
        if extension in writerclasses:
            outputclass = writerclasses[extension]
        else:
            raise UnknownOutputFormatError(extension)

    return outputclass

def path_leaf(path):
    """
    Splits the path to give the filename. Works irrespective of '\'
    or '/' appearing in the path and also with path ending with '/' or '\'.

    Inputs:
      path - a string path of a logfile.
    Returns:
      tail - 'directory/subdirectory/logfilename' will return 'logfilename'.
      ntpath.basename(head) - 'directory/subdirectory/logfilename/' will return 'logfilename'.
    """
    head, tail = os.path.split(path)
    return tail or os.path.basename(head)

def sort_turbomole_outputs(filelist):
    """
    Sorts a list of inputs (or list of log files) according to the order
    defined below. Just appends the unknown files in the end of the sorted list.

    Inputs:
      filelist - a list of Turbomole log files needed to be parsed.
    Returns:
      sorted_list - a sorted list of Turbomole files needed for proper parsing.
    """
    sorting_order = {
        'basis' : 0,
        'control' : 1,
        'mos' : 2,
        'alpha' : 3,
        'beta' : 4,
        'job.last' : 5,
        'coord' : 6,
        'gradient' : 7,
        'aoforce' : 8,
    }
    
    known_files = []
    unknown_files = []
    sorted_list = []
    for fname in filelist:
        filename = path_leaf(fname)
        if filename in sorting_order:
            known_files.append([fname, sorting_order[filename]])
        elif re.match(r"^job\.[0-9]+$", filename):
            # Calling 'jobex -keep' will also write job.n files, where n ranges from 0 to inf.
            # Numbered job files are inserted before job.last.
            job_number = int(filename[4:]) +1
            job_order = float(f"{sorting_order['job.last'] - 1}.{job_number}")
            known_files.append([fname, job_order])
        else:
            unknown_files.append(fname)
    for i in sorted(known_files, key=lambda x: x[1]):
        sorted_list.append(i[0])
    if unknown_files:
        sorted_list.extend(unknown_files)
    return sorted_list


def _check_pandas(found_pandas):
    if not found_pandas:
        raise ImportError("You must install `pandas` to use this function")


def ccframe(ccobjs, *args, **kwargs):
    """Returns a pandas.DataFrame of data attributes parsed by cclib from one
    or more logfiles.

    Inputs:
        ccobjs - an iterable of either cclib jobs (from ccopen) or data (from
        job.parse()) objects

    Returns:
        a pandas.DataFrame
    """
    _check_pandas(_has_pandas)
    logfiles = []
    for ccobj in ccobjs:
        # Is ccobj an job object (unparsed), or is it a ccdata object (parsed)?
        if isinstance(ccobj, logfileparser.Logfile):
            jobfilename = ccobj.filename
            ccdata = ccobj.parse()
        elif isinstance(ccobj, data.ccData):
            jobfilename = None
            ccdata = ccobj
        else:
            raise ValueError

        attributes = ccdata.getattributes()
        attributes.update({
            'jobfilename': jobfilename
        })

        logfiles.append(pd.Series(attributes))
    return pd.DataFrame(logfiles)


del find_package

