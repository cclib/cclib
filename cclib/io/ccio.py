# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
"""Tools for identifying, reading and writing files and streams."""

import io
import logging
import os
import typing
import warnings
from typing import Optional, Union

from cclib.attribute_parsers.data import ccData
from cclib.driver.ccdriver import triggers_on as triggers
from cclib.driver import ccDriver
from cclib.file_handler import FileHandler
from cclib.file_handler.utils import find_package
from cclib.io import (
    cjsonreader,
    cjsonwriter,
    cmlwriter,
    moldenwriter,
    wfxwriter,
    xyzreader,
    xyzwriter,
)

FileWrapper = FileHandler

_has_cclib2openbabel = find_package("openbabel")
if _has_cclib2openbabel:
    from cclib.bridge import cclib2openbabel

_has_pandas = find_package("pandas")
if _has_pandas:
    import pandas as pd

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

readerclasses = {"cjson": cjsonreader.CJSON, "json": cjsonreader.CJSON, "xyz": xyzreader.XYZ}

writerclasses = {
    "cjson": cjsonwriter.CJSON,
    "json": cjsonwriter.CJSON,
    "cml": cmlwriter.CML,
    "molden": moldenwriter.MOLDEN,
    "wfx": wfxwriter.WFXWriter,
    "xyz": xyzwriter.XYZ,
}


class UnknownOutputFormatError(Exception):
    """Raised when an unknown output format is encountered."""

def ccread(
    source: Union[str, typing.IO, FileHandler, typing.List[Union[str, typing.IO]]], *args, **kwargs
):
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
    if not isinstance(source, list):
        source = [source]

    a = ccDriver(source, **kwargs)
    a.process_combinator()
    return a._ccCollection._parsed_data


def ccopen(
    source: Union[str, typing.IO, FileHandler, typing.List[Union[str, typing.IO]]],
    *args,
    quiet: bool = False,
    cjson: bool = False,
    **kwargs,
):
    """Guess the identity of a particular log file and return an instance of it.

    Inputs:
        source - a single logfile, a list of logfiles (for a single job),
                 an input stream, or an URL pointing to a log file.
        *args, **kwargs - arguments and keyword arguments passed to filetype

    Returns:
        ccCollection: by default a single point combinator, eventually dynamically determined.
    """
    if not isinstance(source, list):
        source = [source]

    ccdriver_inst = ccDriver(source)
    return ccdriver_inst


def ccwrite(
    ccobj,
    outputtype=None,
    outputdest=None,
    indices=None,
    terse=False,
    returnstr=False,
    *args,
    **kwargs,
):
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

    # Is ccobj an unparsed driver or a data object (parsed)?
    if isinstance(ccobj, ccData):
        jobfilename = None
        ccdata = ccobj
    elif hasattr(ccobj, "process_combinator"):
        jobfilename = None
        parsed_data = ccobj.process_combinator().parsed_data
        ccdata = parsed_data[0] if parsed_data else None
    else:
        raise ValueError
    if ccdata is None:
        raise ValueError

    # If the logfile name has been passed in through kwargs (such as
    # in the ccwrite script), make sure it has precedence.
    if "jobfilename" in kwargs:
        jobfilename = kwargs["jobfilename"]
        # Avoid passing multiple times into the main call.
        del kwargs["jobfilename"]

    outputobj = outputclass(
        ccdata, jobfilename=jobfilename, indices=indices, terse=terse, *args, **kwargs
    )
    output = outputobj.generate_repr()

    # If outputdest isn't None, write the output to disk.
    if outputdest is not None:
        if isinstance(outputdest, str):
            with open(outputdest, "w") as outputobj:
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
        if isinstance(ccobj, ccData):
            jobfilename = None
            ccdata = ccobj
        elif hasattr(ccobj, "process_combinator"):
            jobfilename = None
            parsed_data = ccobj.process_combinator().parsed_data
            ccdata = parsed_data[0] if parsed_data else None
        else:
            raise ValueError

        if ccdata is None:
            raise ValueError

        attributes = ccdata.getattributes()
        attributes.update({"jobfilename": jobfilename})

        logfiles.append(pd.Series(attributes))
    return pd.DataFrame(logfiles)


del find_package
