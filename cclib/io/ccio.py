# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
"""Tools for identifying, reading and writing files and streams."""

import io
import os
import typing
import warnings
from typing import Optional, Union

from cclib.attribute_parsers.data import ccData
from cclib.collection import ccCollection
from cclib.driver import ccDriver
from cclib.driver.ccdriver import triggers_on as triggers
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
from cclib.io.filereader import Reader
from cclib.tree import Tree

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
) -> Optional[ccCollection]:
    """Attempt to open and read computational chemistry data from a file.

    Inputs:
        source - a single logfile, a list of logfiles (for a single job),
                 an input stream, or an URL pointing to a log file.
        *args, **kwargs - arguments and keyword arguments passed to ccopen
    Returns:
        a ccCollection containing parsed cclib data attributes
    """
    inputobj = ccopen(source, *args, **kwargs)
    if inputobj is None:
        return None

    try:
        if isinstance(inputobj, Reader):
            tree = Tree()
            tree.add_root()
            collection = ccCollection(tree=tree)
            collection.parsed_data[0] = inputobj.parse()
            return collection

        if isinstance(inputobj, ccDriver):
            collection = inputobj.process_combinator()
            if not any(data.getattributes() for data in collection.parsed_data):
                return None
            return collection

        raise TypeError(f"Unsupported input object: {type(inputobj).__name__}")
    finally:
        if isinstance(inputobj, Reader):
            inputobj.inputfile.close()
        elif isinstance(inputobj, ccDriver):
            inputobj.fileHandler.close()


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
        A CJSON/XYZ reader or a ccDriver for computational output.
    """
    if source is None or source == "" or source == []:
        return None

    try:
        inputfile = source if isinstance(source, FileHandler) else FileHandler(source)
    except Exception:
        if not quiet:
            raise
        return None

    if cjson:
        return readerclasses["cjson"](inputfile, *args, **kwargs)

    if len(inputfile.filenames) == 1:
        extension = os.path.splitext(inputfile.filenames[0])[1][1:].lower()
        if extension in readerclasses:
            return readerclasses[extension](inputfile, *args, **kwargs)

    return ccDriver(inputfile, *args, **kwargs)


def _first_ccdata(ccobj: Union[ccData, ccCollection, ccDriver]) -> ccData:
    """Return the first data object represented by a supported v2 object."""
    if isinstance(ccobj, ccData):
        return ccobj

    if isinstance(ccobj, ccDriver):
        ccobj = ccobj.process_combinator()

    if isinstance(ccobj, ccCollection):
        if not ccobj.parsed_data:
            raise ValueError("Cannot select data from an empty ccCollection")
        return ccobj.parsed_data[0]

    raise ValueError(f"Unsupported object type: {type(ccobj).__name__}")


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
        ccobj - A ccDriver, ccCollection, or v2 ccData object
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

    jobfilename = None
    ccdata = _first_ccdata(ccobj)

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
        ccobjs - an iterable of ccDriver, ccCollection, or v2 ccData objects

    Returns:
        a pandas.DataFrame
    """
    _check_pandas(_has_pandas)
    logfiles = []
    for ccobj in ccobjs:
        jobfilename = None
        ccdata = _first_ccdata(ccobj)

        attributes = ccdata.getattributes()
        attributes.update({"jobfilename": jobfilename})

        logfiles.append(pd.Series(attributes))
    return pd.DataFrame(logfiles)


del find_package
