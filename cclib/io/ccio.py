# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
"""Tools for identifying, reading and writing files and streams."""

import io
import logging
import os
import pathlib
import typing
import warnings
from typing import Optional, Union

# from cclib.io import (
#    cjsonreader,
#    cjsonwriter,
#    cmlwriter,
#    moldenwriter,
#    wfxwriter,
#    xyzreader,
#    xyzwriter,
# )
# from cclib.parser import data, logfileparser
from cclib.driver import ccDriver
from cclib.file_handler import FileHandler

# _has_cclib2openbabel = find_package("openbabel")
# if _has_cclib2openbabel:
#    from cclib.bridge import cclib2openbabel

# _has_pandas = find_package("pandas")
# if _has_pandas:
#    import pandas as pd

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

# readerclasses = {"cjson": cjsonreader.CJSON, "json": cjsonreader.CJSON, "xyz": xyzreader.XYZ}

# writerclasses = {
#     "cjson": cjsonwriter.CJSON,
#     "json": cjsonwriter.CJSON,
#     "cml": cmlwriter.CML,
#     "molden": moldenwriter.MOLDEN,
#     "wfx": wfxwriter.WFXWriter,
#     "xyz": xyzwriter.XYZ,
# }


class UnknownOutputFormatError(Exception):
    """Raised when an unknown output format is encountered."""


# def is_xyz(inputfile: FileWrapper) -> bool:
#     """Is the given inputfile actually an XYZ file?

#     The only way to determine this without reading the entire file is to
#     inspect the file extension.
#     """
#     return (
#         len(inputfile.filenames) == 1
#         and os.path.splitext(inputfile.filenames[0])[1].lower() == ".xyz"
#     )


# def guess_filetype(inputfile) -> Optional[logfileparser.Logfile]:
#     """Try to guess the filetype by searching for trigger strings."""
#     filetype = None
#     logger = logging.getLogger("cclib")
#     try:
#         if isinstance(inputfile, FileWrapper) and is_xyz(inputfile):
#             logger.info("Found XYZ file based on file extension")
#             return filetype
#         for line in inputfile:
#             for parser, phrases, do_break in triggers:
#                 if all([line.lower().find(p.lower()) >= 0 for p in phrases]):
#                     filetype = parser
#                     if do_break:
#                         return filetype
#     except Exception:
#         # guess_filetype() is expected to be quiet by default...
#         logger.error("Failed to determine log file type", exc_info=True)

#     return filetype


# def sort_turbomole_outputs(fileinputs):
#     """
#     Sorts a list of inputs (or list of log files) according to the order
#     required by the Turbomole parser for correct parsing. Unrecognised
#     files are appended to the end of the list in the same order they are
#     given.

#     This function has been deprecated as of version 1.8; use:
#     cclib.parser.turbomoleparser.Turbomole.sort_input() instead

#     Inputs:
#       filelist - a list of Turbomole log files needed to be parsed.
#     Returns:
#       sorted_list - a sorted list of Turbomole files needed for proper parsing.
#     """
#     warnings.warn(
#         "sort_turbomole_outputs() has been deprecated as of v1.8; use: "
#         + "cclib.parser.turbomoleparser.Turbomole.sort_input() instead"
#     )
#     return Turbomole.sort_input(fileinputs)


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


# def fallback(source):
#     """Attempt to read standard molecular formats using other libraries.

#     Currently this will read XYZ files with OpenBabel, but this can easily
#     be extended to other formats and libraries, too.
#     """

#     if isinstance(source, str):
#         ext = os.path.splitext(source)[1][1:].lower()
#         if _has_cclib2openbabel:
#             # From OB 3.0 onward, Pybel is contained inside the OB module.
#             try:
#                 import openbabel.pybel as pb
#             except:
#                 import pybel as pb
#             if ext in pb.informats:
#                 return cclib2openbabel.readfile(source, ext)
#         else:
#             # This should be a warning, but warnings are currently disabled by default.
#             logging.getLogger("cclib").error(
#                 "Could not import `openbabel`, fallback mechanism might not work."
#             )


# def ccwrite(
#     ccobj,
#     outputtype=None,
#     outputdest=None,
#     indices=None,
#     terse=False,
#     returnstr=False,
#     *args,
#     **kwargs,
# ):
#     """Write the parsed data from an outputfile to a standard chemical
#     representation.

#     Inputs:
#         ccobj - Either a job (from ccopen) or a data (from job.parse()) object
#         outputtype - The output format (should be a string)
#         outputdest - A filename or file object for writing
#         indices - One or more indices for extracting specific geometries/etc. (zero-based)
#         terse -  This option is currently limited to the cjson/json format. Whether to indent the cjson/json or not
#         returnstr - Whether or not to return a string representation.

#     The different writers may take additional arguments, which are
#     documented in their respective docstrings.

#     Returns:
#         the string representation of the chemical datatype
#           requested, or None.
#     """

#     # Determine the correct output format.
#     outputclass = _determine_output_format(outputtype, outputdest)

#     # Is ccobj an job object (unparsed), or is it a ccdata object (parsed)?
#     if isinstance(ccobj, logfileparser.Logfile):
#         jobfilename = ccobj.filename
#         ccdata = ccobj.parse()
#     elif isinstance(ccobj, data.ccData):
#         jobfilename = None
#         ccdata = ccobj
#     else:
#         raise ValueError

#     # If the logfile name has been passed in through kwargs (such as
#     # in the ccwrite script), make sure it has precedence.
#     if "jobfilename" in kwargs:
#         jobfilename = kwargs["jobfilename"]
#         # Avoid passing multiple times into the main call.
#         del kwargs["jobfilename"]

#     outputobj = outputclass(
#         ccdata, jobfilename=jobfilename, indices=indices, terse=terse, *args, **kwargs
#     )
#     output = outputobj.generate_repr()

#     # If outputdest isn't None, write the output to disk.
#     if outputdest is not None:
#         if isinstance(outputdest, str):
#             with open(outputdest, "w") as outputobj:
#                 outputobj.write(output)
#         elif isinstance(outputdest, io.IOBase):
#             outputdest.write(output)
#         else:
#             raise ValueError
#     # If outputdest is None, return a string representation of the output.
#     else:
#         return output

#     if returnstr:
#         return output


# def _determine_output_format(outputtype, outputdest):
#     """
#     Determine the correct output format.

#     Inputs:
#       outputtype - a string corresponding to the file type
#       outputdest - a filename string or file handle
#     Returns:
#       outputclass - the class corresponding to the correct output format
#     Raises:
#       UnknownOutputFormatError for unsupported file writer extensions
#     """

#     # Priority for determining the correct output format:
#     #  1. outputtype
#     #  2. outputdest

#     outputclass = None
#     # First check outputtype.
#     if isinstance(outputtype, str):
#         extension = outputtype.lower()
#         if extension in writerclasses:
#             outputclass = writerclasses[extension]
#         else:
#             raise UnknownOutputFormatError(extension)
#     else:
#         # Then checkout outputdest.
#         if isinstance(outputdest, str):
#             extension = os.path.splitext(outputdest)[1].lower()
#         elif isinstance(outputdest, io.IOBase):
#             extension = os.path.splitext(outputdest.name)[1].lower()
#         else:
#             raise UnknownOutputFormatError
#         if extension in writerclasses:
#             outputclass = writerclasses[extension]
#         else:
#             raise UnknownOutputFormatError(extension)

#     return outputclass


# def _check_pandas(found_pandas):
#     if not found_pandas:
#         raise ImportError("You must install `pandas` to use this function")


# def ccframe(ccobjs, *args, **kwargs):
#     """Returns a pandas.DataFrame of data attributes parsed by cclib from one
#     or more logfiles.

#     Inputs:
#         ccobjs - an iterable of either cclib jobs (from ccopen) or data (from
#         job.parse()) objects

#     Returns:
#         a pandas.DataFrame
#     """
#     _check_pandas(_has_pandas)
#     logfiles = []
#     for ccobj in ccobjs:
#         # Is ccobj an job object (unparsed), or is it a ccdata object (parsed)?
#         if isinstance(ccobj, logfileparser.Logfile):
#             jobfilename = ccobj.filename
#             ccdata = ccobj.parse()
#         elif isinstance(ccobj, data.ccData):
#             jobfilename = None
#             ccdata = ccobj
#         else:
#             raise ValueError

#         attributes = ccdata.getattributes()
#         attributes.update({"jobfilename": jobfilename})

#         logfiles.append(pd.Series(attributes))
#     return pd.DataFrame(logfiles)


# del find_package
