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

from cclib.io import (
    cjsonreader,
    cjsonwriter,
    cmlwriter,
    moldenwriter,
    wfxwriter,
    xyzreader,
    xyzwriter,
)
from cclib.parser import data, logfileparser
from cclib.parser.adfparser import ADF
from cclib.parser.daltonparser import DALTON
from cclib.parser.fchkparser import FChk
from cclib.parser.gamessdatparser import GAMESSDAT
from cclib.parser.gamessparser import GAMESS
from cclib.parser.gamessukparser import GAMESSUK
from cclib.parser.gaussianparser import Gaussian
from cclib.parser.jaguarparser import Jaguar
from cclib.parser.logfilewrapper import FileWrapper
from cclib.parser.molcasparser import Molcas
from cclib.parser.molproparser import Molpro
from cclib.parser.mopacparser import MOPAC
from cclib.parser.nboparser import NBO
from cclib.parser.nwchemparser import NWChem
from cclib.parser.orcaparser import ORCA
from cclib.parser.psi3parser import Psi3
from cclib.parser.psi4parser import Psi4
from cclib.parser.qchemparser import QChem
from cclib.parser.turbomoleparser import Turbomole
from cclib.parser.utils import find_package
from cclib.parser.xtbparser import XTB

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
triggers = [
    (ADF, ["Amsterdam Density Functional"], True),
    (DALTON, ["Dalton - An Electronic Structure Program"], True),
    (FChk, ["Number of atoms", "I"], True),
    (GAMESS, ["GAMESS"], False),
    (GAMESS, ["Firefly (PC GAMESS)"], True),
    (GAMESS, ["GAMESS VERSION"], True),
    (GAMESSUK, ["G A M E S S - U K"], True),
    (GAMESSDAT, ["$DATA"], True),
    (Gaussian, ["Gaussian, Inc."], True),
    (Jaguar, ["Jaguar"], True),
    (Molcas, ["MOLCAS"], True),
    (Molpro, ["PROGRAM SYSTEM MOLPRO"], True),
    (Molpro, ["1PROGRAM"], False),
    (MOPAC, ["MOPAC20"], True),
    (NBO, ["N A T U R A L   A T O M I C   O R B I T A L   A N D"], True),
    (NWChem, ["Northwest Computational Chemistry Package"], True),
    (ORCA, ["O   R   C   A"], True),
    (Psi3, ["PSI3: An Open-Source Ab Initio Electronic Structure Package"], True),
    (Psi4, ["Psi4: An Open-Source Ab Initio Electronic Structure Package"], True),
    (QChem, ["A Quantum Leap Into The Future Of Chemistry"], True),
    (Turbomole, ["TURBOMOLE"], True),
    (XTB, ["x T B"], True),
]

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


def is_xyz(inputfile: FileWrapper) -> bool:
    """Is the given inputfile actually an XYZ file?

    The only way to determine this without reading the entire file is to
    inspect the file extension.
    """
    return (
        len(inputfile.filenames) == 1
        and os.path.splitext(inputfile.filenames[0])[1].lower() == ".xyz"
    )


def guess_filetype(inputfile) -> Optional[logfileparser.Logfile]:
    """Try to guess the filetype by searching for trigger strings."""
    filetype = None
    logger = logging.getLogger("cclib")
    try:
        if isinstance(inputfile, FileWrapper) and is_xyz(inputfile):
            logger.info("Found XYZ file based on file extension")
            return filetype
        for line in inputfile:
            for parser, phrases, do_break in triggers:
                if all([line.lower().find(p.lower()) >= 0 for p in phrases]):
                    filetype = parser
                    if do_break:
                        return filetype
    except Exception:
        # guess_filetype() is expected to be quiet by default...
        logger.error("Failed to determine log file type", exc_info=True)

    return filetype


def sort_turbomole_outputs(fileinputs):
    """
    Sorts a list of inputs (or list of log files) according to the order
    required by the Turbomole parser for correct parsing. Unrecognised
    files are appended to the end of the list in the same order they are
    given.

    This function has been deprecated as of version 1.8; use:
    cclib.parser.turbomoleparser.Turbomole.sort_input() instead

    Inputs:
      filelist - a list of Turbomole log files needed to be parsed.
    Returns:
      sorted_list - a sorted list of Turbomole files needed for proper parsing.
    """
    warnings.warn(
        "sort_turbomole_outputs() has been deprecated as of v1.8; use: "
        + "cclib.parser.turbomoleparser.Turbomole.sort_input() instead"
    )
    return Turbomole.sort_input(fileinputs)


def ccread(
    source: Union[str, typing.IO, FileWrapper, typing.List[Union[str, typing.IO]]], *args, **kwargs
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
    log = None
    try:
        log = ccopen(source, *args, **kwargs)
        logger = logging.getLogger("cclib")
        if log:
            logger.info("Identified logfile to be in {} format".format(type(log).__name__))

            return log.parse()
        else:
            logger.info("Attempting to use fallback mechanism to read file")
            return fallback(source)

    finally:
        if log:
            log.inputfile.close()


def ccopen(
    source: Union[str, typing.IO, FileWrapper, typing.List[Union[str, typing.IO]]],
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
      one of ADF, DALTON, GAMESS, GAMESS UK, Gaussian, Jaguar,
      Molpro, MOPAC, NWChem, ORCA, Psi3, Psi/Psi4, QChem, CJSON or None
      (if it cannot figure it out or the file does not exist).
    """
    if not isinstance(source, list):
        source = [source]

    inputfile = None

    logger = logging.getLogger("cclib")

    try:
        # Wrap our input with custom file object.
        inputfile = FileWrapper(*source)

        if cjson:
            filetype = readerclasses["cjson"]

        else:
            # Try and guess the parser we need.
            filetype = guess_filetype(inputfile)

            # Reset our position back to 0.
            inputfile.reset()

        # If the input file isn't a standard compchem log file, try one of
        # the readers, falling back to Open Babel.
        if not filetype:
            # TODO: This assumes we only got a single file...
            filename = list(inputfile.filenames)[0]
            ext = pathlib.Path(filename).name[1:].lower()

            for extension in readerclasses:
                if ext == extension:
                    filetype = readerclasses[extension]

        # Proceed to return an instance of the logfile parser only if the filetype
        # could be guessed.
        if filetype:
            return filetype(inputfile, *args, **kwargs)

        elif inputfile is not None:
            inputfile.close()
            # Stop us closing twice in the except block.
            inputfile = None

        logger.warning(
            "Unable to determine the type of logfile %s, try the fallback mechanism", source
        )

    except Exception:
        if inputfile is not None:
            inputfile.close()

        if not quiet:
            raise

        # We're going to swallow this exception if quiet is True.
        # This can hide a lot of errors, so we'll make sure to log it.
        logger.error("Failed to open logfile", exc_info=True)


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
                try:
                    import pybel as pb
                except:
                    return
            if ext in pb.informats:
                return cclib2openbabel.readfile(source, ext)
        else:
            # This should be a warning, but warnings are currently disabled by default.
            logging.getLogger("cclib").error(
                "Could not import `openbabel`, fallback mechanism might not work."
            )


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
        attributes.update({"jobfilename": jobfilename})

        logfiles.append(pd.Series(attributes))
    return pd.DataFrame(logfiles)


del find_package
