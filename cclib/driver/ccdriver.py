# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
"""Tools for identifying, reading and writing files and streams."""

import io
import logging
import os
import pathlib
import sys
import typing
import warnings
from enum import Enum, auto
from typing import TYPE_CHECKING, List, Optional, Union

from cclib.collection import ccCollection
from cclib.combinator import AutoCombinator
from cclib.file_handler import FileHandler
from cclib.tree import Tree

if TYPE_CHECKING:
    from cclib.combinator import Combinator

# from cclib.parser.utils import find_package

# todo from cclib.io import cjsonreader
# todo from cclib.io import cjsonwriter
# todo from cclib.io import cmlwriter
# todo from cclib.io import moldenwriter
# todo from cclib.io import wfxwriter
# todo from cclib.io import xyzreader
# todo from cclib.io import xyzwriter

from dataclasses import dataclass


class TriggerType(Enum):
    ON = auto()
    OFF = auto()


# TODO possibly replace with attrs so we can always use kw_only
@dataclass
class Trigger:
    program: str
    trigger_strings: List[str]
    do_break: bool
    type: TriggerType

    def matches(self, line: str) -> bool:
        return all([line.lower().find(phrase.lower()) >= 0 for phrase in self.trigger_strings])


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
TRIGGERS = [
    Trigger(
        program="adf",
        trigger_strings=["Amsterdam Density Functional"],
        do_break=True,
        type=TriggerType.ON,
    ),
    Trigger(
        program="dalton",
        trigger_strings=["Dalton - An Electronic Structure Program"],
        do_break=True,
        type=TriggerType.ON,
    ),
    Trigger(
        program="fchk", trigger_strings=["Number of atoms", "I"], do_break=True, type=TriggerType.ON
    ),
    Trigger(program="gamess", trigger_strings=["GAMESS"], do_break=False, type=TriggerType.ON),
    Trigger(
        program="gamess",
        trigger_strings=["Firefly (PC GAMESS)"],
        do_break=True,
        type=TriggerType.ON,
    ),
    Trigger(
        program="gamess", trigger_strings=["GAMESS VERSION"], do_break=True, type=TriggerType.ON
    ),
    Trigger(
        program="gamessuk",
        trigger_strings=["G A M E S S - U K"],
        do_break=True,
        type=TriggerType.ON,
    ),
    Trigger(program="gamessdat", trigger_strings=["$DATA"], do_break=True, type=TriggerType.ON),
    Trigger(
        program="gaussian", trigger_strings=["Gaussian, Inc."], do_break=True, type=TriggerType.ON
    ),
    Trigger(program="jaguar", trigger_strings=["Jaguar"], do_break=True, type=TriggerType.ON),
    Trigger(program="molcas", trigger_strings=["MOLCAS"], do_break=True, type=TriggerType.ON),
    Trigger(
        program="molpro",
        trigger_strings=["PROGRAM SYSTEM MOLPRO"],
        do_break=True,
        type=TriggerType.ON,
    ),
    Trigger(program="molpro", trigger_strings=["1PROGRAM"], do_break=False, type=TriggerType.ON),
    Trigger(program="mopac", trigger_strings=["MOPAC20"], do_break=True, type=TriggerType.ON),
    Trigger(
        program="nbo",
        trigger_strings=["N A T U R A L   A T O M I C   O R B I T A L   A N D"],
        do_break=True,
        type=TriggerType.ON,
    ),
    Trigger(
        program="nwchem",
        trigger_strings=["Northwest Computational Chemistry Package"],
        do_break=True,
        type=TriggerType.ON,
    ),
    Trigger(program="orca", trigger_strings=["O   R   C   A"], do_break=True, type=TriggerType.ON),
    Trigger(
        program="psi4",
        trigger_strings=["Psi4: An Open-Source Ab Initio Electronic Structure Package"],
        do_break=True,
        type=TriggerType.ON,
    ),
    Trigger(
        program="qchem",
        trigger_strings=["A Quantum Leap Into The Future Of Chemistry"],
        do_break=True,
        type=TriggerType.ON,
    ),
    Trigger(program="turbomole", trigger_strings=["TURBOMOLE"], do_break=True, type=TriggerType.ON),
    Trigger(
        program="gaussian",
        trigger_strings=["Normal Termination of Gaussian"],
        do_break=True,
        type=TriggerType.OFF,
    ),
    Trigger(
        program="orca",
        trigger_strings=["****ORCA TERMINATED NORMALLY****"],
        do_break=True,
        type=TriggerType.OFF,
    ),
    Trigger(
        program="psi4",
        trigger_strings=["Psi4 exiting successfully. Buy a developer a beer!"],
        do_break=True,
        type=TriggerType.OFF,
    ),
]


# todo readerclasses = {
# todo     'cjson': cjsonreader.CJSON,
# todo     'json': cjsonreader.CJSON,
# todo     'xyz': xyzreader.XYZ,
# todo }
# todo
# todo writerclasses = {
# todo     'cjson': cjsonwriter.CJSON,
# todo     'json': cjsonwriter.CJSON,
# todo     'cml': cmlwriter.CML,
# todo     'molden': moldenwriter.MOLDEN,
# todo     'wfx': wfxwriter.WFXWriter,
# todo     'xyz': xyzwriter.XYZ,
# todo }


# todo class UnknownOutputFormatError(Exception):
# todo     """Raised when an unknown output format is encountered."""
# todo
# todo
# todo def is_xyz(inputfile: FileWrapper) -> bool:
# todo     """Is the given inputfile actually an XYZ file?
# todo
# todo     The only way to determine this without reading the entire file is to
# todo     inspect the file extension.
# todo     """
# todo     return len(inputfile.filenames) == 1 and \
# todo         os.path.splitext(inputfile.filenames[0])[1].lower() == ".xyz"


# todo def guess_filetype(inputfile) -> Optional[logfileparser.Logfile]:
# todo     """Try to guess the filetype by searching for trigger strings."""
# todo     filetype = None
# todo     logger = logging.getLogger("cclib")
# todo     try:
# todo         if isinstance(inputfile, FileWrapper) and is_xyz(inputfile):
# todo             logger.info("Found XYZ file based on file extension")
# todo             return filetype
# todo         for line in inputfile:
# todo             for parser, phrases, do_break in triggers:
# todo                 if all([line.lower().find(p.lower()) >= 0 for p in phrases]):
# todo                     filetype = parser
# todo                     if do_break:
# todo                         return filetype
# todo     except Exception:
# todo         # guess_filetype() is expected to be quiet by default...
# todo         logger.error("Failed to determine log file type", exc_info = True)
# todo
# todo     return filetype
# todo
# todo def sort_turbomole_outputs(fileinputs):
# todo     """
# todo     Sorts a list of inputs (or list of log files) according to the order
# todo     required by the Turbomole parser for correct parsing. Unrecognised
# todo     files are appended to the end of the list in the same order they are
# todo     given.
# todo
# todo     This function has been deprecated as of version 1.8; use:
# todo     cclib.parser.turbomoleparser.Turbomole.sort_input() instead
# todo
# todo     Inputs:
# todo       filelist - a list of Turbomole log files needed to be parsed.
# todo     Returns:
# todo       sorted_list - a sorted list of Turbomole files needed for proper parsing.
# todo     """
# todo     warnings.warn(
# todo         "sort_turbomole_outputs() has been deprecated as of v1.8; use: "+ \
# todo         "cclib.parser.turbomoleparser.Turbomole.sort_input() instead")
# todo     return Turbomole.sort_input(fileinputs)
# todo
# todo def ccread(
# todo         source: Union[str, typing.IO, FileWrapper, typing.List[Union[str, typing.IO]]],
# todo         *args,
# todo         **kwargs
# todo     ):
# todo     """Attempt to open and read computational chemistry data from a file.
# todo
# todo     If the file is not appropriate for cclib parsers, a fallback mechanism
# todo     will try to recognize some common chemistry formats and read those using
# todo     the appropriate bridge such as Open Babel.
# todo
# todo     Inputs:
# todo         source - a single logfile, a list of logfiles (for a single job),
# todo                  an input stream, or an URL pointing to a log file.
# todo         *args, **kwargs - arguments and keyword arguments passed to ccopen
# todo     Returns:
# todo         a ccData object containing cclib data attributes
# todo     """
# todo     log = None
# todo     try:
# todo         log = ccopen(source, *args, **kwargs)
# todo         logger = logging.getLogger("cclib")
# todo         if log:
# todo             logger.info("Identified logfile to be in {} format".format(type(log).__name__))
# todo
# todo             return log.parse()
# todo         else:
# todo             logger.info('Attempting to use fallback mechanism to read file')
# todo             return fallback(source)
# todo
# todo     finally:
# todo         if log:
# todo             log.inputfile.close()
# todo
# todo
# todo def ccopen(
# todo         source: Union[str, typing.IO, FileWrapper, typing.List[Union[str, typing.IO]]],
# todo         *args,
# todo         quiet: bool = False,
# todo         cjson: bool = False,
# todo         **kwargs
# todo     ):
# todo     """Guess the identity of a particular log file and return an instance of it.
# todo
# todo     Inputs:
# todo         source - a single logfile, a list of logfiles (for a single job),
# todo                  an input stream, or an URL pointing to a log file.
# todo         *args, **kwargs - arguments and keyword arguments passed to filetype
# todo
# todo     Returns:
# todo       one of ADF, DALTON, GAMESS, GAMESS UK, Gaussian, Jaguar,
# todo       Molpro, MOPAC, NWChem, ORCA, Psi3, Psi/Psi4, QChem, CJSON or None
# todo       (if it cannot figure it out or the file does not exist).
# todo     """
# todo     if not isinstance(source, list):
# todo         source = [source]
# todo
# todo     inputfile = None
# todo
# todo     logger = logging.getLogger("cclib")
# todo
# todo     try:
# todo         # Wrap our input with custom file object.
# todo         inputfile = FileWrapper(*source)
# todo
# todo         if cjson:
# todo             filetype = readerclasses['cjson']
# todo
# todo         else:
# todo             # Try and guess the parser we need.
# todo             filetype = guess_filetype(inputfile)
# todo
# todo             # Reset our position back to 0.
# todo             inputfile.reset()
# todo
# todo         # If the input file isn't a standard compchem log file, try one of
# todo         # the readers, falling back to Open Babel.
# todo         if not filetype:
# todo             # TODO: This assumes we only got a single file...
# todo             filename = list(inputfile.filenames)[0]
# todo             ext = pathlib.Path(filename).name[1:].lower()
# todo
# todo             for extension in readerclasses:
# todo                 if ext == extension:
# todo                     filetype = readerclasses[extension]
# todo
# todo         # Proceed to return an instance of the logfile parser only if the filetype
# todo         # could be guessed.
# todo         if filetype:
# todo                 return filetype(inputfile, *args, **kwargs)
# todo
# todo         elif inputfile is not None:
# todo             inputfile.close()
# todo             # Stop us closing twice in the except block.
# todo             inputfile = None
# todo
# todo         logger.warning(
# todo             "Unable to determine the type of logfile %s, try the fallback mechanism",
# todo             source
# todo         )
# todo
# todo     except Exception:
# todo         if inputfile is not None:
# todo             inputfile.close()
# todo
# todo         if not quiet:
# todo             raise
# todo
# todo         # We're going to swallow this exception if quiet is True.
# todo         # This can hide a lot of errors, so we'll make sure to log it.
# todo         logger.error("Failed to open logfile", exc_info = True)


# todo def ccwrite(ccobj, outputtype=None, outputdest=None,
# todo             indices=None, terse=False, returnstr=False,
# todo             *args, **kwargs):
# todo     """Write the parsed data from an outputfile to a standard chemical
# todo     representation.
# todo
# todo     Inputs:
# todo         ccobj - Either a job (from ccopen) or a data (from job.parse()) object
# todo         outputtype - The output format (should be a string)
# todo         outputdest - A filename or file object for writing
# todo         indices - One or more indices for extracting specific geometries/etc. (zero-based)
# todo         terse -  This option is currently limited to the cjson/json format. Whether to indent the cjson/json or not
# todo         returnstr - Whether or not to return a string representation.
# todo
# todo     The different writers may take additional arguments, which are
# todo     documented in their respective docstrings.
# todo
# todo     Returns:
# todo         the string representation of the chemical datatype
# todo           requested, or None.
# todo     """
# todo
# todo     # Determine the correct output format.
# todo     outputclass = _determine_output_format(outputtype, outputdest)
# todo
# todo     # Is ccobj an job object (unparsed), or is it a ccdata object (parsed)?
# todo     if isinstance(ccobj, logfileparser.Logfile):
# todo         jobfilename = ccobj.filename
# todo         ccdata = ccobj.parse()
# todo     elif isinstance(ccobj, data.ccData):
# todo         jobfilename = None
# todo         ccdata = ccobj
# todo     else:
# todo         raise ValueError
# todo
# todo     # If the logfile name has been passed in through kwargs (such as
# todo     # in the ccwrite script), make sure it has precedence.
# todo     if 'jobfilename' in kwargs:
# todo         jobfilename = kwargs['jobfilename']
# todo         # Avoid passing multiple times into the main call.
# todo         del kwargs['jobfilename']
# todo
# todo     outputobj = outputclass(ccdata, jobfilename=jobfilename,
# todo                             indices=indices, terse=terse,
# todo                             *args, **kwargs)
# todo     output = outputobj.generate_repr()
# todo
# todo     # If outputdest isn't None, write the output to disk.
# todo     if outputdest is not None:
# todo         if isinstance(outputdest, str):
# todo             with open(outputdest, 'w') as outputobj:
# todo                 outputobj.write(output)
# todo         elif isinstance(outputdest, io.IOBase):
# todo             outputdest.write(output)
# todo         else:
# todo             raise ValueError
# todo     # If outputdest is None, return a string representation of the output.
# todo     else:
# todo         return output
# todo
# todo     if returnstr:
# todo         return output
# todo
# todo
# todo def _determine_output_format(outputtype, outputdest):
# todo     """
# todo     Determine the correct output format.
# todo
# todo     Inputs:
# todo       outputtype - a string corresponding to the file type
# todo       outputdest - a filename string or file handle
# todo     Returns:
# todo       outputclass - the class corresponding to the correct output format
# todo     Raises:
# todo       UnknownOutputFormatError for unsupported file writer extensions
# todo     """
# todo
# todo     # Priority for determining the correct output format:
# todo     #  1. outputtype
# todo     #  2. outputdest
# todo
# todo     outputclass = None
# todo     # First check outputtype.
# todo     if isinstance(outputtype, str):
# todo         extension = outputtype.lower()
# todo         if extension in writerclasses:
# todo             outputclass = writerclasses[extension]
# todo         else:
# todo             raise UnknownOutputFormatError(extension)
# todo     else:
# todo         # Then checkout outputdest.
# todo         if isinstance(outputdest, str):
# todo             extension = os.path.splitext(outputdest)[1].lower()
# todo         elif isinstance(outputdest, io.IOBase):
# todo             extension = os.path.splitext(outputdest.name)[1].lower()
# todo         else:
# todo             raise UnknownOutputFormatError
# todo         if extension in writerclasses:
# todo             outputclass = writerclasses[extension]
# todo         else:
# todo             raise UnknownOutputFormatError(extension)
# todo
# todo     return outputclass
# todo
# todo
# todo def ccframe(ccobjs, *args, **kwargs):
# todo     """Returns a pandas.DataFrame of data attributes parsed by cclib from one
# todo     or more logfiles.
# todo
# todo     Inputs:
# todo         ccobjs - an iterable of either cclib jobs (from ccopen) or data (from
# todo         job.parse()) objects
# todo
# todo     Returns:
# todo         a pandas.DataFrame
# todo     """
# todo     _check_pandas(_has_pandas)
# todo     logfiles = []
# todo     for ccobj in ccobjs:
# todo         # Is ccobj an job object (unparsed), or is it a ccdata object (parsed)?
# todo         if isinstance(ccobj, logfileparser.Logfile):
# todo             jobfilename = ccobj.filename
# todo             ccdata = ccobj.parse()
# todo         elif isinstance(ccobj, data.ccData):
# todo             jobfilename = None
# todo             ccdata = ccobj
# todo         else:
# todo             raise ValueError
# todo
# todo         attributes = ccdata.getattributes()
# todo         attributes.update({
# todo             'jobfilename': jobfilename
# todo         })
# todo
# todo         logfiles.append(pd.Series(attributes))
# todo     return pd.DataFrame(logfiles)


class ccDriver:
    """Driver for cclib. This takes in one or many file or file-like objects. Creates the datastructure ccCollection. Creates the stateful file handler. Handle the processing of a parsing combinator."""

    def __init__(
        self,
        source: typing.Union[
            str, typing.IO, FileHandler, typing.List[typing.Union[str, typing.IO]]
        ],
        tree: Optional[Tree] = None,
        combinator: Optional["Combinator"] = None,
        loglevel: int = logging.ERROR,
        logname: str = "Log",
        logstream=sys.stderr,
        **kwds,
    ):
        """Initialise the ccDriver object.

        This should be called by a subclass in its own __init__ method.

        Inputs:
            source - a logfile, list of logfiles, or stream with at least a read method
            combinator - a way to parse data and structure ccCollection
            loglevel - integer corresponding to a log level from the logging module
            logstream - where to output the logging information
        """
        if not isinstance(source, FileHandler):
            source = FileHandler(source)

        self._combinator = combinator
        if tree is None:
            # default to single job
            tree = Tree()
            tree.add_root()
        self._tree = tree

        if combinator is None:
            combinator = AutoCombinator(self._tree)
        self._combinator = combinator

        self._ccCollection = ccCollection(self._combinator, self._tree)
        self._fileHandler = source
        self._identified_program: List[str] = []

    @property
    def cccollection(self) -> ccCollection:
        return self._ccCollection

    @property
    def fileHandler(self):
        return self._fileHandler

    @property
    def combinator(self) -> "Combinator":
        return self._combinator

    @property
    def tree(self) -> Tree:
        return self._tree

    @property
    def identified_program(self) -> Optional[str]:
        if not self._identified_program:
            return None
        else:
            return self._identified_program[-1]

    @identified_program.setter
    def identified_program(self, in_prog: Optional[str]) -> None:
        if in_prog is None:
            if self._identified_program:
                self._identified_program.pop()
        else:
            self._identified_program.append(in_prog)

    def process_combinator(self) -> ccCollection:
        """Process the combinator and populate the ccData objects in the ccCollection"""
        # Current:
        #
        # The object is already initialized with the tree and corresponding number of ccData instances in the collection.
        # Loop over all the lines in the (combined) file.  For each line...
        #   Look for triggers that start programs,  For each that matches, set the program on the ccData instance for that node, increasing the node counter.
        #   Look for triggers that stop programs.  For each that matches, increase the node counter.
        #   Run each parser in the combinator for the node id on the discovered program.
        #     parser -> parse() -> program_parse_method() -> returns data to go on collection's current ccData
        #
        # Future:
        #
        # Create an empty tree.
        # [maybe] Create a (empty) stack for indices into the tree. (TODO what indices?)
        # Loop over all the lines in the (combined) file.  For each line...
        #   Look for a program trigger.
        #     If the trigger matches and is to start a program,
        #       Create new ccData instance and append it to collection
        #       <update ccData reference to new>
        #     Elif the trigger matches and is to stop a program,
        #       <update ccData reference to previous>
        #     Else,
        #       continue
        #   Run each parser in the combinator for the current ccData
        #
        # At this point we have a linear sequence of ccData instances that need to be arranged into a hierarchical tree.
        # The structure of the tree is determined by the job run and the
        # program (Psi4 BSSE, Psi4 SAPT, Q-Chem EDA, Psi4 MBE, Q-Chem MBE,
        # Gaussian ONIOM, ...).
        self.identified_program = None
        line = next(self._fileHandler)
        current_idx = self._tree.get_next_idx()
        # TODO for line in self._fileHandler:
        while line:
            for trigger in TRIGGERS:
                if trigger.matches(line):
                    if trigger.type == TriggerType.ON:
                        if self.identified_program is None:
                            current_idx = 0
                        else:
                            # if a program is within a program this might mean things are ok but we proceed to a child node.. think about how to handle this?
                            current_idx = self._tree.get_next_idx()
                            assert current_idx is not None
                        self.identified_program = trigger.program
                        self._ccCollection.parsed_data[current_idx].setattributes(
                            {"metadata": {"identified_program": trigger.program}}
                        )
                        if trigger.do_break:
                            break
                    elif trigger.type == TriggerType.OFF:
                        self.identified_program = None
                        current_idx = self._tree.get_next_idx()
                        if trigger.do_break:
                            break
            if self.identified_program is None:
                line = next(self._fileHandler)
                continue
            assert current_idx is not None
            # right now combinator is just a list of list of subparsers, tree idxs access what list of parsers we are using
            for subparser in self._combinator.job_list[current_idx]:
                parsed_data = subparser.parse(
                    self._fileHandler,
                    self.identified_program,
                    self._ccCollection.parsed_data[current_idx],
                )
                if parsed_data is not None:
                    self._ccCollection.parsed_data[current_idx].setattributes(parsed_data)
            line = next(self._fileHandler)
        return self._ccCollection
