# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Run data tests for cclib."""

import logging
import os
import sys
from typing import Optional, Type, Union

import cclib

__filedir__ = os.path.realpath(os.path.dirname(__file__))


# We need this in Python3 for importing things from the same directory
# within the unit test files.
sys.path.insert(1, os.path.join(__filedir__, "data"))

[  # "ADF",
    # "DALTON",
    # "FChk",
    # "GAMESS",
    # "GAMESSDAT",
    # "GAMESSUK",
    "Gaussian",
    # "Jaguar",
    # "Molpro",
    # "Molcas",
    # "MOPAC",
    # "Molcas",
    # "NWChem",
    "NBO",
    # "NWChem",
    # "QChem",
    # "Turbomole",
    # "XTB",
    # "Turbomole",
]
# all_parsers = {name: getattr(cclib.parser, name) for name in parser_names}


def get_program_dir(parser_name: str) -> str:
    """Return a directory name given a parser name.

    In at least one case (GAMESS-UK) the directory is named differently.
    """

    if parser_name == "GAMESSUK":
        return "GAMESS-UK"
    return parser_name


def getdatafile(
    parser: Union[str, Type[cclib.file_handler.FileHandler]],
    subdir,
    files,
    loglevel: int = logging.ERROR,
    datatype: Optional[Type[cclib.parser_properties.data.ccData]] = None,
):
    """Returns a parsed logfile.

    Inputs:
        parser   - a logfile parser class (subclass of LogFile)
        subdir   - subdirectory containing data files (program version)
        files    - data filename(s)
        stream   - where to log to (sys.stdout by default)
        loglevel - what level to log at
        datatype - ccData or child class

    Outputs:
        data - the resulting data object
        logfile - the parser object used for parsing
    """

    # Convert any string into the parser object we will be using.
    if isinstance(parser, str):
        parser = all_parsers[parser]

    datadir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
    programdir = os.path.join(get_program_dir(parser.__name__), subdir)
    inputs = [os.path.join(datadir, programdir, fn) for fn in files]

    # We should be able to pass a list of length one here, but for some reason
    # this does not work with some parsers and we get errors.
    if len(inputs) == 1:
        inputs = inputs[0]

    stream = stream or sys.stdout
    logfile = parser(
        inputs, logstream=stream, loglevel=loglevel, datatype=datatype or cclib.parser.data.ccData
    )

    data = logfile.parse()
    return data, logfile
