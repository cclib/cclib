# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""The top-level interface for writing parsed quantum chemical data to
standardized chemical representations, such as chemical JSON, CML, and
Cartesian coordinates."""

import os.path
import sys
# Python 2->3 changes the default file object hierarchy.
if sys.version_info[0] == 2:
    fileclass = file
else:
    import io
    fileclass = io.IOBase

from .. import parser

from . import cjsonwriter
from . import cmlwriter
from . import xyzwriter


def ccwrite(ccobj, outputtype=None, outputdest=None, returnstr=False,
            *args, **kwargs):
    """Write the parsed data from an outputfile to a standard chemical
    representation.

    Inputs:
        ccobj - Either a job (from ccopen) or a data (from job.parse()) object
        outputtype - The output format (should be one of 'cjson', 'cml', 'xyz')
        outputdest - A filename or file object for writing
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
    if isinstance(ccobj, parser.logfileparser.Logfile):
        jobfilename = ccobj.filename
        ccdata = ccobj.parse()
    elif isinstance(ccobj, parser.data.ccData):
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

    outputobj = outputclass(ccdata, jobfilename=jobfilename, *args, **kwargs)
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


if __name__ == "__main__":
    pass
