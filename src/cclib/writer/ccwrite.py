# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""..."""

import os.path

from . import cjsonwriter
from . import cmlwriter
from . import xyzwriter


def ccwrite(ccjob, outputtype, outputdest=None, *args, **kwargs):
    """Write the parsed data from an outputfile to a standard chemical
    representation.

    Inputs:
        ccjob - ...
        outputtype - ...
        outputdest - ...

    Returns:
        the string representation of the chemical datatype
          requested, or None.
    """

    # Determine the correct output format based on outputtype.
    if outputtype.lower() == 'cjson':
        outputclass = cjsonwriter.CJSON
    elif outputtype.lower() == 'cml':
        outputclass = cmlwriter.CML
    elif outputtype.lower() == 'xyz':
        outputclass = xyzwriter.XYZ
    else:
        raise ValueError

    jobfilename = ccjob.filename
    ccdata = ccjob.parse()
    outputobj = outputclass(ccdata, jobfilename=jobfilename, *args, **kwargs)
    output = outputobj.generate_repr()

    # If outputdest isn't None, write the output to disk.
    if outputdest is not None:
        if isinstance(outputdest, str):
            with open(outputdest, 'w') as outputobj:
                outputobj.write(output)
        elif isinstance(outputdest, file):
            outputdest.write(output)
        else:
            raise ValueError

    # Always return a string representation of the output.
    return output
