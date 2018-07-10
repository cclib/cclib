# -*- coding: utf-8 -*-

"""Generate the attributes.rst and attributes_dev.rst files from the
ccData docstring that describes attributes."""

from __future__ import print_function

from docs_common import check_cclib

import cclib
check_cclib(cclib)


def generate_attributes():
    """Generate a string containing a reStructuredText table
    representation of the ccData docstring, which contains a list of
    all supported attributes with
    1. the name of each attribute,
    2. the text definition of each attribute,
    3. the attribute's container data type, shape (if relevant), and
    4. the physical units for each attribute.
    """
    lines = []

    # Need to parse the ccData docstring, since only that currently
    # contains all the information needed for this table.
    data_doc = cclib.parser.data.ccData.__doc__
    attributes = [line for line in data_doc.split('\n') if line[:8].strip() == '']
    attributes = [line for line in attributes if "--" in line]

    # These are the widths of the columns in the table
    wattr = 20
    wdesc = 65
    wunit = 28
    wtype = 32

    dashes = "    "
    for w in [wattr, wdesc, wunit, wtype]:
        dashes += "="*(w-1) + " "
    header = "    "
    header += "Name".ljust(wattr)
    header += "Description".ljust(wdesc)
    header += "Units".ljust(wunit)
    header += "Data type".ljust(wtype)
    lines.append(dashes)
    lines.append(header)
    lines.append(dashes)

    names = []
    for line in attributes:

        # There is always a double dash after the name.
        attr, desc = line.strip().split(' -- ')
        names.append(attr)

        # The type and unit are in parentheses, but these
        # are not always the only parentheses on the line.
        other = desc.split('(')[-1]
        desc = desc[:-len(other)-1].strip()
        other = other.split(')')[0]

        # Furthermore, the unit is not always there.
        if "," in other:
            atype, aunit = other.split(", ")
        else:
            atype = other
            aunit = ''

        # Print the line with columns align to the table. Note that
        # the description sometimes contain Unicode characters, so
        # decode-encode when justifying to get the correct length.
        attr = ("`%s`_" % attr).ljust(wattr)
        desc = desc.ljust(wdesc)
        aunit = aunit.ljust(wunit)
        for i in range(1, 4):
            atype = atype.replace('[%i]' % i, ' of rank %i' % i)
        lines.append("    " + attr + desc + aunit + atype)

    lines.append(dashes)
    lines.append("")

    for n in names:
        lines.append(".. _`%s`: data_notes.html#%s" % (n, n))

    return "\n".join(lines)

if __name__ == "__main__":
    print(generate_attributes())
