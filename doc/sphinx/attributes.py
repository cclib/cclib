# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Generate the attributes.rst and attributes_dev.rst files from the
ccData docstring that describes attributes."""

import cclib

from docs_common import check_cclib

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
    attributes = [line for line in data_doc.split("\n") if line[:8].strip() == ""]
    attributes = [line for line in attributes if "--" in line]

    names = []
    descriptions = []
    units = []
    types = []
    for line in attributes:
        # There is always a double dash after the name.
        attr, desc = line.strip().split(" -- ")

        # The type and unit are in parentheses, but these
        # are not always the only parentheses on the line.
        other = desc.split("(")[-1]
        desc = desc[: -len(other) - 1].strip()
        other = other.split(")")[0]

        # Furthermore, the unit is not always there.
        if "," in other:
            atype, aunit = other.split(", ")
        else:
            atype = other
            aunit = ""

        for i in range(1, 4):
            atype = atype.replace(f"[{int(i)}]", f" of rank {int(i)}")

        names.append(attr)
        descriptions.append(desc)
        units.append(aunit)
        types.append(atype)

    # These are the widths of the columns in the table
    wattr = 4 + max(len(attr) for attr in names)
    wdesc = 1 + max(len(desc) for desc in descriptions)
    wunit = 1 + max(len(aunit) for aunit in units)
    wtype = 1 + max(len(atype) for atype in types)

    dashes = "    "
    for w in [wattr, wdesc, wunit, wtype]:
        dashes += f"{'=' * (w - 1)} "
    header = "    "
    header += "Name".ljust(wattr)
    header += "Description".ljust(wdesc)
    header += "Units".ljust(wunit)
    header += "Data type".ljust(wtype)
    lines.append(dashes)
    lines.append(header)
    lines.append(dashes)

    for attr, desc, aunit, atype in zip(names, descriptions, units, types):
        # Print the line with columns align to the table. Note that
        # the description sometimes contain Unicode characters, so
        # decode-encode when justifying to get the correct length.
        attr = f"`{attr}`_".ljust(wattr)
        desc = desc.ljust(wdesc)
        aunit = aunit.ljust(wunit)

        lines.append(f"    {attr}{desc}{aunit}{atype}")

    lines.append(dashes)
    lines.append("")

    for n in names:
        lines.append(f".. _`{n}`: data_notes.html#{n}")

    return "\n".join(lines)


if __name__ == "__main__":
    print(generate_attributes())
