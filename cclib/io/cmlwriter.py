# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A writer for chemical markup language (CML) files."""

import xml.etree.ElementTree as ET

from cclib.io import filewriter
from cclib.parser.utils import find_package

_has_openbabel = find_package("openbabel")


class CML(filewriter.Writer):
    """A writer for chemical markup language (CML) files."""

    def __init__(self, ccdata, *args, **kwargs):
        """Initialize the CML writer object.

        Inputs:
          ccdata - An instance of ccData, parsed from a logfile.
        """
        super().__init__(ccdata, *args, **kwargs)

    def generate_repr(self):
        """Generate the CML representation of the logfile data."""

        # Create the base molecule.
        molecule = ET.Element('molecule')
        d = {
            # Write the namespace directly.
            'xmlns': 'http://www.xml-cml.org/schema',
        }
        if self.jobfilename is not None:
            d['id'] = self.jobfilename
        _set_attrs(molecule, d)

        # Form the listing of all the atoms present.
        atomArray = ET.SubElement(molecule, 'atomArray')
        if hasattr(self.ccdata, 'atomcoords') and hasattr(self.ccdata, 'atomnos'):
            elements = [self.pt.element[Z] for Z in self.ccdata.atomnos]
            for atomid in range(self.ccdata.natom):
                atom = ET.SubElement(atomArray, 'atom')
                x, y, z = self.ccdata.atomcoords[-1][atomid].tolist()
                d = {
                    "id": f"a{atomid + 1}",
                    "elementType": elements[atomid],
                    "x3": f"{x:.10f}",
                    "y3": f"{y:.10f}",
                    "z3": f"{z:.10f}",
                }
                _set_attrs(atom, d)

        # Form the listing of all the bonds present.
        bondArray = ET.SubElement(molecule, 'bondArray')
        if _has_openbabel:
            for bc in self.bond_connectivities:
                bond = ET.SubElement(bondArray, "bond")
                d = {"atomRefs2": f"a{bc[0] + 1} a{bc[1] + 1}", "order": str(bc[2])}
                _set_attrs(bond, d)

        _indent(molecule)

        return _tostring(molecule)


def _set_attrs(element, d):
    """Set all the key-value pairs from a dictionary as element
    attributes.
    """
    for (k, v) in d.items():
        element.set(k, v)
    return


def _indent(elem, level=0):
    """An in-place pretty-print indenter for XML."""
    i = f"\n{level * '  '}"
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = f"{i}  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            _indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def _tostring(element, xml_declaration=True, encoding='utf-8', method='xml'):
    """A reimplementation of tostring() found in ElementTree."""
    class dummy:
        pass
    data = []
    file = dummy()
    file.write = data.append
    ET.ElementTree(element).write(file,
                                  xml_declaration=xml_declaration,
                                  encoding=encoding,
                                  method=method)
    return b''.join(data).decode(encoding)


del find_package
