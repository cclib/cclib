# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""A writer for chemical markup language (CML) files."""

try:
    import openbabel as ob
    has_openbabel = True
except ImportError:
    has_openbabel = False

import xml.etree.cElementTree as ET

from . import filewriter


class CML(filewriter.Writer):
    """A writer for chemical markup language (CML) files."""

    def __init__(self, ccdata, *args, **kwargs):
        """Initialize the CML writer object.

        Inputs:
          ccdata - An instance of ccData, parsed from a logfile.
        """

        # Call the __init__ method of the superclass
        super(CML, self).__init__(ccdata, *args, **kwargs)

        self.generate_repr()

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
        for atomid in range(self.ccdata.natom):
            atom = ET.SubElement(atomArray, 'atom')
            x, y, z = self.ccdata.atomcoords[-1][atomid].tolist()
            d = {
                'id': 'a{}'.format(atomid + 1),
                'elementType': self.elements[atomid],
                'x3': '{:.10f}'.format(x),
                'y3': '{:.10f}'.format(y),
                'z3': '{:.10f}'.format(z),
            }
            _set_attrs(atom, d)

        # Form the listing of all the bonds present.
        bondArray = ET.SubElement(molecule, 'bondArray')
        if has_openbabel:
            for bc in self.bond_connectivities:
                bond = ET.SubElement(bondArray, 'bond')
                d = {
                    'atomRefs2': 'a{} a{}'.format(bc[0] + 1, bc[1] + 1),
                    'order': str(bc[2]),
                }
                _set_attrs(bond, d)

        _indent(molecule)

        return _tostring(molecule)


def _set_attrs(element, d):
    """Set all the key-value pairs from a dictionary as element
    attributes.
    """
    for (k, v) in d.items():
        element.set(k, v)


def _indent(elem, level=0):
    """An in-place pretty-print indenter for XML."""
    i = "\n" + (level * "  ")
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
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


if __name__ == "__main__":
    pass
