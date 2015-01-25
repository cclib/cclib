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

import openbabel as ob

import xml.etree.cElementTree as ET

from collections import OrderedDict

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

        # Generate the Open Babel/Pybel representation of the molecule.
        # Used for calculating SMILES/InChI, formula, MW, etc.
        obmol, pbmol = self._make_openbabel_from_ccdata()
        bond_connectivities = self._make_bond_connectivity_from_openbabel(obmol)

        # Create the base molecule.
        molecule = ET.Element('molecule')
        d = {
            'xmlns': 'http://www.xml-cml.org/schema',
        }
        if self.jobfilename is not None:
            d['id'] = self.jobfilename
        _set_attrs(molecule, d)

        # Form the listing of all the atoms present.
        atomArray = ET.SubElement(molecule, 'atomArray')
        for obatom in ob.OBMolAtomIter(obmol):
            atom = ET.SubElement(atomArray, 'atom')
            atomid = obatom.GetIndex()
            x, y, z = self.ccdata.atomcoords[-1][atomid].tolist()
            d = [
                ('id', 'a{}'.format(atomid + 1)),
                ('elementType', self.elements[atomid]),
                ('x3', '{:.10f}'.format(x)),
                ('y3', '{:.10f}'.format(y)),
                ('z3', '{:.10f}'.format(z)),
            ]
            _set_attrs(atom, OrderedDict(d))

        # Form the listing of all the bonds present.
        bondArray = ET.SubElement(molecule, 'bondArray')
        for bc in bond_connectivities:
            bond = ET.SubElement(bondArray, 'bond')
            d = [
                ('atomRefs2', 'a{} a{}'.format(bc[0] + 1, bc[1] + 1)),
                ('order', str(bc[2]))
            ]
            _set_attrs(bond, OrderedDict(d))

        _indent(molecule)

        return ET.tostring(molecule)


def _set_attrs(element, d):
    """Set all the key-value pairs from a dictionary as element
    attributes.
    """
    for (k, v) in d.iteritems():
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


if __name__ == "__main__":
    pass
