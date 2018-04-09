# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A writer for chemical JSON (CJSON) files."""

try:
    import openbabel as ob
    _has_openbabel = True
except ImportError:
    _has_openbabel = False

import os.path
import json
import numpy as np

from cclib.io import filewriter
from cclib.parser.data import ccData


class CJSON(filewriter.Writer):
    """A writer for chemical JSON (CJSON) files."""
    def __init__(self, ccdata, terse=False, *args, **kwargs):
        """Initialize the chemical JSON writer object.

        Inputs:
          ccdata - An instance of ccData, parsed from a logfile.
        """

        # Call the __init__ method of the superclass
        super(CJSON, self).__init__(ccdata, terse=terse, *args, **kwargs)

    def pathname(self, path):
        """Return filename without extension to be used as name."""
        name = os.path.basename(os.path.splitext(path)[0])
        return name

    def generate_repr(self):
        """Generate the CJSON representation of the logfile data."""

        cjson_dict = dict()
        # Need to decide on a number format.
        cjson_dict['chemical json'] = 0
        if self.jobfilename is not None:
            cjson_dict['name'] = self.pathname(self.jobfilename)

        # These are properties that can be collected using Open Babel.
        if _has_openbabel:
            cjson_dict['smiles'] = self.pbmol.write('smiles')
            cjson_dict['inchi'] = self.pbmol.write('inchi')
            cjson_dict['inchikey'] = self.pbmol.write('inchikey')
            cjson_dict['formula'] = self.pbmol.formula
        # Incorporate Unit Cell into the chemical JSON.

        # Iterate through the attribute list present in ccData. Depending on the
        # availability of the attribute add it at the right 'level'.
        for attributeName, Value in ccData._attributes.items():
            if not hasattr(self.ccdata, attributeName):
                continue

            attributePath = Value.attributePath.split(":")

            # Depth of the attribute in the CJSON.
            levels = len(attributePath)

            # The attributes which haven't been included in the CJSON format.
            if attributePath[0] == 'N/A':
                continue

            if attributePath[0] not in cjson_dict:
                cjson_dict[attributePath[0]] = dict()
            l1_data_object = cjson_dict[attributePath[0]]

            # 'moments' and 'atomcoords' key will contain processed data obtained from the output file.
            if attributeName == 'moments' or attributeName == 'atomcoords' :
                if attributeName == 'moments':
                    cjson_dict['properties'][ccData._attributes['moments'].jsonKey] = self._calculate_total_dipole_moment()
                else:
                    cjson_dict['atoms']['coords'] = dict()
                    cjson_dict['atoms']['coords']['3d'] = self.ccdata.atomcoords[-1].flatten().tolist()
                continue

            if levels == 1:
                self.set_JSON_attribute(l1_data_object, attributeName)
            elif levels >= 2:
                if attributePath[1] not in l1_data_object:
                    l1_data_object[attributePath[1]] = dict()
                l2_data_object = l1_data_object[attributePath[1]]

                if levels == 2:
                    self.set_JSON_attribute(l2_data_object, attributeName)
                elif levels == 3:
                    if attributePath[2] not in l2_data_object:
                        l2_data_object[attributePath[2]] = dict()
                    l3_data_object = l2_data_object[attributePath[2]]
                    self.set_JSON_attribute(l3_data_object, attributeName)

        # Attributes which are not directly obtained from the output files.
        if hasattr(self.ccdata, 'moenergies') and hasattr(self.ccdata, 'homos'):
            if 'energy' not in cjson_dict['properties']:
                cjson_dict['properties']['energy'] = dict()

            cjson_dict['properties']['energy']['alpha'] = dict()
            cjson_dict['properties']['energy']['beta'] = dict()

            homo_idx_alpha = int(self.ccdata.homos[0])
            homo_idx_beta = int(self.ccdata.homos[-1])
            energy_alpha_homo = self.ccdata.moenergies[0][homo_idx_alpha]
            energy_alpha_lumo = self.ccdata.moenergies[0][homo_idx_alpha + 1]
            energy_alpha_gap = energy_alpha_lumo - energy_alpha_homo
            energy_beta_homo = self.ccdata.moenergies[-1][homo_idx_beta]
            energy_beta_lumo = self.ccdata.moenergies[-1][homo_idx_beta + 1]
            energy_beta_gap = energy_beta_lumo - energy_beta_homo

            cjson_dict['properties']['energy']['alpha']['homo'] = energy_alpha_homo
            cjson_dict['properties']['energy']['alpha']['gap'] = energy_alpha_gap
            cjson_dict['properties']['energy']['beta']['homo'] = energy_beta_homo
            cjson_dict['properties']['energy']['beta']['gap'] = energy_beta_gap
            cjson_dict['properties']['energy']['total'] = self.ccdata.scfenergies[-1]

        if hasattr(self.ccdata, 'atomnos'):
            cjson_dict['atoms']['elements']['atom count'] = len(self.ccdata.atomnos)
            cjson_dict['atoms']['elements']['heavy atom count'] = len([x for x in self.ccdata.atomnos if x > 1])

        # Bond attributes:
        if _has_openbabel and (len(self.ccdata.atomnos) > 1):
            cjson_dict['bonds'] = dict()
            cjson_dict['bonds']['connections'] = dict()
            cjson_dict['bonds']['connections']['index'] = []
            for bond in self.bond_connectivities:
                cjson_dict['bonds']['connections']['index'].append(bond[0] + 1)
                cjson_dict['bonds']['connections']['index'].append(bond[1] + 1)
            cjson_dict['bonds']['order'] = [bond[2] for bond in self.bond_connectivities]

        if _has_openbabel:
            cjson_dict['properties']['molecular mass'] = self.pbmol.molwt
            cjson_dict['diagram'] = self.pbmol.write(format='svg')

        if self.terse:
            return json.dumps(cjson_dict, cls=NumpyAwareJSONEncoder)
        else:
            return json.dumps(cjson_dict, cls=JSONIndentEncoder, sort_keys=True, indent=4)

    def set_JSON_attribute(self, object, key):
        """
        Args:
            object: Python dictionary which is being appended with the key value.
            key: cclib attribute name.

        Returns: 
            None. The dictionary is modified to contain the attribute with the
                 cclib keyname as key
        """
        if hasattr(self.ccdata, key):
            object[ccData._attributes[key].jsonKey] = getattr(self.ccdata, key)

class NumpyAwareJSONEncoder(json.JSONEncoder):
    """A encoder for numpy.ndarray's obtained from the cclib attributes.
       For all other types the json default encoder is called.
       Do Not rename the 'default' method as it is required to be implemented
       by any subclass of the json.JSONEncoder
    """
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            if obj.ndim == 1:
                nan_list = obj.tolist()
                return [None if np.isnan(x) else x for x in nan_list]
            else:
                return [self.default(obj[i]) for i in range(obj.shape[0])]
        return json.JSONEncoder.default(self, obj)


class JSONIndentEncoder(json.JSONEncoder):
    def __init__(self, *args, **kwargs):
        super(JSONIndentEncoder, self).__init__(*args, **kwargs)
        self.current_indent = 0
        self.current_indent_str = ""

    def encode(self, o):
        # Special Processing for lists
        if isinstance(o, (list, tuple)):
            primitives_only = True
            for item in o:
                if isinstance(item, (list, tuple, dict)):
                    primitives_only = False
                    break
            output = []
            if primitives_only:
                for item in o:
                    output.append(json.dumps(item,  cls=NumpyAwareJSONEncoder))
                return "[ " + ", ".join(output) + " ]"
            else:
                self.current_indent += self.indent
                self.current_indent_str = "".join([" " for x in range(self.current_indent)])
                for item in o:
                    output.append(self.current_indent_str + self.encode(item))
                self.current_indent -= self.indent
                self.current_indent_str = "".join([" " for x in range(self.current_indent)])
                return "[\n" + ",\n".join(output) + "\n" + self.current_indent_str + "]"
        elif isinstance(o, dict):
            output = []
            self.current_indent += self.indent
            self.current_indent_str = "".join([" " for x in range(self.current_indent)])
            for key, value in o.items():
                output.append(self.current_indent_str + json.dumps(key, cls=NumpyAwareJSONEncoder) + ": " +
                              str(self.encode(value)))
            self.current_indent -= self.indent
            self.current_indent_str = "".join([" " for x in range(self.current_indent)])
            return "{\n" + ",\n".join(output) + "\n" + self.current_indent_str + "}"
        elif isinstance(o, np.generic):
            return json.dumps(np.asscalar(o), cls=NumpyAwareJSONEncoder)
        else:
            return json.dumps(o, cls=NumpyAwareJSONEncoder)
