# -*- coding: utf-8 -*-
#
# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2007-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

from .data import ccData

class cjsonHolder:
    # replicate the structure of a CJSON format found here:
    # https://docs.google.com/document/d/1_RYFXzhxHK525id0A930Pa1y38Ui2X5GgtAo68iE5Oc/edit?usp=sharing
    CONST_VALUE = None

    def construct_cjson(self):
        # Initial values
        self.cjson['chemical json'] = 0
        self.cjson['name'] = self.CONST_VALUE
        self.cjson['smiles'] = self.CONST_VALUE
        self.cjson['inchi'] = self.CONST_VALUE
        self.cjson['inchikey'] = self.CONST_VALUE
        self.cjson['formula'] = self.CONST_VALUE

        self.construct_properties()
        self.construct_atoms()
        self.construct_optimization()
        self.construct_vibrations()
        self.construct_bonds()
        self.construct_transitions()
        self.construct_fragments()

    def set_default_value(self, input_dict, input_list):
        for key in input_list:
            input_dict[ccData._attributes[key].jsonKey] = self.CONST_VALUE

    def construct_properties(self):
        # Generate Properties
        self.cjson['properties'] = dict()

        attr_list = ['charge', 'mult', 'enthalpy', 'entropy', 'natom', 'temperature', 'atomcharges']
        self.set_default_value(self.cjson["properties"], attr_list)

        self.cjson["properties"]["energy"] = dict()
        self.cjson["properties"]["energy"]["alpha"] = dict()
        self.cjson["properties"]["energy"]["beta"] = dict()
        self.cjson['properties']['energy']['alpha']['homo'] = self.CONST_VALUE
        self.cjson['properties']['energy']['alpha']['gap'] = self.CONST_VALUE
        self.cjson['properties']['energy']['beta']['homo'] = self.CONST_VALUE
        self.cjson['properties']['energy']['beta']['gap'] = self.CONST_VALUE
        self.cjson['properties']['energy']['total'] = self.CONST_VALUE
        energy_attr = ['freeenergy', 'mpenergies', 'ccenergies']
        self.set_default_value(self.cjson["properties"]["energy"], energy_attr)

        self.cjson["properties"]["orbitals"] = dict()
        orbital_attr = ['homos', 'moenergies', 'aooverlaps', 'mosyms', 'mocoeffs']
        self.set_default_value(self.cjson["properties"]["orbitals"], orbital_attr)

    def construct_atoms(self):
        self.cjson['atoms'] = dict()

        self.cjson['atoms']['elements'] = dict()
        self.cjson['atoms']['elements'][ccData._attributes['atomnos'].jsonKey] = self.CONST_VALUE
        self.cjson['atoms']['elements']['atom count'] = self.CONST_VALUE
        self.cjson['atoms']['elements']['heavy atom count'] = self.CONST_VALUE

        self.cjson['atoms']['coords'] = dict()
        self.cjson['atoms']['coords']['3d'] = self.CONST_VALUE

        orbital_list = ['aonames', 'atombasis']
        self.cjson['atoms']['orbitals'] = dict()
        self.set_default_value(self.cjson['atoms']['orbitals'], orbital_list)

        self.set_default_value(self.cjson['atoms'], ['coreelectrons', 'atommasses', 'atomspins'])

    def construct_optimization(self):
        self.cjson['optimization'] = dict()
        attr_list = ['optdone', 'optstatus', 'geotargets', 'geovalues', 'nbasis', 'nmo']
        self.set_default_value(self.cjson['optimization'], attr_list)

        self.cjson['optimization']['scf'] = dict()
        scf_list = ['scfenergies', 'scftargets', 'scfvalues']
        self.set_default_value(self.cjson['optimization']['scf'], scf_list)

        self.cjson['optimization']['scan'] = dict()
        attr_list = ['scancoords', 'scanenergies', 'scannames', 'scanparm']
        self.set_default_value(self.cjson['optimization']['scan'], attr_list)

    def construct_vibrations(self):
        self.cjson['vibrations'] = dict()
        attr_list = ['vibanharms', 'vibfreqs', 'vibsyms', 'hessian', 'vibdisps']
        self.set_default_value(self.cjson['vibrations'], attr_list)

        self.cjson['vibrations']['intensities'] = dict()
        intensities_list = ['vibirs', 'vibramans']
        self.set_default_value(self.cjson['vibrations']['intensities'], intensities_list)

    def construct_bonds(self):
        self.cjson['bonds'] = dict()
        self.cjson['bonds']['connections'] = dict()
        self.cjson['bonds']['connections']['index'] = self.CONST_VALUE
        self.cjson['bonds']['order'] = self.CONST_VALUE

    def construct_transitions(self):
        attr_list = ['etenergies', 'etoscs', 'etrotats', 'etsecs', 'etsyms']
        self.cjson['transitions'] = dict()
        self.set_default_value(self.cjson['transitions'], attr_list)

    def construct_fragments(self):
        attr_list = ['fragnames', 'frags', 'fonames', 'fooverlaps']
        self.cjson['fragments'] = dict()
        self.set_default_value(self.cjson['fragments'], attr_list)

    def purge_cjson(self, input_dict):
        if not isinstance(input_dict, (dict, list)):
            return input_dict
        if isinstance(input_dict, list):
            return [v for v in (self.purge_cjson(v) for v in input_dict) if v]
        return {k: v for k, v in ((k, self.purge_cjson(v)) for k, v in input_dict.items()) if v}

    def compress(self):
        self.cjson = self.purge_cjson(self.cjson)
        #print(self.cjson)

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, item):
        return setattr(self, key, item)

    def __finditem__(self, obj, key):
        if key in obj: return obj[key]
        for k, v in obj.items():
            if isinstance(v, dict):
                item = self.__finditem__(v, key)
                if item is not None:
                    return item

    def __getattr__(self, name):
        return self.__finditem__(self.cjson, ccData._attributes[name].jsonKey)

    def __init__(self):
        self.cjson = dict()
        self.construct_cjson()

