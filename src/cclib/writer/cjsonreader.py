"""
******************************************************************************

    License information

******************************************************************************/
"""
import json
from ..parser.cjsonHolder import cjsonHolder
from ..parser.data import ccData


class CJSON:
    """ CJSON log file"""

    def __init__(self, source, *args, **kwargs):

        # Set the filename to source if it is a string or a list of strings, which are
        # assumed to be filenames. Otherwise, assume the source is a file-like object
        # if it has a read method, and we will try to use it like a stream.
        if isinstance(source, str):
            self.filename = source
        else:
            raise ValueError

        self.datatype = cjsonHolder()

    def open_cjson(self):
        inputfile = self.filename

        # Actual update of cjsonHolder happens here
        json_data = open(inputfile).read()
        self.construct(json_data)

        # Removal of the keys with the Placeholder values
        #self.datatype = cjsonHolder().purge_cjson(self.datatype)

        self.datatype.compress()
        # Debugging statement. To-Do: Remove
        #print(self.datatype.cjson)

        return self.datatype



    def construct(self,input_file):
        # if input_file is None:
        #     raise TypeError
        data = json.loads(input_file)

        first_level = ['chemical json', 'name', 'smiles', 'inchi', 'inchikey', 'formula']
        for name in first_level:
            if name in data:
                self.datatype.cjson[name] = data[name]

        attributes = ['properties', 'atoms', 'optimization', 'vibrations', 'bonds', 'transitions', 'fragments']

        for name in attributes:
            if name in data:
                method_name = 'self.generate_' + name + '(data["' + name + '"])'
                eval(method_name)

    def generate_properties(self, input_dict):
        attr_list = ['charge', 'mult', 'enthalpy', 'entropy', 'natom', 'temperature', 'atomcharges']
        self.fill_cjson_dict(self.datatype.cjson['properties'], input_dict, attr_list)

        if "orbitals" in input_dict:
            orbital_attr = ['homos', 'moenergies', 'aooverlaps', 'mosyms', 'mocoeffs']
            self.fill_cjson_dict(self.datatype.cjson['properties']['orbitals'], input_dict["orbitals"], orbital_attr)

        if "energy" in input_dict:
            self.datatype.cjson['properties']['energy']['alpha']['homo'] = input_dict['energy']['alpha']['homo']
            self.datatype.cjson['properties']['energy']['alpha']['gap'] = input_dict['energy']['alpha']['gap']
            self.datatype.cjson['properties']['energy']['beta']['homo'] = input_dict['energy']['beta']['homo']
            self.datatype.cjson['properties']['energy']['beta']['gap'] = input_dict['energy']['beta']['gap']
            self.datatype.cjson['properties']['energy']['total'] = input_dict['energy']['total']
            energy_attr = ['freeenergy', 'mpenergies', 'ccenergies']
            self.fill_cjson_dict(self.datatype.cjson['properties']['energy'], input_dict['energy'], energy_attr)

    def generate_atoms(self, input_dict):
        if 'elements' in input_dict:
            attr_list = ['heavy atom count', 'atom count', ccData._attributes['atomnos'].jsonKey]
            for name in attr_list:
                if name in input_dict['elements']:
                    self.datatype.cjson['atoms']['elements'][name] = input_dict['elements'][name]

        if 'coords' in input_dict:
            self.datatype.cjson['atoms']['coords']['3d'] = input_dict['coords']['3d']

        if 'orbitals' in input_dict:
            orbital_list = ['aonames', 'atombasis']
            self.fill_cjson_dict(self.datatype.cjson['atoms']['orbitals'], input_dict['orbitals'], orbital_list)

        key_list = ['coreelectrons', 'atommasses', 'atomspins']
        self.fill_cjson_dict(self.datatype.cjson['atoms'], input_dict, key_list)

    def generate_optimization(self, input_dict):
        attr_list = ['optdone', 'optstatus', 'geotargets', 'geovalues', 'nbasis', 'nmo']
        self.fill_cjson_dict(self.datatype.cjson['optimization'], input_dict, attr_list)

        if 'scf' in input_dict:
            scf_list = ['scfenergies', 'scftargets', 'scfvalues']
            self.fill_cjson_dict(self.datatype.cjson['optimization']['scf'], input_dict['scf'], scf_list)

        if 'scan' in input_dict:
            attr_list = ['scancoords', 'scanenergies', 'scannames', 'scanparm']
            self.fill_cjson_dict(self.datatype.cjson['optimization']['scan'], input_dict['scan'], attr_list)

    def generate_vibrations(self, input_dict):
        attr_list = ['vibanharms', 'vibfreqs', 'vibsyms', 'hessian', 'vibdisps']
        self.fill_cjson_dict(self.datatype.cjson['vibrations'], input_dict, attr_list)

        if 'intensities' in input_dict:
            intensities_list = ['vibirs', 'vibramans']
            self.fill_cjson_dict(self.datatype.cjson['vibrations']['intensities'], input_dict['intensities'],
                                 intensities_list)

    def generate_bonds(self, input_dict):
        if 'connections' in input_dict:
            self.datatype.cjson['bonds']['connections']['index'] = input_dict['connections']['index']
        if 'order' in input_dict:
            self.datatype.cjson['bonds']['order'] = input_dict['order']

    def generate_transitions(self, input_dict):
        attr_list = ['etenergies', 'etoscs', 'etrotats', 'etsecs', 'etsyms']
        self.fill_cjson_dict(self.datatype.cjson['transitions'], input_dict, attr_list)

    def generate_fragments(self, input_dict):
        attr_list = ['fragnames', 'frags', 'fonames', 'fooverlaps']
        self.fill_cjson_dict(self.datatype.cjson['fragments'], input_dict, attr_list)


    def fill_cjson_dict(self, cjson_dict, input_dict, attr_list):
        key_list = [ccData._attributes[key].jsonKey for key in attr_list]
        for name in key_list:
            if name in input_dict:
                cjson_dict[name] = input_dict[name]
