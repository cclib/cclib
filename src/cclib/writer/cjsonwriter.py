# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2014-2016, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.
from scipy.constants.constants import alpha

"""A writer for chemical JSON (CJSON) files."""

try:
    import openbabel as ob
    has_openbabel = True
except ImportError:
    has_openbabel = False

import os.path
import json
import base64
import numpy as np

from . import filewriter


class CJSON(filewriter.Writer):
    """A writer for chemical JSON (CJSON) files."""
    
    # The expected Key names for all supported attributes.
    _attrkeynames = {
        "aonames":        'names',
        "aooverlaps":     'overlaps',
        "atombasis":      'indices',
        "atomcharges":    'dict',
        "atomcoords":     'coords',
        "atommasses":     'mass',
        "atomnos":        'number',
        "atomspins":      'spins',
        "ccenergies":     'coupledCluster',
        "charge":         'charge',
        "coreelectrons":  'coreElectrons',
        "enthalpy":       'enthalpy',
        "entropy":        'entropy',
        "etenergies":     'electronicTransitions',
        "etoscs":         'oscillatorStrength',
        "etrotats":       'rotatoryStrength',
        "etsecs":         'oneExcitedConfig',
        "etsyms":         'symmetry',
        "freeenergy":     'freeEnergy',
        "fonames":        'orbitalNames',
        "fooverlaps":     'orbitalOverlap',
        "fragnames":      'names',
        "frags":          'atomIndices',
        'gbasis':         'TBD',
        "geotargets":     'geometricTargets',
        "geovalues":      'geometricValues',
        "grads":          'TBD',
        "hessian":        'hessianMatrix',
        "homos":          'homos',
        "mocoeffs":       'coeffs',
        "moenergies":     'energies',
        "moments":        'totalDipoleMoment',
        "mosyms":         'symmetry',
        "mpenergies":     'mollerPlesset',
        "mult":           'multiplicity',
        "natom":          'numberOfAtoms',
        "nbasis":         'basisNumber',
        "nmo":            'MO-number',
        "nocoeffs":       'TBD',
        "nooccnos":       'TBD',
        "optdone":        'done',
        "optstatus":      'Status',
        "scancoords":     'stepGeometry',
        "scanenergies":   'PES-energies',
        "scannames":      'variableNames',
        "scanparm":       'PES-parameterValues',
        "scfenergies":    'energies',
        "scftargets":     'targets',
        "scfvalues":      'values',
        "temperature":    'temperature',
        "vibanharms":     'anharmonicityConstants',
        "vibdisps":       'displacement',
        "vibfreqs":       'frequencies',
        "vibirs":         'IR',
        "vibramans":      'raman',
        "vibsyms":        'symmetry',
    }

    def __init__(self, ccdata, *args, **kwargs):
        """Initialize the chemical JSON writer object.

        Inputs:
          ccdata - An instance of ccData, parsed from a logfile.
        """

        # Call the __init__ method of the superclass
        super(CJSON, self).__init__(ccdata, *args, **kwargs)

        self.generate_repr()
    
    
    def generate_repr(self):
        """Generate the CJSON representation of the logfile data.
           Naming Convention followed:
              Dictionary object have names first letter Capitalized
              Attribute names have all lower characters
        """

        cjson_dict = dict()
        cjson_dict['chemical json'] = 0

        if self.jobfilename is not None:
            cjson_dict['name'] = os.path.splitext(self.jobfilename)[0]

        # These are properties that can be collected using Open Babel.
        if has_openbabel:
            cjson_dict['smiles'] = self.pbmol.write('smiles')
            cjson_dict['inchi'] = self.pbmol.write('inchi')
            cjson_dict['inchikey'] = self.pbmol.write('inchikey')
            cjson_dict['formula'] = self.pbmol.formula
        
        #Helpers functions which use properties provided by cclib      
        self.generateProperties(cjson_dict)
        self.generateAtoms(cjson_dict)
        self.generateOptimization(cjson_dict)
        self.generateVibrations(cjson_dict)
        self.generateBonds(cjson_dict)
        self.generateTransitions(cjson_dict)
        self.generateFragments(cjson_dict)
     
        if has_openbabel:
            cjson_dict['diagram'] = self.pbmol.write(format='svg')
            
                 
        return json.dumps(cjson_dict, cls=NumpyEncoder)

    
    def hasData(self, attrNames):
        """Returns true if any of the attributes 
           exist in ccData
        """ 
        for name in attrNames:
            if name == 'moenergies':
                if (hasattr(self.ccdata, 'moenergies') and hasattr(self.ccdata, 'homos')):
                    return True
            elif hasattr(self.ccdata, name):
                return True
        return False
    
    def getAttrList(self, attrList):
        """Returns a list of attributes which can 
           be parsed by cclib
        """
        list = []
        for name in attrList:
            if hasattr(self.ccdata, name):
                list.append(name)
        return list
        
        
            
    def generateProperties(self, cjson_dict):
        """ Appends the Properties object into the cjson
        Properties table:
            1) Molecular Mass
            2) Charge
            3) Multiplicity
            4) Energy
                 i) alpha
                     a) Homo  
                     b) Gap   
                 ii) Beta        
                     a) Homo        
                     b) Gap         
                iii) Total       
                 iv) Free Energy 
                  v) Moller - Plesset 
                 vi) Coupled Cluster 
            5) Enthalpy 
            6) Entropy 
            7) numberOfAtoms    
            8) Temperature 
            9) totalDipoleMoment 
            10) Partial Charges
                 i) Mulliken 
            11) Orbitals 
                 i) Homos   
                ii) Energies  
               iii) Overlaps 
                iv) Symmetry 
                 v) Coeffs
        """
        
        cjson_dict['Properties'] = dict()
        
        if has_openbabel:
            cjson_dict['Properties']['molecularMass'] = self.pbmol.molwt
            
        if hasattr(self.ccdata, 'charge'):
            cjson_dict['Properties']['charge'] = self.ccdata.charge
        if hasattr(self.ccdata, 'mult'):
            cjson_dict['Properties']['mult'] = self.ccdata.mult
        
        energyAttr = ['moenergies', 'freeenergy', 'mpenergies', 'ccenergies' ]
        if self.hasData(energyAttr):
            cjson_dict['Properties']['Energy'] = dict()
            
            if (hasattr(self.ccdata, 'moenergies') and hasattr(self.ccdata, 'homos')):            
                cjson_dict['Properties']['Energy']['Alpha'] = dict()
                cjson_dict['Properties']['Energy']['Beta'] = dict()
                
                homo_idx_alpha = int(self.ccdata.homos[0])
                homo_idx_beta = int(self.ccdata.homos[-1])
                energy_alpha_homo = self.ccdata.moenergies[0][homo_idx_alpha]
                energy_alpha_lumo = self.ccdata.moenergies[0][homo_idx_alpha + 1]
                energy_alpha_gap = energy_alpha_lumo - energy_alpha_homo
                energy_beta_homo = self.ccdata.moenergies[-1][homo_idx_beta]
                energy_beta_lumo = self.ccdata.moenergies[-1][homo_idx_beta + 1]
                energy_beta_gap = energy_beta_lumo - energy_beta_homo
                
                cjson_dict['Properties']['Energy']['Alpha']['homo'] = energy_alpha_homo
                cjson_dict['Properties']['Energy']['Alpha']['gap'] = energy_alpha_gap
                cjson_dict['Properties']['Energy']['Beta']['homo'] = energy_beta_homo
                cjson_dict['Properties']['Energy']['Beta']['gap'] = energy_beta_gap
                cjson_dict['Properties']['Energy']['total'] = self.ccdata.scfenergies[-1]
            
            if hasattr(self.ccdata, 'freeenergy'):
                cjson_dict['Properties']['Energy']['freeEnergy'] = self.ccdata.freeenergy
            
            # Check if we need to pass the entire ndarray or can we just pass the last value of the nD array
            if hasattr(self.ccdata, 'mpenergies'):
                cjson_dict['Properties']['Energy']['mollerPlesset'] = self.ccdata.mpenergies
            
            #Same point as above
            if hasattr(self.ccdata, 'ccenergies'):
                cjson_dict['Properties']['Energy']['coupledCluster'] = self.ccdata.mpenergies
        
        if hasattr(self.ccdata, 'natom'):
            cjson_dict['Properties']['number of atoms'] = self.ccdata.natom
            
        propertyAttr = ['enthalpy', 'entropy','temperature', 'atomcharges']
        existAttr = self.getAttrList(propertyAttr)
        
        if hasattr(self.ccdata, 'moments'):
            cjson_dict['totalDipoleMoment'] = self._calculate_total_dipole_moment()
        
        for attribute in existAttr:
            cjson_dict['Properties'][attribute] = getattr(self.ccdata, attribute)
        
        orbitalAttr = ['homos', 'moenergies', 'aooverlaps', 'mosyms', 'mocoeffs']
        if self.hasData(orbitalAttr):
            cjson_dict['Properties']['Orbitals'] = dict()
            
            # I will be sacrificing descriptive attribute names for code brevity - Check
            existOrbitalAttr = self.getAttrList(orbitalAttr)
            for attribute in existOrbitalAttr:
                cjson_dict['Properties'][attribute] = getattr(self.ccdata, attribute)
                
                
    def generateAtoms(self,cjson_dict):
        """ Appends the Atoms object into the cjson
        Atoms Table:
            1) Elements
                a) Number
                b) atomCount                
                c) heavyAtomCount                
            2) Coords
                a) 3d                    
            3) Orbitals
                a) Names
                b) Indices
            4) Coreelectrons
            5) Mass
            6) Spins
        """
        
        cjson_dict['Atoms'] = dict()
        
        if hasattr(self.ccdata, 'atomnos'):
            cjson_dict['Atoms']['Elements'] = dict()
            cjson_dict['Atoms']['Elements']['atomicNumber'] = self.ccdata.atomnos.tolist()
            cjson_dict['Atoms']['Elements']['atomCount'] = len(self.ccdata.atomnos)
            cjson_dict['Atoms']['Elements']['heavyAtomCount'] = len([x for x in self.ccdata.atomnos if x > 1])
        
        if hasattr(self.ccdata, 'atomcoords'):
            cjson_dict['Atoms']['Coords'] = dict()
            cjson_dict['Atoms']['Coords']['3d'] = self.ccdata.atomcoords[-1].flatten().tolist()
            
        orbitalList = ['aonames', 'atombasis']
        if self.hasData(orbitalList):
            cjson_dict['Atoms']['Orbitals'] = dict()
            if hasattr(self.ccdata, 'aonames'):
                cjson_dict['Atoms']['Orbitals']['names'] = self.ccdata.aonames
            if hasattr(self.ccdata, 'atombasis'):
                cjson_dict['Atoms']['Orbitals']['Indices'] = self.ccdata.atombasis
                
        if hasattr(self.ccdata, 'coreelectrons'):
            cjson_dict['Atoms']['coreElectrons'] = self.ccdata.coreelectrons.tolist()
            
        if hasattr(self.ccdata, 'atommasses'):
            cjson_dict['Atoms']['atomicMasses'] = self.ccdata.atommasses.tolist()
            
        if hasattr(self.ccdata, 'atomspins'):
            cjson_dict['Atoms']['atomicSpins'] = self.ccdata.atomspins
            
        
    def generateOptimization(self,cjson_dict):
        """ Appends the Optimization object into the cjson
            Optimization table:
                1) Done 
                2) Status  
                3) Geometric Targets 
                4) Geometric Values 
                5) Basis number 
                6) MO number 
                7) SCF 
                    a) Energies 
                    b) Targets 
                    c) Values 
                8) Scan 
                    a) Step Geometry 
                    b) Potential Energy Surface - energies     
                    c) Variable names 
                    d) PES Parameter Values 
        
        """
        optiAttr = ['optdone', 'geotargets','nbasis', 'nmo', 'scfenergies', 'scancoords', 'scannames']
        if self.hasData(optiAttr):
            cjson_dict['Optimization'] = dict()
            attrList = ['optdone', 'optstatus', 'geotargets', 'geovalues', 'nbasis', 'nmo' ]
            existAttr = self.getAttrList(attrList)
            
            for attribute in existAttr:
                cjson_dict['Optimization'][attribute] = getattr(self.ccdata, attribute)
            
            #assumption: If SCFenergies exist, then scftargets will also exist
            if hasattr(self.ccdata, 'scfenergies'):
                cjson_dict['Optimization']['SCF'] = dict()
                cjson_dict['Optimization']['SCF']['energies'] = self.ccdata.scfenergies
                cjson_dict['Optimization']['SCF']['targets'] = self.ccdata.scftargets
            if hasattr(self.ccdata, 'scfvalues'):
                cjson_dict['Optimization']['SCF']['values'] = self.ccdata.scfvalues
            
            #Same assumption as above
            if hasattr(self.ccdata, 'scanenergies'):
                cjson_dict['Optimization']['Scan'] = dict()
                cjson_dict['Optimization']['Scan']['stepGeometry'] = self.ccdata.scancoords.scancoords
                cjson_dict['Optimization']['Scan']['PES-energies'] = self.ccdata.scancoords.scanenergies
                cjson_dict['Optimization']['Scan']['variableNames'] = self.ccdata.scancoords.scannames
                cjson_dict['Optimization']['Scan']['PES-params'] = self.ccdata.scancoords.scanparm
                
             
    def generateVibrations(self,cjson_dict):
        """ Appends the Vibrations object into the cjson
            Vibrations table:
                1) Anharmonicity constants 
                2) Frequencies 
                3) Intensities 
                    a)IR 
                    b) Raman 
                4) Symmetry 
                5) Hessian matrix 
                6) Displacement               
        """
        vibAttr = ['vibanharms', 'vibanharms', 'vibirs', 'vibramans', 'vibsyms', 'hessian', 'vibdisps' ]
        if self.hasData(vibAttr):
            cjson_dict['Vibrations'] = dict()
            
            if hasattr(self.ccdata, 'vibanharms'):
                cjson_dict['Vibrations']['anharmonicityConstants'] = self.ccdata.vibanharms
            
            if hasattr(self.ccdata, 'vibfreqs'):
                cjson_dict['Vibrations']['frequencies'] = self.ccdata.vibfreqs
            
            if hasattr(self.ccdata, 'vibirs') or hasattr(self.ccdata, 'vibramans'):
                cjson_dict['Vibrations']['Intensities'] = dict()
                if hasattr(self.ccdata, 'vibirs'):
                    cjson_dict['Vibrations']['Intensities']['IR'] = self.ccdata.vibirs
                else:
                    cjson_dict['Vibrations']['Intensities']['Raman'] = self.ccdata.vibramans
                
            attrList = ['vibsyms', 'hessian', 'vibdisps']
            availableAttrList = self.getAttrList(attrList)
            
            for attribute in availableAttrList:
                cjson_dict['Vibrations'][attribute] = getattr(self.ccdata, attribute)
            
    def generateBonds(self,cjson_dict):
        """ Appends the Bonds object into the cjson
            Bonds table:
                1) Connections 
                    a) Index
                2) Order    
        """
        if has_openbabel:
            cjson_dict['Bonds'] = dict()
            cjson_dict['Bonds']['Connections'] = dict()
            cjson_dict['Bonds']['Connections']['index'] = []
            for bond in self.bond_connectivities:
                cjson_dict['bonds']['connections']['index'].append(bond[0] + 1)
                cjson_dict['bonds']['connections']['index'].append(bond[1] + 1)
            cjson_dict['bonds']['order'] = [bond[2] for bond in self.bond_connectivities]
            
    def generateTransitions(self,cjson_dict):
        """ Appends the Transition object into the cjson
            Transitions table:
               1) Electronic Transitions 
               2) Oscillator Strength 
               3) Rotatory Strength 
               4) 1-excited-config  
               5) Symmetry 
        """
        attrList = ['etenergies', 'etoscs', 'etrotats', 'etsecs', 'etsyms' ]
        if self.hasData(attrList):
            cjson_dict['Transitions'] = dict()
            availableAttr = self.getAttrList(attrList)
            
            for attribute in availableAttr:
                cjson_dict['Transitions'][attribute] = getattr(self.ccdata, attribute)
                
    def generateFragments(self,cjson_dict):
        """ Appends the Fragments object into the cjson
            Fragments table:
               1) Names 
               2) Atom Indices 
               3) Orbital Names 
               4) Orbital Overlap 
        """
        attrList = ['fragnames', 'frags', 'fonames', 'fooverlaps' ]
        if self.hasData(attrList):
            cjson_dict['Fragments'] = dict()
            availableAttr = self.getAttrList(attrList)
            
            for attribute in availableAttr:
                cjson_dict['Fragments'][attribute] = getattr(self.ccdata, attribute)
    

class NumpyEncoder(json.JSONEncoder):

    def default(self, obj):
        """If input object is an ndarray it will be converted into a dict 
        holding dtype, shape and the data, base64 encoded.
        """
        if isinstance(obj, np.ndarray):
            if obj.flags['C_CONTIGUOUS']:
                obj_data = obj.data
            else:
                cont_obj = np.ascontiguousarray(obj)
                assert(cont_obj.flags['C_CONTIGUOUS'])
                obj_data = cont_obj.data
            data_b64 = base64.b64encode(obj_data)
            return dict(__ndarray__=data_b64,
                        dtype=str(obj.dtype),
                        shape=obj.shape)
        # Let the base class default method raise the TypeError
        return json.JSONEncoder(self, obj)


    def json_numpy_obj_hook(dct):
        """Decodes a previously encoded numpy ndarray with proper shape and dtype.
    
        :param dct: (dict) json encoded ndarray
        :return: (ndarray) if input was an encoded ndarray
        """
        if isinstance(dct, dict) and '__ndarray__' in dct:
            data = base64.b64decode(dct['__ndarray__'])
            return np.frombuffer(data, dct['dtype']).reshape(dct['shape'])
        return dct
                
if __name__ == "__main__":
    pass
