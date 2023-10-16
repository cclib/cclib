# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

# import sys
# import inspect
from abc import ABC
from collections import namedtuple

import numpy as np


class _Property(ABC):

    def __init__(self, name: str, main_type, json_key: str, json_path: str, *args, **kwargs) -> None:

        self.name = __name__
        self.main_type = main_type
        self.json_key = json_key
        self.json_path = json_path

        # backwards compatibility
        self.type = main_type
        self.attribute_path = json_path

        self._impl = None

    def __call__(self, val=None):
        if val:
            self._impl = val
            self.typecheck()
        else:
            if not self._impl:
                raise PropertyError("{} doesn't contain any data".format(type(self)))
            else:
                return self._impl

    def __str__(self):
        return str(self.__class__.__name__)

    def typecheck(self):
        """Check the types of all attributes.

        If an attribute does not match the expected type, then attempt to
        convert; if that fails, only then raise a TypeError.
        """

        if type(self._impl) == self.main_type:
            return

        try:
            if self.main_type == np.ndarray:
                self._impl = np.asarray(self._impl)
            else:
                self._impl = self.main_type(self._impl)
        except ValueError:
            args = (self.name, type(self._impl), self.main_type)
            raise TypeError("attribute {} is {} instead of {} and could not be converted".format(*args))


class aonames(_Property):
    """atomic orbital names (list of strings)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('aonames', list, 'names', 'atoms:orbitals',
                                         *args, **kwargs)


class aooverlaps(_Property):
    """atomic orbital overlap matrix (array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('aooverlaps', np.ndarray, 'overlaps', 'properties:orbitals',
                                         *args, **kwargs)


class atombasis(_Property):
    """indices of atomic orbitals on each atom (list of lists)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('atombasis', list, 'indices', 'atoms:orbitals',
                                         *args, **kwargs)


class atomcharges(_Property):
    """atomic partial charges (dict of arrays[1])"""

    # Propertys that should be dictionaries of arrays (double precision).
    _dictsofarrays = ["atomcharges", "atomspins"]
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('atomcharges', dict, 'partial charges', 'properties',
                                         *args, **kwargs)


class atomcoords(_Property):
    """atom coordinates (array[3], angstroms)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('atomcoords', np.ndarray, 'coords', 'atoms:coords:3d',
                                         *args, **kwargs)


class atommasses(_Property):
    """atom masses (array[1], daltons)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('atommasses', np.ndarray, 'mass', 'atoms',
                                         *args, **kwargs)


class atomnos(_Property):
    """atomic numbers (array[1])"""

    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ['atomnos', 'coreelectrons', 'homos', 'optstatus']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('atomnos', np.ndarray, 'number', 'atoms:elements',
                                         *args, **kwargs)


class atomspins(_Property):
    """atomic spin densities (dict of arrays[1])"""

    # Propertys that should be dictionaries of arrays (double precision).
    _dictsofarrays = ["atomcharges", "atomspins"]
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('atomspins', dict, 'spins', 'atoms',
                                         *args, **kwargs)


class ccenergies(_Property):
    """molecular energies with Coupled-Cluster corrections (array[2], eV)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('ccenergies', np.ndarray, 'coupled cluster', 'properties:energy',
                                         *args, **kwargs)


class charge(_Property):
    """net charge of the system (integer)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('charge', int, 'charge', 'properties',
                                         *args, **kwargs)


class coreelectrons(_Property):
    """number of core electrons in atom pseudopotentials (array[1])"""

    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ['atomnos', 'coreelectrons', 'homos', 'optstatus']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('coreelectrons', np.ndarray, 'core electrons', 'atoms',
                                         *args, **kwargs)


class dispersionenergies(_Property):
    """a molecular dispersion energy corrections (array[1], eV)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('dispersionenergies', np.ndarray, 'dispersion correction', 'properties:energy',
                                         *args, **kwargs)

class enthalpy(_Property):
    """sum of electronic and thermal enthalpies (float, hartree/particle)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('enthalpy', float, 'enthalpy', 'properties',
                                         *args, **kwargs)


class entropy(_Property):
    """entropy (float, hartree/particle)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('entropy', float, 'entropy', 'properties',
                                         *args, **kwargs)


class etenergies(_Property):
    """energies of electronic transitions (array[1], 1/cm)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('etenergies', np.ndarray, 'electronic transitions', 'transitions',
                                         *args, **kwargs)


class etoscs(_Property):
    """oscillator strengths of electronic transitions (array[1])"""
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('etoscs', np.ndarray, 'oscillator strength', 'transitions',
                                         *args, **kwargs)


class etdips(_Property):
    """electric transition dipoles of electronic transitions (array[2], ebohr)"""
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('etdips', np.ndarray, 'electric transition dipoles', 'transitions',
                                         *args, **kwargs)


class etveldips(_Property):
    """velocity-gauge electric transition dipoles of electronic transitions (array[2], ebohr)"""
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('etveldips', np.ndarray, 'velocity-gauge electric transition dipoles', 'transitions',
                                         *args, **kwargs)


class etmagdips(_Property):
    """magnetic transition dipoles of electronic transitions (array[2], ebohr)"""
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('etmagdips', np.ndarray, 'magnetic transition dipoles', 'transitions',
                                         *args, **kwargs)


class etrotats(_Property):
    """rotatory strengths of electronic transitions (array[1], ??)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('etrotats', np.ndarray, 'rotatory strength', 'transitions',
                                         *args, **kwargs)


class etsecs(_Property):
    """singly-excited configurations for electronic transitions (list of lists)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('etsecs', list, 'one excited config', 'transitions',
                                         *args, **kwargs)


class etsyms(_Property):
    """symmetries of electronic transitions (list of string)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('etsyms', list, 'symmetry', 'transitions',
                                         *args, **kwargs)


class freeenergy(_Property):
    """sum of electronic and thermal free energies (float, hartree/particle)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('freeenergy', float, 'free energy', 'properties:energy',
                                         *args, **kwargs)


class fonames(_Property):
    """fragment orbital names (list of strings)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('fonames', list, 'orbital names', 'fragments',
                                         *args, **kwargs)


class fooverlaps(_Property):
    """fragment orbital overlap matrix (array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('fooverlaps', np.ndarray, 'orbital overlap', 'fragments',
                                         *args, **kwargs)


class fragnames(_Property):
    """names of fragments (list of strings)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('fragnames', list, 'fragment names', 'fragments',
                                         *args, **kwargs)


class frags(_Property):
    """indices of atoms in a fragment (list of lists)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('frags', list, 'atom indices', 'fragments',
                                         *args, **kwargs)


class gbasis(_Property):
    """coefficients and exponents of Gaussian basis functions (PyQuante format)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('gbasis', list, 'basis functions', 'atoms:orbitals',
                                         *args, **kwargs)


class geotargets(_Property):
    """targets for convergence of geometry optimization (array[1])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('geotargets', np.ndarray, 'geometric targets', 'optimization',
                                         *args, **kwargs)


class geovalues(_Property):
    """current values for convergence of geometry optmization (array[1])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('geovalues', np.ndarray, 'geometric values', 'optimization',
                                         *args, **kwargs)


class grads(_Property):
    """current values of forces (gradients) in geometry optimization (array[3])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('grads', np.ndarray, 'TBD', 'N/A',
                                         *args, **kwargs)


class hessian(_Property):
    """elements of the force constant matrix (array[1])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('hessian', np.ndarray, 'hessian matrix', 'vibrations',
                                         *args, **kwargs)


class homos(_Property):
    """molecular orbital indices of HOMO(s) (array[1])"""

    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ['atomnos', 'coreelectrons', 'homos', 'optstatus']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('homos', np.ndarray, 'homos', 'properties:orbitals',
                                         *args, **kwargs)


class metadata(_Property):
    """various metadata about the package and computation (dict)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('metadata', dict, 'TBD', 'N/A',
                                         *args, **kwargs)


class mocoeffs(_Property):
    """molecular orbital coefficients (list of arrays[2])"""

    # Propertys that should be lists of arrays (double precision).
    _listsofarrays = ['mocoeffs', 'moenergies', 'moments', 'polarizabilities', 'scfvalues']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('mocoeffs', list, 'coeffs', 'properties:orbitals',
                                         *args, **kwargs)


class moenergies(_Property):
    """molecular orbital energies (list of arrays[1], eV)"""

    # Propertys that should be lists of arrays (double precision).
    _listsofarrays = ['mocoeffs', 'moenergies', 'moments', 'polarizabilities', 'scfvalues']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('moenergies', list, 'energies', 'properties:orbitals',
                                         *args, **kwargs)


class moments(_Property):
    """molecular multipole moments (list of arrays[], a.u.)"""

    # Propertys that should be lists of arrays (double precision).
    _listsofarrays = ['mocoeffs', 'moenergies', 'moments', 'polarizabilities', 'scfvalues']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('moments', list, 'total dipole moment', 'properties',
                                         *args, **kwargs)


class mosyms(_Property):
    """orbital symmetries (list of lists)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('mosyms', list, 'molecular orbital symmetry', 'properties:orbitals',
                                         *args, **kwargs)


class mpenergies(_Property):
    """molecular electronic energies with Møller-Plesset corrections (array[2], eV)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('mpenergies', np.ndarray, 'moller plesset', 'properties:energy',
                                         *args, **kwargs)


class mult(_Property):
    """multiplicity of the system (integer)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('mult', int, 'multiplicity', 'properties',
                                         *args, **kwargs)


class natom(_Property):
    """number of atoms (integer)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('natom', int, 'number of atoms', 'properties',
                                         *args, **kwargs)


class nbasis(_Property):
    """number of basis functions (integer)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('nbasis', int, 'basis number', 'properties:orbitals',
                                         *args, **kwargs)


class nmo(_Property):
    """number of molecular orbitals (integer)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('nmo', int, 'MO number', 'properties:orbitals',
                                         *args, **kwargs)


class nmrtensors(_Property):
    """Nuclear magnetic resonance chemical shielding tensors (dict of dicts of array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('nmr', int, 'NMR chemical shielding tensors', 'properties:nmr',
                                         *args, **kwargs)


class nocoeffs(_Property):
    """natural orbital coefficients (array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('nocoeffs', np.ndarray, 'TBD', 'N/A',
                                         *args, **kwargs)


class nooccnos(_Property):
    """natural orbital occupation numbers (array[1])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('nooccnos', np.ndarray, 'TBD', 'N/A',
                                         *args, **kwargs)


class nsocoeffs(_Property):
    """natural spin orbital coefficients (list of array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('nsocoeffs', list, 'TBD', 'N/A',
                                         *args, **kwargs)


class nsooccnos(_Property):
    """natural spin orbital occupation numbers (list of array[1])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('nsooccnos', list, 'TBD', 'N/A',
                                         *args, **kwargs)


class optdone(_Property):
    """flags whether an optimization has converged (Boolean)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('optdone', bool, 'done', 'optimization',
                                         *args, **kwargs)


class optstatus(_Property):
    """optimization status for each set of atomic coordinates (array[1])"""

    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ['atomnos', 'coreelectrons', 'homos', 'optstatus']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('optstatus', np.ndarray, 'status', 'optimization',
                                         *args, **kwargs)


class polarizabilities(_Property):
    """(dipole) polarizabilities, static or dynamic (list of arrays[2])"""

    # Propertys that should be lists of arrays (double precision).
    _listsofarrays = ['mocoeffs', 'moenergies', 'moments', 'polarizabilities', 'scfvalues']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('polarizabilities', list, 'polarizabilities', 'N/A',
                                         *args, **kwargs)


class pressure(_Property):
    """temperature used for Thermochemistry (float, atm)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('pressure', float, 'pressure', 'properties',
                                         *args, **kwargs)


class rotconsts(_Property):
    """rotational constants (array[2], GHz)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('rotconsts', np.ndarray, 'rotational constants', 'atoms:coords:rotconsts',
                                         *args, **kwargs)


class scancoords(_Property):
    """geometries of each scan step (array[3], angstroms)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('scancoords', np.ndarray, 'step geometry', 'optimization:scan',
                                         *args, **kwargs)


class scanenergies(_Property):
    """energies of potential energy surface (list)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('scanenergies', list, 'PES energies', 'optimization:scan',
                                         *args, **kwargs)


class scannames(_Property):
    """names of varaibles scanned (list of strings)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('scannames', list, 'variable names', 'optimization:scan',
                                         *args, **kwargs)


class scanparm(_Property):
    """values of parameters in potential energy surface (list of tuples)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('scanparm', list, 'PES parameter values', 'optimization:scan',
                                         *args, **kwargs)


class scfenergies(_Property):
    """molecular electronic energies after SCF (Hartree-Fock, DFT) (array[1], eV)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('scfenergies', np.ndarray, 'scf energies', 'optimization:scf',
                                         *args, **kwargs)


class scftargets(_Property):
    """targets for convergence of the SCF (array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('scftargets', np.ndarray, 'targets', 'optimization:scf',
                                         *args, **kwargs)


class scfvalues(_Property):
    """current values for convergence of the SCF (list of arrays[2])"""

    # Propertys that should be lists of arrays (double precision).
    _listsofarrays = ['mocoeffs', 'moenergies', 'moments', 'polarizabilities', 'scfvalues']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('scfvalues', list, 'values', 'optimization:scf',
                                         *args, **kwargs)


class temperature(_Property):
    """temperature used for Thermochemistry (float, kelvin)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('temperature', float, 'temperature', 'properties',
                                         *args, **kwargs)


class time(_Property):
    """time in molecular dynamics and other trajectories (array[1], fs)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('time', np.ndarray, 'time', 'N/A',
                                         *args, **kwargs)


class transprop(_Property):
    """all absorption and emission spectra (dictionary {name:(etenergies, etoscs)})

    WARNING: this attribute is not standardized and is liable to change in cclib 2.0
    """

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('transprop', dict, 'electronic transitions', 'transitions',
                                         *args, **kwargs)


class vibanharms(_Property):
    """vibrational anharmonicity constants (array[2], 1/cm)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('vibanharms', np.ndarray, 'anharmonicity constants', 'vibrations',
                                         *args, **kwargs)


class vibdisps(_Property):
    """cartesian displacement vectors (array[3], delta angstrom)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('vibdisps', np.ndarray, 'displacement', 'vibrations',
                                         *args, **kwargs)


class vibfreqs(_Property):
    """vibrational frequencies (array[1], 1/cm)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('vibfreqs', np.ndarray, 'frequencies', 'vibrations',
                                         *args, **kwargs)


class vibirs(_Property):
    """IR intensities (array[1], km/mol)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('vibirs', np.ndarray, 'IR', 'vibrations:intensities',
                                         *args, **kwargs)


class vibramans(_Property):
    """Raman intensities (array[1], A^4/Da)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('vibramans', np.ndarray, 'raman', 'vibrations:intensities',
                                         *args, **kwargs)


class vibrmasses(_Property):
    """reduced masses of vibrations (array[1], daltons)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('vibrmasses', np.ndarray, 'reduced masses', 'vibrations',
                                         *args, **kwargs)


class vibsyms(_Property):
    """symmetries of vibrations (list of strings)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('vibsyms', list, 'vibration symmetry', 'vibrations',
                                         *args, **kwargs)


class zpve(_Property):
    """zero-point vibrational energy correction (float, hartree/particle)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('zpve', float, 'zero-point correction', 'properties:energies',
                                         *args, **kwargs)


# https://stackoverflow.com/questions/1796180/how-can-i-get-a-list-of-all-classes-within-current-module-in-python#1796247
# attributes = dict()
# for name, obj in inspect.getmembers(sys.modules[__name__]):
#     if inspect.isclass(obj):
#         if name != '_Property':
#             attributes[name] = obj

Property = namedtuple('Property', ['type', 'json_key', 'attribute_path'])

# The expected types for all supported attributes.
# The json_key is the key name used for attributes in the CJSON/JSON format
# 'TBD' - To Be Decided are the key names of attributes which haven't been included in the cjson format
_properties = {
   "aonames":            Property(list,             'names',                       'atoms:orbitals'),
   "aooverlaps":         Property(np.ndarray,       'overlaps',                    'properties:orbitals'),
   "atombasis":          Property(list,             'indices',                     'atoms:orbitals'),
   "atomcharges":        Property(dict,             'partial charges',             'properties'),
   "atomcoords":         Property(np.ndarray,       'coords',                      'atoms:coords:3d'),
   "atommasses":         Property(np.ndarray,       'mass',                        'atoms'),
   "atomnos":            Property(np.ndarray,       'number',                      'atoms:elements'),
   "atomspins":          Property(dict,             'spins',                       'atoms'),
   "ccenergies":         Property(np.ndarray,       'coupled cluster',             'properties:energy'),
   "charge":             Property(int,              'charge',                      'properties'),
   "coreelectrons":      Property(np.ndarray,       'core electrons',              'atoms'),
   "dispersionenergies": Property(np.ndarray,       'dispersion correction',       'properties:energy'),
   "enthalpy":           Property(float,            'enthalpy',                    'properties'),
   "entropy":            Property(float,            'entropy',                     'properties'),
   "etenergies":         Property(np.ndarray,       'electronic transitions',      'transitions'),
   "etoscs":             Property(np.ndarray,       'oscillator strength',         'transitions'),
   "etdips":             Property(np.ndarray,       'electic transition dipoles',  'transitions'),
   "etveldips":          Property(np.ndarray,       'velocity-gauge electric transition dipoles', 'transitions'),
   "etmagdips":          Property(np.ndarray,       'magnetic transition dipoles', 'transitions'),
   "etrotats":           Property(np.ndarray,       'rotatory strength',           'transitions'),
   "etsecs":             Property(list,             'one excited config',          'transitions'),
   "etsyms":             Property(list,             'symmetry',                    'transitions'),
   "freeenergy":         Property(float,            'free energy',                 'properties:energy'),
   "fonames":            Property(list,             'orbital names',               'fragments'),
   "fooverlaps":         Property(np.ndarray,       'orbital overlap',             'fragments'),
   "fragnames":          Property(list,             'fragment names',              'fragments'),
   "frags":              Property(list,             'atom indices',                'fragments'),
   "gbasis":             Property(list,             'basis functions',             'atoms:orbitals'),
   "geotargets":         Property(np.ndarray,       'geometric targets',           'optimization'),
   "geovalues":          Property(np.ndarray,       'geometric values',            'optimization'),
   "grads":              Property(np.ndarray,       'TBD',                         'N/A'),
   "hessian":            Property(np.ndarray,       'hessian matrix',              'vibrations'),
   "homos":              Property(np.ndarray,       'homos',                       'properties:orbitals'),
   "metadata":           Property(dict,             'TBD',                         'N/A'),
   "mocoeffs":           Property(list,             'coeffs',                      'properties:orbitals'),
   "moenergies":         Property(list,             'energies',                    'properties:orbitals'),
   "moments":            Property(list,             'total dipole moment',         'properties'),
   "mosyms":             Property(list,             'molecular orbital symmetry',  'properties:orbitals'),
   "mpenergies":         Property(np.ndarray,       'moller plesset',              'properties:energy'),
   "mult":               Property(int,              'multiplicity',                'properties'),
   "natom":              Property(int,              'number of atoms',             'properties'),
   "nbasis":             Property(int,              'basis number',                'properties:orbitals'),
   "nmo":                Property(int,              'MO number',                   'properties:orbitals'),
   "nmrtensors":         Property(dict,             'NMR chemical shielding tensors', 'properties:nmr'),
   "nocoeffs":           Property(np.ndarray,       'TBD',                         'N/A'),
   "nooccnos":           Property(np.ndarray,       'TBD',                         'N/A'),
   "nsocoeffs":          Property(list,             'TBD',                         'N/A'),
   "nsooccnos":          Property(list,             'TBD',                         'N/A'),
   "optdone":            Property(list,             'done',                        'optimization'),
   "optstatus":          Property(np.ndarray,       'status',                      'optimization'),
   "polarizabilities":   Property(list,             'polarizabilities',            'N/A'),
   "pressure":           Property(float,            'pressure',                    'properties'),
   "rotconsts":          Property(np.ndarray,       'rotational constants',        'atoms:coords:rotconsts'),
   "scancoords":         Property(np.ndarray,       'step geometry',               'optimization:scan'),
   "scanenergies":       Property(list,             'PES energies',                'optimization:scan'),
   "scannames":          Property(list,             'variable names',              'optimization:scan'),
   "scanparm":           Property(list,             'PES parameter values',        'optimization:scan'),
   "scfenergies":        Property(np.ndarray,       'scf energies',                'optimization:scf'),
   "scftargets":         Property(np.ndarray,       'targets',                     'optimization:scf'),
   "scfvalues":          Property(list,             'values',                      'optimization:scf'),
   "temperature":        Property(float,            'temperature',                 'properties'),
   "time":               Property(np.ndarray,       'time',                        'N/A'),
   "transprop":          Property(dict,             'electronic transitions',      'transitions'),
   "vibanharms":         Property(np.ndarray,       'anharmonicity constants',     'vibrations'),
   "vibdisps":           Property(np.ndarray,       'displacement',                'vibrations'),
   "vibfreqs":           Property(np.ndarray,       'frequencies',                 'vibrations'),
   "vibfconsts":         Property(np.ndarray,       'force constants',             'vibrations'),
   "vibirs":             Property(np.ndarray,       'IR',                          'vibrations:intensities'),
   "vibramans":          Property(np.ndarray,       'raman',                       'vibrations:intensities'),
   "vibrmasses":         Property(np.ndarray,       'reduced masses',              'vibrations'),
   "vibsyms":            Property(list,             'vibration symmetry',          'vibrations'),
   "zpve":               Property(float,            'zero-point correction',       'properties:energies')
}
