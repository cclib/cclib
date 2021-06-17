# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import sys
import inspect

import numpy as np


class Attribute:

    def __init__(self, name, main_type, json_key, json_path, *args, **kwargs):

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
                raise AttributeError("{} doesn't contain any data".format(type(self)))
            else:
                return self._impl

    def __str__(self):
        return str(self.__class__.__name__)

    def typecheck(self):
        """Check the types of all attributes.

        If an attribute does not match the expected type, then attempt to
        convert; if that fails, only then raise a TypeError.
        """

        if type(self._impl) == self.type:
            return

        try:
            self._impl = self.main_type(self._impl)
        except ValueError:
            args = (self.name, type(self._impl), self.main_type)
            raise TypeError("attribute {} is {} instead of {} and could not be converted".format(*args))


class aonames(Attribute):
    """atomic orbital names (list of strings)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('aonames', list, 'names', 'atoms:orbitals',
                                         *args, **kwargs)


class aooverlaps(Attribute):
    """atomic orbital overlap matrix (array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('aooverlaps', np.ndarray, 'overlaps', 'properties:orbitals',
                                         *args, **kwargs)


class atombasis(Attribute):
    """indices of atomic orbitals on each atom (list of lists)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('atombasis', list, 'indices', 'atoms:orbitals',
                                         *args, **kwargs)


class atomcharges(Attribute):
    """atomic partial charges (dict of arrays[1])"""

    # Attributes that should be dictionaries of arrays (double precision).
    _dictsofarrays = ["atomcharges", "atomspins"]
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('atomcharges', dict, 'partial charges', 'properties',
                                         *args, **kwargs)


class atomcoords(Attribute):
    """atom coordinates (array[3], angstroms)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('atomcoords', np.ndarray, 'coords', 'atoms:coords:3d',
                                         *args, **kwargs)


class atommasses(Attribute):
    """atom masses (array[1], daltons)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('atommasses', np.ndarray, 'mass', 'atoms',
                                         *args, **kwargs)


class atomnos(Attribute):
    """atomic numbers (array[1])"""

    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ['atomnos', 'coreelectrons', 'homos', 'optstatus']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('atomnos', np.ndarray, 'number', 'atoms:elements',
                                         *args, **kwargs)


class atomspins(Attribute):
    """atomic spin densities (dict of arrays[1])"""

    # Attributes that should be dictionaries of arrays (double precision).
    _dictsofarrays = ["atomcharges", "atomspins"]
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('atomspins', dict, 'spins', 'atoms',
                                         *args, **kwargs)


class ccenergies(Attribute):
    """molecular energies with Coupled-Cluster corrections (array[2], eV)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('ccenergies', np.ndarray, 'coupled cluster', 'properties:energy',
                                         *args, **kwargs)


class charge(Attribute):
    """net charge of the system (integer)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('charge', int, 'charge', 'properties',
                                         *args, **kwargs)


class coreelectrons(Attribute):
    """number of core electrons in atom pseudopotentials (array[1])"""

    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ['atomnos', 'coreelectrons', 'homos', 'optstatus']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('coreelectrons', np.ndarray, 'core electrons', 'atoms',
                                         *args, **kwargs)


class dispersionenergies(Attribute):
    """a molecular dispersion energy corrections (array[1], eV)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('dispersionenergies', np.ndarray, 'dispersion correction', 'properties:energy',
                                         *args, **kwargs)

class enthalpy(Attribute):
    """sum of electronic and thermal enthalpies (float, hartree/particle)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('enthalpy', float, 'enthalpy', 'properties',
                                         *args, **kwargs)


class entropy(Attribute):
    """entropy (float, hartree/particle)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('entropy', float, 'entropy', 'properties',
                                         *args, **kwargs)


class etenergies(Attribute):
    """energies of electronic transitions (array[1], 1/cm)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('etenergies', np.ndarray, 'electronic transitions', 'transitions',
                                         *args, **kwargs)


class etoscs(Attribute):
    """oscillator strengths of electronic transitions (array[1])"""
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('etoscs', np.ndarray, 'oscillator strength', 'transitions',
                                         *args, **kwargs)


class etdips(Attribute):
    """electric transition dipoles of electronic transitions (array[2], ebohr)"""
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('etdips', np.ndarray, 'electric transition dipoles', 'transitions',
                                         *args, **kwargs)


class etveldips(Attribute):
    """velocity-gauge electric transition dipoles of electronic transitions (array[2], ebohr)"""
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('etveldips', np.ndarray, 'velocity-gauge electric transition dipoles', 'transitions',
                                         *args, **kwargs)


class etmagdips(Attribute):
    """magnetic transition dipoles of electronic transitions (array[2], ebohr)"""
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('etmagdips', np.ndarray, 'magnetic transition dipoles', 'transitions',
                                         *args, **kwargs)


class etrotats(Attribute):
    """rotatory strengths of electronic transitions (array[1], ??)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('etrotats', np.ndarray, 'rotatory strength', 'transitions',
                                         *args, **kwargs)


class etsecs(Attribute):
    """singly-excited configurations for electronic transitions (list of lists)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('etsecs', list, 'one excited config', 'transitions',
                                         *args, **kwargs)


class etsyms(Attribute):
    """symmetries of electronic transitions (list of string)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('etsyms', list, 'symmetry', 'transitions',
                                         *args, **kwargs)


class freeenergy(Attribute):
    """sum of electronic and thermal free energies (float, hartree/particle)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('freeenergy', float, 'free energy', 'properties:energy',
                                         *args, **kwargs)


class fonames(Attribute):
    """fragment orbital names (list of strings)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('fonames', list, 'orbital names', 'fragments',
                                         *args, **kwargs)


class fooverlaps(Attribute):
    """fragment orbital overlap matrix (array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('fooverlaps', np.ndarray, 'orbital overlap', 'fragments',
                                         *args, **kwargs)


class fragnames(Attribute):
    """names of fragments (list of strings)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('fragnames', list, 'fragment names', 'fragments',
                                         *args, **kwargs)


class frags(Attribute):
    """indices of atoms in a fragment (list of lists)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('frags', list, 'atom indices', 'fragments',
                                         *args, **kwargs)


class gbasis(Attribute):
    """coefficients and exponents of Gaussian basis functions (PyQuante format)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('gbasis', list, 'basis functions', 'atoms:orbitals',
                                         *args, **kwargs)


class geotargets(Attribute):
    """targets for convergence of geometry optimization (array[1])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('geotargets', np.ndarray, 'geometric targets', 'optimization',
                                         *args, **kwargs)


class geovalues(Attribute):
    """current values for convergence of geometry optmization (array[1])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('geovalues', np.ndarray, 'geometric values', 'optimization',
                                         *args, **kwargs)


class grads(Attribute):
    """current values of forces (gradients) in geometry optimization (array[3])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('grads', np.ndarray, 'TBD', 'N/A',
                                         *args, **kwargs)


class hessian(Attribute):
    """elements of the force constant matrix (array[1])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('hessian', np.ndarray, 'hessian matrix', 'vibrations',
                                         *args, **kwargs)


class homos(Attribute):
    """molecular orbital indices of HOMO(s) (array[1])"""

    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ['atomnos', 'coreelectrons', 'homos', 'optstatus']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('homos', np.ndarray, 'homos', 'properties:orbitals',
                                         *args, **kwargs)


class metadata(Attribute):
    """various metadata about the package and computation (dict)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('metadata', dict, 'TBD', 'N/A',
                                         *args, **kwargs)


class mocoeffs(Attribute):
    """molecular orbital coefficients (list of arrays[2])"""

    # Attributes that should be lists of arrays (double precision).
    _listsofarrays = ['mocoeffs', 'moenergies', 'moments', 'polarizabilities', 'scfvalues']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('mocoeffs', list, 'coeffs', 'properties:orbitals',
                                         *args, **kwargs)


class moenergies(Attribute):
    """molecular orbital energies (list of arrays[1], eV)"""

    # Attributes that should be lists of arrays (double precision).
    _listsofarrays = ['mocoeffs', 'moenergies', 'moments', 'polarizabilities', 'scfvalues']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('moenergies', list, 'energies', 'properties:orbitals',
                                         *args, **kwargs)


class moments(Attribute):
    """molecular multipole moments (list of arrays[], a.u.)"""

    # Attributes that should be lists of arrays (double precision).
    _listsofarrays = ['mocoeffs', 'moenergies', 'moments', 'polarizabilities', 'scfvalues']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('moments', list, 'total dipole moment', 'properties',
                                         *args, **kwargs)


class mosyms(Attribute):
    """orbital symmetries (list of lists)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('mosyms', list, 'molecular orbital symmetry', 'properties:orbitals',
                                         *args, **kwargs)


class mpenergies(Attribute):
    """molecular electronic energies with Møller-Plesset corrections (array[2], eV)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('mpenergies', np.ndarray, 'moller plesset', 'properties:energy',
                                         *args, **kwargs)


class mult(Attribute):
    """multiplicity of the system (integer)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('mult', int, 'multiplicity', 'properties',
                                         *args, **kwargs)


class natom(Attribute):
    """number of atoms (integer)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('natom', int, 'number of atoms', 'properties',
                                         *args, **kwargs)


class nbasis(Attribute):
    """number of basis functions (integer)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('nbasis', int, 'basis number', 'properties:orbitals',
                                         *args, **kwargs)


class nmo(Attribute):
    """number of molecular orbitals (integer)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('nmo', int, 'MO number', 'properties:orbitals',
                                         *args, **kwargs)


class nmrtensors(Attribute):
    """Nuclear magnetic resonance chemical shielding tensors (dict of dicts of array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('nmr', int, 'NMR chemical shielding tensors', 'properties:nmr',
                                         *args, **kwargs)


class nocoeffs(Attribute):
    """natural orbital coefficients (array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('nocoeffs', np.ndarray, 'TBD', 'N/A',
                                         *args, **kwargs)


class nooccnos(Attribute):
    """natural orbital occupation numbers (array[1])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('nooccnos', np.ndarray, 'TBD', 'N/A',
                                         *args, **kwargs)


class nsocoeffs(Attribute):
    """natural spin orbital coefficients (list of array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('nsocoeffs', list, 'TBD', 'N/A',
                                         *args, **kwargs)


class nsooccnos(Attribute):
    """natural spin orbital occupation numbers (list of array[1])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('nsooccnos', list, 'TBD', 'N/A',
                                         *args, **kwargs)


class optdone(Attribute):
    """flags whether an optimization has converged (Boolean)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('optdone', bool, 'done', 'optimization',
                                         *args, **kwargs)


class optstatus(Attribute):
    """optimization status for each set of atomic coordinates (array[1])"""

    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ['atomnos', 'coreelectrons', 'homos', 'optstatus']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('optstatus', np.ndarray, 'status', 'optimization',
                                         *args, **kwargs)


class polarizabilities(Attribute):
    """(dipole) polarizabilities, static or dynamic (list of arrays[2])"""

    # Attributes that should be lists of arrays (double precision).
    _listsofarrays = ['mocoeffs', 'moenergies', 'moments', 'polarizabilities', 'scfvalues']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('polarizabilities', list, 'polarizabilities', 'N/A',
                                         *args, **kwargs)


class pressure(Attribute):
    """temperature used for Thermochemistry (float, atm)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('pressure', float, 'pressure', 'properties',
                                         *args, **kwargs)


class rotconsts(Attribute):
    """rotational constants (array[2], GHz)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('rotconsts', np.ndarray, 'rotational constants', 'atoms:coords:rotconsts',
                                         *args, **kwargs)


class scancoords(Attribute):
    """geometries of each scan step (array[3], angstroms)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('scancoords', np.ndarray, 'step geometry', 'optimization:scan',
                                         *args, **kwargs)


class scanenergies(Attribute):
    """energies of potential energy surface (list)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('scanenergies', list, 'PES energies', 'optimization:scan',
                                         *args, **kwargs)


class scannames(Attribute):
    """names of varaibles scanned (list of strings)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('scannames', list, 'variable names', 'optimization:scan',
                                         *args, **kwargs)


class scanparm(Attribute):
    """values of parameters in potential energy surface (list of tuples)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('scanparm', list, 'PES parameter values', 'optimization:scan',
                                         *args, **kwargs)


class scfenergies(Attribute):
    """molecular electronic energies after SCF (Hartree-Fock, DFT) (array[1], eV)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('scfenergies', np.ndarray, 'scf energies', 'optimization:scf',
                                         *args, **kwargs)


class scftargets(Attribute):
    """targets for convergence of the SCF (array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('scftargets', np.ndarray, 'targets', 'optimization:scf',
                                         *args, **kwargs)


class scfvalues(Attribute):
    """current values for convergence of the SCF (list of arrays[2])"""

    # Attributes that should be lists of arrays (double precision).
    _listsofarrays = ['mocoeffs', 'moenergies', 'moments', 'polarizabilities', 'scfvalues']
    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('scfvalues', list, 'values', 'optimization:scf',
                                         *args, **kwargs)


class temperature(Attribute):
    """temperature used for Thermochemistry (float, kelvin)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('temperature', float, 'temperature', 'properties',
                                         *args, **kwargs)


class time(Attribute):
    """time in molecular dynamics and other trajectories (array[1], fs)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('time', np.ndarray, 'time', 'N/A',
                                         *args, **kwargs)


class transprop(Attribute):
    """all absorption and emission spectra (dictionary {name:(etenergies, etoscs)})

    WARNING: this attribute is not standardized and is liable to change in cclib 2.0
    """

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('transprop', dict, 'electronic transitions', 'transitions',
                                         *args, **kwargs)


class vibanharms(Attribute):
    """vibrational anharmonicity constants (array[2], 1/cm)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('vibanharms', np.ndarray, 'anharmonicity constants', 'vibrations',
                                         *args, **kwargs)


class vibdisps(Attribute):
    """cartesian displacement vectors (array[3], delta angstrom)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('vibdisps', np.ndarray, 'displacement', 'vibrations',
                                         *args, **kwargs)


class vibfreqs(Attribute):
    """vibrational frequencies (array[1], 1/cm)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('vibfreqs', np.ndarray, 'frequencies', 'vibrations',
                                         *args, **kwargs)


class vibirs(Attribute):
    """IR intensities (array[1], km/mol)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('vibirs', np.ndarray, 'IR', 'vibrations:intensities',
                                         *args, **kwargs)


class vibramans(Attribute):
    """Raman intensities (array[1], A^4/Da)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('vibramans', np.ndarray, 'raman', 'vibrations:intensities',
                                         *args, **kwargs)


class vibrmasses(Attribute):
    """reduced masses of vibrations (array[1], daltons)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('vibrmasses', np.ndarray, 'reduced masses', 'vibrations',
                                         *args, **kwargs)


class vibsyms(Attribute):
    """symmetries of vibrations (list of strings)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('vibsyms', list, 'vibration symmetry', 'vibrations',
                                         *args, **kwargs)


class zpve(Attribute):
    """zero-point vibrational energy correction (float, hartree/particle)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__('zpve', float, 'zero-point correction', 'properties:energies',
                                         *args, **kwargs)


# https://stackoverflow.com/questions/1796180/how-can-i-get-a-list-of-all-classes-within-current-module-in-python#1796247
attributes = dict()
for name, obj in inspect.getmembers(sys.modules[__name__]):
    if inspect.isclass(obj):
        if name != 'Attribute':
            attributes[name] = obj
