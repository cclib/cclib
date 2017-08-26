import sys
import inspect

import numpy


class Attribute(object):

    def __init__(self, name, main_type, json_key, json_path, *args, **kwargs):

        self.name = __name__
        self.main_type = main_type
        self.json_key = json_key
        self.json_path = json_path

        self._impl = None

    def typeify(self):
        """Check the types of all attributes.

        If an attribute does not match the expected type, then attempt to
        convert; if that fails, only then raise a TypeError.
        """

        try:
            self._impl = self.main_type(self._impl)
        except ValueError:
            args = (self.name, type(self._impl), self.main_type)
            raise TypeError("attribute {} is {} instead of {} and could not be converted".format(*args))


class aonames(Attribute):
    def __init__(self):
        super().__init__('aonames', list, 'names', 'atoms:orbitals')


class aooverlaps(Attribute):
    def __init__(self):
        super().__init__('aooverlaps', numpy.ndarray, 'overlaps', 'properties:orbitals')


class atombasis(Attribute):
    def __init__(self):
        super().__init__('atombasis', list, 'indices', 'atoms:orbitals')


class atomcharges(Attribute):
    # Attributes that should be dictionaries of arrays (double precision).
    _dictsofarrays = ["atomcharges", "atomspins"]
    def __init__(self):
        super().__init__('atomcharges', dict, 'partial charges', 'properties')


class atomcoords(Attribute):
    def __init__(self):
        super().__init__('atomcoords', numpy.ndarray, 'coords', 'atoms:coords:3d')


class atommasses(Attribute):
    def __init__(self):
        super().__init__('atommasses', numpy.ndarray, 'mass', 'atoms')


class atomnos(Attribute):
    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ['atomnos', 'coreelectrons', 'homos', 'optstatus']
    def __init__(self):
        super().__init__('atomnos', numpy.ndarray, 'number', 'atoms:elements')


class atomspins(Attribute):
    # Attributes that should be dictionaries of arrays (double precision).
    _dictsofarrays = ["atomcharges", "atomspins"]
    def __init__(self):
        super().__init__('atomspins', dict, 'spins', 'atoms')


class ccenergies(Attribute):
    def __init__(self):
        super().__init__('ccenergies', numpy.ndarray, 'coupled cluster', 'properties:energy')


class charge(Attribute):
    def __init__(self):
        super().__init__('charge', int, 'charge', 'properties')


class coreelectrons(Attribute):
    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ['atomnos', 'coreelectrons', 'homos', 'optstatus']
    def __init__(self):
        super().__init__('coreelectrons', numpy.ndarray, 'core electrons', 'atoms')


class enthalpy(Attribute):
    def __init__(self):
        super().__init__('enthalpy', float, 'enthalpy', 'properties')


class entropy(Attribute):
    def __init__(self):
        super().__init__('entropy', float, 'entropy', 'properties')


class etenergies(Attribute):
    def __init__(self):
        super().__init__('etenergies', numpy.ndarray, 'electronic transitions', 'transitions')


class etoscs(Attribute):
    def __init__(self):
        super().__init__('etoscs', numpy.ndarray, 'oscillator strength', 'transitions')


class etrotats(Attribute):
    def __init__(self):
        super().__init__('etrotats', numpy.ndarray, 'rotatory strength', 'transitions')


class etsecs(Attribute):
    def __init__(self):
        super().__init__('etsecs', list, 'one excited config', 'transitions')


class etsyms(Attribute):
    def __init__(self):
        super().__init__('etsyms', list, 'symmetry', 'transitions')


class freeenergy(Attribute):
    def __init__(self):
        super().__init__('freeenergy', float, 'free energy', 'properties:energy')


class fonames(Attribute):
    def __init__(self):
        super().__init__('fonames', list, 'orbital names', 'fragments')


class fooverlaps(Attribute):
    def __init__(self):
        super().__init__('fooverlaps', numpy.ndarray, 'orbital overlap', 'fragments')


class fragnames(Attribute):
    def __init__(self):
        super().__init__('fragnames', list, 'fragment names', 'fragments')


class frags(Attribute):
    def __init__(self):
        super().__init__('frags', list, 'atom indices', 'fragments')


class gbasis(Attribute):
    def __init__(self):
        super().__init__('gbasis', list, 'basis functions', 'atoms:orbitals')


class geotargets(Attribute):
    def __init__(self):
        super().__init__('geotargets', numpy.ndarray, 'geometric targets', 'optimization')


class geovalues(Attribute):
    def __init__(self):
        super().__init__('geovalues', numpy.ndarray, 'geometric values', 'optimization')


class grads(Attribute):
    def __init__(self):
        super().__init__('grads', numpy.ndarray, 'TBD', 'N/A')


class hessian(Attribute):
    def __init__(self):
        super().__init__('hessian', numpy.ndarray, 'hessian matrix', 'vibrations')


class homos(Attribute):
    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ['atomnos', 'coreelectrons', 'homos', 'optstatus']
    def __init__(self):
        super().__init__('homos', numpy.ndarray, 'homos', 'properties:orbitals')


class metadata(Attribute):
    def __init__(self):
        super().__init__('metadata', dict, 'TBD', 'N/A')


class mocoeffs(Attribute):
    # Attributes that should be lists of arrays (double precision).
    _listsofarrays = ['mocoeffs', 'moenergies', 'moments', 'polarizabilities', 'scfvalues']
    def __init__(self):
        super().__init__('mocoeffs', list, 'coeffs', 'properties:orbitals')


class moenergies(Attribute):
    # Attributes that should be lists of arrays (double precision).
    _listsofarrays = ['mocoeffs', 'moenergies', 'moments', 'polarizabilities', 'scfvalues']
    def __init__(self):
        super().__init__('moenergies', list, 'energies', 'properties:orbitals')


class moments(Attribute):
    # Attributes that should be lists of arrays (double precision).
    _listsofarrays = ['mocoeffs', 'moenergies', 'moments', 'polarizabilities', 'scfvalues']
    def __init__(self):
        super().__init__('moments', list, 'total dipole moment', 'properties')


class mosyms(Attribute):
    def __init__(self):
        super().__init__('mosyms', list, 'molecular orbital symmetry', 'properties:orbitals')


class mpenergies(Attribute):
    def __init__(self):
        super().__init__('mpenergies', numpy.ndarray, 'moller plesset', 'properties:energy')


class mult(Attribute):
    def __init__(self):
        super().__init__('mult', int, 'multiplicity', 'properties')


class natom(Attribute):
    def __init__(self):
        super().__init__('natom', int, 'number of atoms', 'properties')


class nbasis(Attribute):
    def __init__(self):
        super().__init__('nbasis', int, 'basis number', 'properties:orbitals')


class nmo(Attribute):
    def __init__(self):
        super().__init__('nmo', int, 'MO number', 'properties:orbitals')


class nocoeffs(Attribute):
    def __init__(self):
        super().__init__('nocoeffs', numpy.ndarray, 'TBD', 'N/A')


class nooccnos(Attribute):
    def __init__(self):
        super().__init__('nooccnos', numpy.ndarray, 'TBD', 'N/A')


class optdone(Attribute):
    def __init__(self):
        super().__init__('optdone', bool, 'done', 'optimization')


class optstatus(Attribute):
    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ['atomnos', 'coreelectrons', 'homos', 'optstatus']
    def __init__(self):
        super().__init__('optstatus', numpy.ndarray, 'status', 'optimization')


class polarizabilities(Attribute):
    # Attributes that should be lists of arrays (double precision).
    _listsofarrays = ['mocoeffs', 'moenergies', 'moments', 'polarizabilities', 'scfvalues']
    def __init__(self):
        super().__init__('polarizabilities', list, 'polarizabilities', 'N/A')


class pressure(Attribute):
    def __init__(self):
        super().__init__('pressure', float, 'pressure', 'properties')


class scancoords(Attribute):
    def __init__(self):
        super().__init__('scancoords', numpy.ndarray, 'step geometry', 'optimization:scan')


class scanenergies(Attribute):
    def __init__(self):
        super().__init__('scanenergies', list, 'PES energies', 'optimization:scan')


class scannames(Attribute):
    def __init__(self):
        super().__init__('scannames', list, 'variable names', 'optimization:scan')


class scanparm(Attribute):
    def __init__(self):
        super().__init__('scanparm', list, 'PES parameter values', 'optimization:scan')


class scfenergies(Attribute):
    def __init__(self):
        super().__init__('scfenergies', numpy.ndarray, 'scf energies', 'optimization:scf')


class scftargets(Attribute):
    def __init__(self):
        super().__init__('scftargets', numpy.ndarray, 'targets', 'optimization:scf')


class scfvalues(Attribute):
    # Attributes that should be lists of arrays (double precision).
    _listsofarrays = ['mocoeffs', 'moenergies', 'moments', 'polarizabilities', 'scfvalues']
    def __init__(self):
        super().__init__('scfvalues', list, 'values', 'optimization:scf')


class temperature(Attribute):
    def __init__(self):
        super().__init__('temperature', float, 'temperature', 'properties')


class time(Attribute):
    def __init__(self):
        super().__init__('time', numpy.ndarray, 'time', 'N/A')


class vibanharms(Attribute):
    def __init__(self):
        super().__init__('vibanharms', numpy.ndarray, 'anharmonicity constants', 'vibrations')


class vibdisps(Attribute):
    def __init__(self):
        super().__init__('vibdisps', numpy.ndarray, 'displacement', 'vibrations')


class vibfreqs(Attribute):
    def __init__(self):
        super().__init__('vibfreqs', numpy.ndarray, 'frequencies', 'vibrations')


class vibirs(Attribute):
    def __init__(self):
        super().__init__('vibirs', numpy.ndarray, 'IR', 'vibrations:intensities')


class vibramans(Attribute):
    def __init__(self):
        super().__init__('vibramans', numpy.ndarray, 'raman', 'vibrations:intensities')


class vibsyms(Attribute):
    def __init__(self):
        super().__init__('vibsyms', list, 'vibration symmetry', 'vibrations')


# https://stackoverflow.com/questions/1796180/how-can-i-get-a-list-of-all-classes-within-current-module-in-python#1796247
_attributes = dict()
for name, obj in inspect.getmembers(sys.modules[__name__]):
    if inspect.isclass(obj):
        if name != 'Attribute':
            _attributes[name] = obj

_attrlist = sorted(_attributes.keys())
