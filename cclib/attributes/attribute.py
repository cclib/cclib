# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

# import sys
# import inspect
from abc import ABC
from collections import namedtuple

import numpy as np


class _Attribute(ABC):
    def __init__(
        self, name: str, main_type, json_key: str, json_path: str, *args, **kwargs
    ) -> None:
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
                raise AttributeError(f"{type(self)} doesn't contain any data")
            else:
                return self._impl

    def __str__(self) -> str:
        return str(self.__class__.__name__)

    def typecheck(self) -> None:
        """Check the types of all attributes.

        If an attribute does not match the expected type, then attempt to
        convert; if that fails, only then raise a TypeError.
        """

        if isinstance(self._impl, self.main_type):
            return

        try:
            if self.main_type == np.ndarray:
                self._impl = np.asarray(self._impl)
            else:
                self._impl = self.main_type(self._impl)
        except ValueError:
            args = (self.name, type(self._impl), self.main_type)
            raise TypeError(
                "attribute {} is {} instead of {} and could not be converted".format(*args)
            )


class Aonames(_Attribute):
    """atomic orbital names (list of strings)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "aonames", list, "names", "atoms:orbitals", *args, **kwargs
        )


class Aooverlaps(_Attribute):
    """atomic orbital overlap matrix (array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "aooverlaps", np.ndarray, "overlaps", "properties:orbitals", *args, **kwargs
        )


class Atombasis(_Attribute):
    """indices of atomic orbitals on each atom (list of lists)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "atombasis", list, "indices", "atoms:orbitals", *args, **kwargs
        )


class Atomcharges(_Attribute):
    """atomic partial charges (dict of arrays[1])"""

    # Attributes that should be dictionaries of arrays (double precision).
    _dictsofarrays = ["atomcharges", "atomspins"]

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "atomcharges", dict, "partial charges", "properties", *args, **kwargs
        )


class Atomcoords(_Attribute):
    """atom coordinates (array[3], angstroms)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "atomcoords", np.ndarray, "coords", "atoms:coords:3d", *args, **kwargs
        )


class Atommasses(_Attribute):
    """atom masses (array[1], daltons)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__("atommasses", np.ndarray, "mass", "atoms", *args, **kwargs)


class Atomnos(_Attribute):
    """atomic numbers (array[1])"""

    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ["atomnos", "coreelectrons", "homos", "optstatus"]

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "atomnos", np.ndarray, "number", "atoms:elements", *args, **kwargs
        )


class Atomspins(_Attribute):
    """atomic spin densities (dict of arrays[1])"""

    # Attributes that should be dictionaries of arrays (double precision).
    _dictsofarrays = ["atomcharges", "atomspins"]

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__("atomspins", dict, "spins", "atoms", *args, **kwargs)


class Ccenergies(_Attribute):
    """molecular energies with Coupled-Cluster corrections (array[2], hartree)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "ccenergies", np.ndarray, "coupled cluster", "properties:energy", *args, **kwargs
        )


class Charge(_Attribute):
    """net charge of the system (integer)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__("charge", int, "charge", "properties", *args, **kwargs)


class Coreelectrons(_Attribute):
    """number of core electrons in atom pseudopotentials (array[1])"""

    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ["atomnos", "coreelectrons", "homos", "optstatus"]

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "coreelectrons", np.ndarray, "core electrons", "atoms", *args, **kwargs
        )


class Dispersionenergies(_Attribute):
    """a molecular dispersion energy corrections (array[1], hartree)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "dispersionenergies",
            np.ndarray,
            "dispersion correction",
            "properties:energy",
            *args,
            **kwargs,
        )


class Enthalpy(_Attribute):
    """sum of electronic and thermal enthalpies (float, hartree/particle)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "enthalpy", float, "enthalpy", "properties", *args, **kwargs
        )


class Entropy(_Attribute):
    """entropy (float, hartree/particle)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__("entropy", float, "entropy", "properties", *args, **kwargs)


class Etenergies(_Attribute):
    """energies of electronic transitions (array[1], hartree)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "etenergies", np.ndarray, "electronic transitions", "transitions", *args, **kwargs
        )


class Etoscs(_Attribute):
    """oscillator strengths of electronic transitions (array[1])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "etoscs", np.ndarray, "oscillator strength", "transitions", *args, **kwargs
        )


class Etdips(_Attribute):
    """electric transition dipoles of electronic transitions (array[2], ebohr)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "etdips", np.ndarray, "electric transition dipoles", "transitions", *args, **kwargs
        )


class Etveldips(_Attribute):
    """velocity-gauge electric transition dipoles of electronic transitions (array[2], ebohr)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "etveldips",
            np.ndarray,
            "velocity-gauge electric transition dipoles",
            "transitions",
            *args,
            **kwargs,
        )


class Etmagdips(_Attribute):
    """magnetic transition dipoles of electronic transitions (array[2], ebohr)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "etmagdips", np.ndarray, "magnetic transition dipoles", "transitions", *args, **kwargs
        )


class Etrotats(_Attribute):
    """rotatory strengths of electronic transitions (array[1], ??)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "etrotats", np.ndarray, "rotatory strength", "transitions", *args, **kwargs
        )


class Etsecs(_Attribute):
    """singly-excited configurations for electronic transitions (list of lists)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "etsecs", list, "one excited config", "transitions", *args, **kwargs
        )


class Etsyms(_Attribute):
    """symmetries of electronic transitions (list of string)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__("etsyms", list, "symmetry", "transitions", *args, **kwargs)


class Freeenergy(_Attribute):
    """sum of electronic and thermal free energies (float, hartree/particle)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "freeenergy", float, "free energy", "properties:energy", *args, **kwargs
        )


class Fonames(_Attribute):
    """fragment orbital names (list of strings)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "fonames", list, "orbital names", "fragments", *args, **kwargs
        )


class Fooverlaps(_Attribute):
    """fragment orbital overlap matrix (array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "fooverlaps", np.ndarray, "orbital overlap", "fragments", *args, **kwargs
        )


class Fragnames(_Attribute):
    """names of fragments (list of strings)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "fragnames", list, "fragment names", "fragments", *args, **kwargs
        )


class Frags(_Attribute):
    """indices of atoms in a fragment (list of lists)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "frags", list, "atom indices", "fragments", *args, **kwargs
        )


class Gbasis(_Attribute):
    """coefficients and exponents of Gaussian basis functions (PyQuante format)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "gbasis", list, "basis functions", "atoms:orbitals", *args, **kwargs
        )


class Geotargets(_Attribute):
    """targets for convergence of geometry optimization (array[1])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "geotargets", np.ndarray, "geometric targets", "optimization", *args, **kwargs
        )


class Geovalues(_Attribute):
    """current values for convergence of geometry optmization (array[1])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "geovalues", np.ndarray, "geometric values", "optimization", *args, **kwargs
        )


class Grads(_Attribute):
    """current values of forces (gradients) in geometry optimization (array[3])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__("grads", np.ndarray, "TBD", "N/A", *args, **kwargs)


class Hessian(_Attribute):
    """elements of the force constant matrix (array[1])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "hessian", np.ndarray, "hessian matrix", "vibrations", *args, **kwargs
        )


class Homos(_Attribute):
    """molecular orbital indices of HOMO(s) (array[1])"""

    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ["atomnos", "coreelectrons", "homos", "optstatus"]

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "homos", np.ndarray, "homos", "properties:orbitals", *args, **kwargs
        )


class Metadata(_Attribute):
    """various metadata about the package and computation (dict)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__("metadata", dict, "TBD", "N/A", *args, **kwargs)


class Mocoeffs(_Attribute):
    """molecular orbital coefficients (list of arrays[2])"""

    # Attributes that should be lists of arrays (double precision).
    _listsofarrays = ["mocoeffs", "moenergies", "moments", "polarizabilities", "scfvalues"]

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "mocoeffs", list, "coeffs", "properties:orbitals", *args, **kwargs
        )


class Moenergies(_Attribute):
    """molecular orbital energies (list of arrays[1], hartree)"""

    # Attributes that should be lists of arrays (double precision).
    _listsofarrays = ["mocoeffs", "moenergies", "moments", "polarizabilities", "scfvalues"]

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "moenergies", list, "energies", "properties:orbitals", *args, **kwargs
        )


class Moments(_Attribute):
    """molecular multipole moments (list of arrays[], a.u.)"""

    # Attributes that should be lists of arrays (double precision).
    _listsofarrays = ["mocoeffs", "moenergies", "moments", "polarizabilities", "scfvalues"]

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "moments", list, "total dipole moment", "properties", *args, **kwargs
        )


class Mosyms(_Attribute):
    """orbital symmetries (list of lists)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "mosyms", list, "molecular orbital symmetry", "properties:orbitals", *args, **kwargs
        )


class Mpenergies(_Attribute):
    """molecular electronic energies with Møller-Plesset corrections (array[2], hartree)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "mpenergies", np.ndarray, "moller plesset", "properties:energy", *args, **kwargs
        )


class Mult(_Attribute):
    """multiplicity of the system (integer)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__("mult", int, "multiplicity", "properties", *args, **kwargs)


class Natom(_Attribute):
    """number of atoms (integer)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "natom", int, "number of atoms", "properties", *args, **kwargs
        )


class Nbasis(_Attribute):
    """number of basis functions (integer)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "nbasis", int, "basis number", "properties:orbitals", *args, **kwargs
        )


class Nmo(_Attribute):
    """number of molecular orbitals (integer)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "nmo", int, "MO number", "properties:orbitals", *args, **kwargs
        )


class Nmrtensors(_Attribute):
    """Nuclear magnetic resonance chemical shielding tensors (dict of dicts of array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "nmr", int, "NMR chemical shielding tensors", "properties:nmr", *args, **kwargs
        )


class Nmrcouplingtensors(_Attribute):
    """Nuclear magnetic resonance spin-spin coupling tensors (dict of dicts of array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "nmr", int, "NMR spin-spin coupling tensors", "properties:nmr", *args, **kwargs
        )


class Nocoeffs(_Attribute):
    """natural orbital coefficients (array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__("nocoeffs", np.ndarray, "TBD", "N/A", *args, **kwargs)


class Nooccnos(_Attribute):
    """natural orbital occupation numbers (array[1])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__("nooccnos", np.ndarray, "TBD", "N/A", *args, **kwargs)


class Nsocoeffs(_Attribute):
    """natural spin orbital coefficients (list of array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__("nsocoeffs", list, "TBD", "N/A", *args, **kwargs)


class Nsooccnos(_Attribute):
    """natural spin orbital occupation numbers (list of array[1])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__("nsooccnos", list, "TBD", "N/A", *args, **kwargs)


class Optdone(_Attribute):
    """flags whether an optimization has converged (Boolean)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__("optdone", bool, "done", "optimization", *args, **kwargs)


class Optstatus(_Attribute):
    """optimization status for each set of atomic coordinates (array[1])"""

    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ["atomnos", "coreelectrons", "homos", "optstatus"]

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "optstatus", np.ndarray, "status", "optimization", *args, **kwargs
        )


class Polarizabilities(_Attribute):
    """(dipole) polarizabilities, static or dynamic (list of arrays[2])"""

    # Attributes that should be lists of arrays (double precision).
    _listsofarrays = ["mocoeffs", "moenergies", "moments", "polarizabilities", "scfvalues"]

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "polarizabilities", list, "polarizabilities", "N/A", *args, **kwargs
        )


class Pressure(_Attribute):
    """temperature used for Thermochemistry (float, atm)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "pressure", float, "pressure", "properties", *args, **kwargs
        )


class Rotconsts(_Attribute):
    """rotational constants (array[2], GHz)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "rotconsts",
            np.ndarray,
            "rotational constants",
            "atoms:coords:rotconsts",
            *args,
            **kwargs,
        )


class Scancoords(_Attribute):
    """geometries of each scan step (array[3], angstroms)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "scancoords", np.ndarray, "step geometry", "optimization:scan", *args, **kwargs
        )


class Scanenergies(_Attribute):
    """energies of potential energy surface (list)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "scanenergies", list, "PES energies", "optimization:scan", *args, **kwargs
        )


class Scannames(_Attribute):
    """names of varaibles scanned (list of strings)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "scannames", list, "variable names", "optimization:scan", *args, **kwargs
        )


class Scanparm(_Attribute):
    """values of parameters in potential energy surface (list of tuples)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "scanparm", list, "PES parameter values", "optimization:scan", *args, **kwargs
        )


class Scfenergies(_Attribute):
    """molecular electronic energies after SCF (Hartree-Fock, DFT) (array[1], eV)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "scfenergies", np.ndarray, "scf energies", "optimization:scf", *args, **kwargs
        )


class Scftargets(_Attribute):
    """targets for convergence of the SCF (array[2])"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "scftargets", np.ndarray, "targets", "optimization:scf", *args, **kwargs
        )


class Scfvalues(_Attribute):
    """current values for convergence of the SCF (list of arrays[2])"""

    # Attributes that should be lists of arrays (double precision).
    _listsofarrays = ["mocoeffs", "moenergies", "moments", "polarizabilities", "scfvalues"]

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "scfvalues", list, "values", "optimization:scf", *args, **kwargs
        )


class Temperature(_Attribute):
    """temperature used for Thermochemistry (float, kelvin)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "temperature", float, "temperature", "properties", *args, **kwargs
        )


class Time(_Attribute):
    """time in molecular dynamics and other trajectories (array[1], fs)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__("time", np.ndarray, "time", "N/A", *args, **kwargs)


class Transprop(_Attribute):
    """all absorption and emission spectra (dictionary {name:(etenergies, etoscs)})

    WARNING: this attribute is not standardized and is liable to change in cclib 2.0
    """

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "transprop", dict, "electronic transitions", "transitions", *args, **kwargs
        )


class Vibanharms(_Attribute):
    """vibrational anharmonicity constants (array[2], 1/cm)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "vibanharms", np.ndarray, "anharmonicity constants", "vibrations", *args, **kwargs
        )


class Vibdisps(_Attribute):
    """cartesian displacement vectors (array[3], delta angstrom)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "vibdisps", np.ndarray, "displacement", "vibrations", *args, **kwargs
        )


class Vibfreqs(_Attribute):
    """vibrational frequencies (array[1], 1/cm)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "vibfreqs", np.ndarray, "frequencies", "vibrations", *args, **kwargs
        )


class Vibfconsts(_Attribute):
    """force constants of vibrations (array[1], mDyne/angstrom)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "vibfreqs", np.ndarray, "force constants", "vibrations", *args, **kwargs
        )


class Vibirs(_Attribute):
    """IR intensities (array[1], km/mol)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "vibirs", np.ndarray, "IR", "vibrations:intensities", *args, **kwargs
        )


class Vibramans(_Attribute):
    """Raman intensities (array[1], A^4/Da)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "vibramans", np.ndarray, "raman", "vibrations:intensities", *args, **kwargs
        )


class Vibrmasses(_Attribute):
    """reduced masses of vibrations (array[1], daltons)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "vibrmasses", np.ndarray, "reduced masses", "vibrations", *args, **kwargs
        )


class Vibsyms(_Attribute):
    """symmetries of vibrations (list of strings)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "vibsyms", list, "vibration symmetry", "vibrations", *args, **kwargs
        )


class Zpve(_Attribute):
    """zero-point vibrational energy correction (float, hartree/particle)"""

    def __init__(self, *args, **kwargs):
        super(type(self), self).__init__(
            "zpve", float, "zero-point correction", "properties:energies", *args, **kwargs
        )


# https://stackoverflow.com/questions/1796180/how-can-i-get-a-list-of-all-classes-within-current-module-in-python#1796247
# attributes = dict()
# for name, obj in inspect.getmembers(sys.modules[__name__]):
#     if inspect.isclass(obj):
#         if name != '_Attribute':
#             attributes[name] = obj

Attribute = namedtuple("Attribute", ["type", "json_key", "attribute_path"])

# The expected types for all supported attributes.
# The json_key is the key name used for attributes in the CJSON/JSON format
# 'TBD' - To Be Decided are the key names of attributes which haven't been included in the cjson format
_attributes = {
    "aonames": Attribute(list, "names", "atoms:orbitals"),
    "aooverlaps": Attribute(np.ndarray, "overlaps", "properties:orbitals"),
    "atombasis": Attribute(list, "indices", "atoms:orbitals"),
    "atomcharges": Attribute(dict, "partial charges", "properties"),
    "atomcoords": Attribute(np.ndarray, "coords", "atoms:coords:3d"),
    "atommasses": Attribute(np.ndarray, "mass", "atoms"),
    "atomnos": Attribute(np.ndarray, "number", "atoms:elements"),
    "atomspins": Attribute(dict, "spins", "atoms"),
    "ccenergies": Attribute(np.ndarray, "coupled cluster", "properties:energy"),
    "charge": Attribute(int, "charge", "properties"),
    "coreelectrons": Attribute(np.ndarray, "core electrons", "atoms"),
    "dispersionenergies": Attribute(np.ndarray, "dispersion correction", "properties:energy"),
    "enthalpy": Attribute(float, "enthalpy", "properties"),
    "entropy": Attribute(float, "entropy", "properties"),
    "etenergies": Attribute(np.ndarray, "electronic transitions", "transitions"),
    "etoscs": Attribute(np.ndarray, "oscillator strength", "transitions"),
    "etdips": Attribute(np.ndarray, "electic transition dipoles", "transitions"),
    "etveldips": Attribute(np.ndarray, "velocity-gauge electric transition dipoles", "transitions"),
    "etmagdips": Attribute(np.ndarray, "magnetic transition dipoles", "transitions"),
    "etrotats": Attribute(np.ndarray, "rotatory strength", "transitions"),
    "etsecs": Attribute(list, "one excited config", "transitions"),
    "etsyms": Attribute(list, "symmetry", "transitions"),
    "freeenergy": Attribute(float, "free energy", "properties:energy"),
    "fonames": Attribute(list, "orbital names", "fragments"),
    "fooverlaps": Attribute(np.ndarray, "orbital overlap", "fragments"),
    "fragnames": Attribute(list, "fragment names", "fragments"),
    "frags": Attribute(list, "atom indices", "fragments"),
    "gbasis": Attribute(list, "basis functions", "atoms:orbitals"),
    "geotargets": Attribute(np.ndarray, "geometric targets", "optimization"),
    "geovalues": Attribute(np.ndarray, "geometric values", "optimization"),
    "grads": Attribute(np.ndarray, "TBD", "N/A"),
    "hessian": Attribute(np.ndarray, "hessian matrix", "vibrations"),
    "homos": Attribute(np.ndarray, "homos", "properties:orbitals"),
    "metadata": Attribute(dict, "TBD", "N/A"),
    "mocoeffs": Attribute(list, "coeffs", "properties:orbitals"),
    "moenergies": Attribute(list, "energies", "properties:orbitals"),
    "moments": Attribute(list, "total dipole moment", "properties"),
    "mosyms": Attribute(list, "molecular orbital symmetry", "properties:orbitals"),
    "mpenergies": Attribute(np.ndarray, "moller plesset", "properties:energy"),
    "mult": Attribute(int, "multiplicity", "properties"),
    "natom": Attribute(int, "number of atoms", "properties"),
    "nbasis": Attribute(int, "basis number", "properties:orbitals"),
    "nmo": Attribute(int, "MO number", "properties:orbitals"),
    "nmrtensors": Attribute(dict, "NMR chemical shielding tensors", "properties:nmr"),
    "nmrcouplingtensors": Attribute(dict, "NMR spin-spin coupling tensors", "properties:nmr"),
    "nocoeffs": Attribute(np.ndarray, "TBD", "N/A"),
    "nooccnos": Attribute(np.ndarray, "TBD", "N/A"),
    "nsocoeffs": Attribute(list, "TBD", "N/A"),
    "nsooccnos": Attribute(list, "TBD", "N/A"),
    "optdone": Attribute(list, "done", "optimization"),
    "optstatus": Attribute(np.ndarray, "status", "optimization"),
    "polarizabilities": Attribute(list, "polarizabilities", "N/A"),
    "pressure": Attribute(float, "pressure", "properties"),
    "rotconsts": Attribute(np.ndarray, "rotational constants", "atoms:coords:rotconsts"),
    "scancoords": Attribute(np.ndarray, "step geometry", "optimization:scan"),
    "scanenergies": Attribute(list, "PES energies", "optimization:scan"),
    "scannames": Attribute(list, "variable names", "optimization:scan"),
    "scanparm": Attribute(list, "PES parameter values", "optimization:scan"),
    "scfenergies": Attribute(np.ndarray, "scf energies", "optimization:scf"),
    "scftargets": Attribute(np.ndarray, "targets", "optimization:scf"),
    "scfvalues": Attribute(list, "values", "optimization:scf"),
    "temperature": Attribute(float, "temperature", "properties"),
    "time": Attribute(np.ndarray, "time", "N/A"),
    "transprop": Attribute(dict, "electronic transitions", "transitions"),
    "vibanharms": Attribute(np.ndarray, "anharmonicity constants", "vibrations"),
    "vibdisps": Attribute(np.ndarray, "displacement", "vibrations"),
    "vibfreqs": Attribute(np.ndarray, "frequencies", "vibrations"),
    "vibfconsts": Attribute(np.ndarray, "force constants", "vibrations"),
    "vibirs": Attribute(np.ndarray, "IR", "vibrations:intensities"),
    "vibramans": Attribute(np.ndarray, "raman", "vibrations:intensities"),
    "vibrmasses": Attribute(np.ndarray, "reduced masses", "vibrations"),
    "vibsyms": Attribute(list, "vibration symmetry", "vibrations"),
    "zpve": Attribute(float, "zero-point correction", "properties:energies"),
}
