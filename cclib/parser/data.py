# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Classes and tools for storing and handling parsed data"""

import logging
from collections import namedtuple
from typing import Any, Dict, List, Mapping, Optional

from cclib.method import Electrons, orbitals
from cclib.properties.property import Property, _properties

import numpy


class ccData:
    """Stores data extracted by cclib parsers

    Description of cclib attributes:
        aonames -- atomic orbital names (list of strings)
        aooverlaps -- atomic orbital overlap matrix (array[2])
        atombasis -- indices of atomic orbitals on each atom (list of lists)
        atomcharges -- atomic partial charges (dict of arrays[1])
        atomcoords -- atom coordinates (array[3], angstroms)
        atommasses -- atom masses (array[1], daltons)
        atomnos -- atomic numbers (array[1])
        atomspins -- atomic spin densities (dict of arrays[1])
        ccenergies -- molecular energies with Coupled-Cluster corrections (array[2], eV)
        charge -- net charge of the system (integer)
        coreelectrons -- number of core electrons in atom pseudopotentials (array[1])
        dispersionenergies -- dispersion energy corrections (array[1], eV)
        enthalpy -- sum of electronic and thermal enthalpies (float, hartree/particle)
        entropy -- entropy (float, hartree/(particle*kelvin))
        etenergies -- energies of electronic transitions (array[1], 1/cm)
        etoscs -- oscillator strengths of electronic transitions (array[1])
        etdips -- electric transition dipoles of electronic transitions (array[2], ebohr)
        etveldips -- velocity-gauge electric transition dipoles of electronic transitions (array[2], ebohr)
        etmagdips -- magnetic transition dipoles of electronic transitions (array[2], ebohr)
        etrotats -- rotatory strengths of electronic transitions (array[1], ??)
        etsecs -- singly-excited configurations for electronic transitions (list of lists)
        etsyms -- symmetries of electronic transitions (list of string)
        freeenergy -- sum of electronic and thermal free energies (float, hartree/particle)
        fonames -- fragment orbital names (list of strings)
        fooverlaps -- fragment orbital overlap matrix (array[2])
        fragnames -- names of fragments (list of strings)
        frags -- indices of atoms in a fragment (list of lists)
        gbasis -- coefficients and exponents of Gaussian basis functions (PyQuante format)
        geotargets -- targets for convergence of geometry optimization (array[1])
        geovalues -- current values for convergence of geometry optmization (array[1])
        grads -- current values of forces (gradients) in geometry optimization (array[3])
        hessian -- elements of the force constant matrix (array[1])
        homos -- molecular orbital indices of HOMO(s) (array[1])
        metadata -- various metadata about the package and computation (dict)
        mocoeffs -- molecular orbital coefficients (list of arrays[2])
        moenergies -- molecular orbital energies (list of arrays[1], eV)
        moments -- molecular multipole moments (list of arrays[], a.u.)
        mosyms -- orbital symmetries (list of lists)
        mpenergies -- molecular electronic energies with Møller-Plesset corrections (array[2], eV)
        mult -- multiplicity of the system (integer)
        natom -- number of atoms (integer)
        nbasis -- number of basis functions (integer)
        nmo -- number of molecular orbitals (integer)
        nmrtensors -- Nuclear magnetic resonance chemical shielding tensors (dict of dicts of array[2])
        nmrcouplingtensors -- Nuclear magnetic resonance spin-spin coupling tensors (dict of dicts of array[2])
        nocoeffs -- natural orbital coefficients (array[2])
        nooccnos -- natural orbital occupation numbers (array[1])
        nsocoeffs -- natural spin orbital coefficients (list of array[2])
        nsooccnos -- natural spin orbital occupation numbers (list of array[1])
        optdone -- flags whether an optimization has converged (Boolean)
        optstatus -- optimization status for each set of atomic coordinates (array[1])
        polarizabilities -- (dipole) polarizabilities, static or dynamic (list of arrays[2])
        pressure -- pressure used for Thermochemistry (float, atm)
        rotconsts -- rotational constants (array[2], GHz)
        scancoords -- geometries of each scan step (array[3], angstroms)
        scanenergies -- energies of potential energy surface (list)
        scannames -- names of variables scanned (list of strings)
        scanparm -- values of parameters in potential energy surface (list of tuples)
        scfenergies -- molecular electronic energies after SCF (Hartree-Fock, DFT) (array[1], eV)
        scftargets -- targets for convergence of the SCF (array[2])
        scfvalues -- current values for convergence of the SCF (list of arrays[2])
        temperature -- temperature used for Thermochemistry (float, kelvin)
        time -- time in molecular dynamics and other trajectories (array[1], fs)
        transprop -- all absorption and emission spectra (dictionary {name:(etenergies, etoscs)})
            WARNING: this attribute is not standardized and is liable to change in cclib 2.0
        vibanharms -- vibrational anharmonicity constants (array[2], 1/cm)
        vibdisps -- cartesian displacement vectors (array[3], delta angstrom)
        vibfreqs -- vibrational frequencies (array[1], 1/cm)
        vibfconsts -- force constants of vibrations (array[1], mDyne/angstrom)
        vibirs -- IR intensities (array[1], km/mol)
        vibramans -- Raman activities (array[1], A^4/Da)
        vibrmasses -- reduced masses of vibrations (array[1], daltons)
        vibsyms -- symmetries of vibrations (list of strings)
        zpve -- zero-point vibrational energy correction (float, hartree/particle)
    (1) The term 'array' refers to a numpy array
    (2) The number of dimensions of an array is given in square brackets
    (3) Python indexes arrays/lists starting at zero, so if homos==[10], then
            the 11th molecular orbital is the HOMO
    """

    # The expected types for all supported attributes.
    # The json_key is the key name used for attributes in the CJSON/JSON format
    # 'TBD' - To Be Decided are the key names of attributes which haven't been included in the cjson format

    # The name of all attributes can be generated from the dictionary above.
    _attrlist = ["scfenergies"]  # sorted(_properties.keys())

    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ["atomnos", "coreelectrons", "homos", "optstatus"]

    # Propertys that should be lists of arrays (double precision).
    _listsofarrays = ["mocoeffs", "moenergies", "moments", "polarizabilities", "scfvalues"]

    # Propertys that should be dictionaries of arrays (double precision).
    _dictsofarrays = ["atomcharges", "atomspins"]

    # Propertys that should be dictionaries of dictionaries.
    _dictsofdicts = ["populations"]

    # Possible statuses for optimization steps.
    # OPT_UNKNOWN is the default and means optimization is in progress.
    # OPT_NEW is set for every new optimization (e.g. PES, IRCs, etc.)
    # OPT_DONE is set for the last step of an optimisation that converged.
    # OPT_UNCONVERGED is set for every unconverged step (e.g. should be mutually exclusive with OPT_DONE)
    # bit value notation allows coding for multiple states: OPT_NEW and OPT_UNCONVERGED or OPT_NEW and OPT_DONE.
    OPT_UNKNOWN = 0b000
    OPT_NEW = 0b001
    OPT_UNCONVERGED = 0b010
    OPT_DONE = 0b100

    def __init__(self, attributes: Mapping[str, Any] = {}) -> None:
        """Initialize the cclibData object.

        Normally called in the parse() method of a Logfile subclass.

        Inputs:
            attributes - optional dictionary of attributes to load as data
        """

        self._parsed_properties = dict()

        if attributes:
            self.setattributes(attributes)

    def listify(self) -> None:
        """Converts all attributes that are arrays or lists/dicts of arrays to lists."""

        attrlist = [k for k in self._attrlist if hasattr(self, k)]
        for k in attrlist:
            v = _properties[k].type
            if v == numpy.ndarray:
                setattr(self, k, getattr(self, k).tolist())
            elif v == list and k in self._listsofarrays:
                setattr(self, k, [x.tolist() for x in getattr(self, k)])
            elif v == dict and k in self._dictsofarrays:
                items = getattr(self, k).items()
                pairs = [(key, val.tolist()) for key, val in items]
                setattr(self, k, dict(pairs))
            elif v == dict and k in self._dictsofdicts:
                items = getattr(self, k).items()
                pairs = [
                    (key, {subkey: subval.tolist()})
                    for key, val in items
                    for subkey, subval in val.items()
                ]
                setattr(self, k, dict(pairs))

    def arrayify(self) -> None:
        """Converts appropriate attributes to arrays or lists/dicts of arrays."""

        attrlist = [k for k in self._attrlist if hasattr(self, k)]
        for k in attrlist:
            v = _properties[k].type
            precision = "d"
            if k in self._intarrays:
                precision = "i"
            if v == numpy.ndarray:
                setattr(self, k, numpy.array(getattr(self, k), precision))
            elif v == list and k in self._listsofarrays:
                setattr(self, k, [numpy.array(x, precision) for x in getattr(self, k)])
            elif v == dict and k in self._dictsofarrays:
                items = getattr(self, k).items()
                pairs = [(key, numpy.array(val, precision)) for key, val in items]
                setattr(self, k, dict(pairs))
            elif v == dict and k in self._dictsofdicts:
                items = getattr(self, k).items()
                pairs = [
                    (
                        key,
                        {
                            subkey: numpy.array(subval, precision)
                            if isinstance(subval, (int, float))
                            else numpy.array(subval)
                        },
                    )
                    for key, val in items
                    for subkey, subval in val.items()
                ]
                setattr(self, k, dict(pairs))

    def getattributes(self, tolists: bool = False) -> Dict[str, Any]:
        """Returns a dictionary of existing data attributes.

        Inputs:
            tolists - flag to convert attributes to lists where applicable
        """

        if tolists:
            self.listify()
        attributes = {}
        for attr in self._attrlist:
            if hasattr(self, attr):
                attributes[attr] = getattr(self, attr)
        if tolists:
            self.arrayify()
        return attributes

    def setattributes(self, attributes: Mapping[str, Any]) -> List[str]:
        """Sets data attributes given in a dictionary.

        Inputs:
            attributes - dictionary of attributes to set
        Outputs:
            invalid - list of attributes names that were not set, which
                      means they are not specified in self._attrlist
        """

        if not isinstance(attributes, dict):
            raise TypeError("attributes must be in a dictionary")

        valid = [a for a in attributes if a in self._attrlist]
        invalid = [a for a in attributes if a not in self._attrlist]

        for attr in valid:
            setattr(self, attr, attributes[attr])

        self.arrayify()
        self.typecheck()

        return invalid

    def typecheck(self) -> None:
        """Check the types of all attributes.

        If an attribute does not match the expected type, then attempt to
        convert; if that fails, only then raise a TypeError.
        """

        self.arrayify()
        for attr in [a for a in self._attrlist if hasattr(self, a)]:
            # attr.typecheck()

            val = getattr(self, attr)
            if isinstance(val, _properties[attr].type):
                continue

            try:
                val = _properties[attr].type(val)
            except ValueError:
                args = (attr, type(val), _properties[attr].type)
                raise TypeError(
                    f"attribute {args[0]} is {args[1]} instead of {args[2]} and could not be converted"
                )

    def check_values(self, logger=logging) -> None:
        """Perform custom checks on the values of attributes."""
        if hasattr(self, "etenergies") and any(e < 0 for e in self.etenergies):
            negative_values = [e for e in self.etenergies if e < 0]
            msg = f"At least one excitation energy is negative. \nNegative values: {negative_values}\nFull etenergies: {self.etenergies}"
            logger.error(msg)

    def write(
        self, filename: Optional[str] = None, indices: Optional = None, *args, **kwargs
    ) -> str:
        """Write parsed attributes to a file.

        Possible extensions:
          .cjson or .json -  output a chemical JSON file
          .cml - output a chemical markup language (CML) file
          .xyz - output a Cartesian XYZ file of the last coordinates available
        """

        from cclib.io import ccwrite

        outputstr = ccwrite(self, outputdest=filename, indices=indices, *args, **kwargs)
        return outputstr

    def writejson(self, filename: Optional[str] = None, indices=None):
        """Write parsed attributes to a JSON file."""
        return self.write(filename=filename, indices=indices, outputtype="cjson")

    def writecml(self, filename: Optional[str] = None, indices=None):
        """Write parsed attributes to a CML file."""
        return self.write(filename=filename, indices=indices, outputtype="cml")

    def writexyz(self, filename: Optional[str] = None, indices=None):
        """Write parsed attributes to an XML file."""
        return self.write(filename=filename, indices=indices, outputtype="xyz")

    @property
    def converged_geometries(self) -> numpy.ndarray:
        """
        Return all converged geometries.

        An array containing only the converged geometries, e.g.:
            - For PES or IRCs, return all geometries for which optstatus matches OPT_DONE
            - The converged geometry for simple optimisations
            - The input geometry for single points
        """
        if hasattr(self, "optstatus"):
            converged_indexes = [x for x, y in enumerate(self.optstatus) if y & self.OPT_DONE > 0]
            return self.atomcoords[converged_indexes]
        else:
            return self.atomcoords

    @property
    def new_geometries(self) -> numpy.ndarray:
        """
        Return all starting geometries.

        An array containing only the starting geometries, e.g.:
            - For PES or IRCs, return all geometries for which optstatus matches OPT_NEW
            - The input geometry for simple optimisations or single points
        """
        if hasattr(self, "optstatus"):
            new_indexes = [x for x, y in enumerate(self.optstatus) if y & self.OPT_NEW > 0]
            return self.atomcoords[new_indexes]
        else:
            return self.atomcoords

    @property
    def unknown_geometries(self) -> numpy.ndarray:
        """
        Return all OPT_UNKNOWN geometries.

        An array containing only the starting geometries, e.g.:
            - For PES or IRCs, return all geometries for which optstatus matches OPT_UNKNOWN
            - The input geometry for simple optimisations or single points
        """
        if hasattr(self, "optstatus"):
            unknown_indexes = [x for x, y in enumerate(self.optstatus) if y == self.OPT_UNKNOWN]
            return self.atomcoords[unknown_indexes]
        else:
            return self.atomcoords

    @property
    def unconverged_geometries(self) -> numpy.ndarray:
        """
        Return all unconverged geometries.

        An array containing only the starting geometries, e.g.:
            - For PES or IRCs, return all geometries for which optstatus matches OPT_UNCONVERGED
            - The input geometry for simple optimisations or single points
        """
        if hasattr(self, "optstatus"):
            unconverged_indexes = [
                x for x, y in enumerate(self.optstatus) if y & self.OPT_UNCONVERGED > 0
            ]
            return self.atomcoords[unconverged_indexes]
        else:
            return self.atomcoords

    @property
    def nelectrons(self) -> int:
        return Electrons(self).count()

    @property
    def closed_shell(self) -> bool:
        return orbitals.Orbitals(self).closed_shell()

    def __setattr__(self, name: str, value: Any) -> None:
        if name in _properties:
            self._parsed_properties[name] = value
        else:
            super().__setattr__(name, value)

    def __getattr__(self, name: str) -> Any:
        # If we couldn't find an attribute directly on the class, which, for
        # an Property, should actually be a property, then it's not
        # implemented as a property yet and is in our special attribute
        # container.
        try:
            return self._parsed_properties[name]
        except KeyError:
            pass
            # raise PropertyError

    @property
    def aonames(self):
        try:
            return self._parsed_properties["aonames"]
        except:
            pass
        # except KeyError:
        #    raise PropertyError

    # @aonames.setter
    # def aonames(self, val):
    #     setattr(self, "aonames", val)
