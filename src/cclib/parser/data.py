# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2007-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

import numpy


class ccData(object):
    """Class for objects containing data from cclib parsers and methods.

    Description of cclib attributes:
        aonames -- atomic orbital names (list)
        aooverlaps -- atomic orbital overlap matrix (array[2])
        atombasis -- indices of atomic orbitals on each atom (list of lists)
        atomcharges -- atomic partial charges (dict of arrays[1])
        atomcoords -- atom coordinates (array[3], angstroms)
        atommasses -- atom masses (array[1], daltons)
        atomnos -- atomic numbers (array[1])
        atomspins -- atomic spin densities (dict of arrays[1])
        charge -- net charge of the system (integer)
        ccenergies -- molecular energies with Coupled-Cluster corrections (array[2], eV)
        coreelectrons -- number of core electrons in atom pseudopotentials (array[1])
        etenergies -- energies of electronic transitions (array[1], 1/cm)
        etoscs -- oscillator strengths of electronic transitions (array[1])
        etrotats -- rotatory strengths of electronic transitions (array[1], ??)
        etsecs -- singly-excited configurations for electronic transitions (list of lists)
        etsyms -- symmetries of electronic transitions (list)
        fonames -- fragment orbital names (list)
        fooverlaps -- fragment orbital overlap matrix (array[2])
        fragnames -- names of fragments (list)
        frags -- indices of atoms in a fragment (list of lists)
        gbasis -- coefficients and exponents of Gaussian basis functions (PyQuante format)
        geotargets -- targets for convergence of geometry optimization (array[1])
        geovalues -- current values for convergence of geometry optmization (array[1])
        grads -- current values of forces (gradients) in geometry optimization (array[3])
        hessian -- elements of the force constant matrix (array[1])
        homos -- molecular orbital indices of HOMO(s) (array[1])
        mocoeffs -- molecular orbital coefficients (list of arrays[2])
        moenergies -- molecular orbital energies (list of arrays[1], eV)
        mosyms -- orbital symmetries (list of lists)
        mpenergies -- molecular electronic energies with Moller-Plesset corrections (array[2], eV)
        mult -- multiplicity of the system (integer)
        natom -- number of atoms (integer)
        nbasis -- number of basis functions (integer)
        nmo -- number of molecular orbitals (integer)
        nocoeffs -- natural orbital coefficients (array[2])
        scfenergies -- molecular electronic energies after SCF (Hartree-Fock, DFT) (array[1], eV)
        scftargets -- targets for convergence of the SCF (array[2])
        scfvalues -- current values for convergence of the SCF (list of arrays[2])
        vibanharms -- vibrational anharmonicity constants (array[2], 1/cm)
        vibdisps -- cartesian displacement vectors (array[3], delta angstrom)
        vibfreqs -- vibrational frequencies (array[1], 1/cm)
        vibirs -- IR intensities (array[1], km/mol)
        vibramans -- Raman intensities (array[1], A^4/Da)
        vibsyms -- symmetries of vibrations (list)
        scannames -- Names of varaibles scanned (list)
        scanenergies -- energies of potential energy surface (list)
        scanparm -- values of parameters in potential energy surface (list of tuples)
        scancoords -- Geometries of each scan step (array[3], angstroms)
        enthaply -- Sum of electronic and thermal Enthalpie (float hartree/particle)
        freeenergy -- Sum of electronic and thermal Free Energies (float hartree/particle)
        temperature -- Tempature used for Thermochemistry (float kelvin)
        entropy -- Entropy (float hartree/particle)
        optdone -- Stores if an optimisation job has completed (boolean)
    (1) The term 'array' refers to a numpy array
    (2) The number of dimensions of an array is given in square brackets
    (3) Python indexes arrays/lists starting at zero, so if homos==[10], then
            the 11th molecular orbital is the HOMO
    """

    # Names of all supported attributes.
    _attrlist = [
        'aonames', 'aooverlaps', 'atombasis', 'atomcharges', 'atomcoords',
        'atommasses', 'atomnos', 'atomspins',
        'ccenergies', 'charge', 'coreelectrons',
        'etenergies', 'etoscs', 'etrotats', 'etsecs', 'etsyms',
        'fonames', 'fooverlaps', 'fragnames', 'frags',
        'gbasis', 'geotargets', 'geovalues', 'grads',
        'hessian', 'homos',
        'mocoeffs', 'moenergies', 'mosyms', 'mpenergies', 'mult',
        'natom', 'nbasis', 'nmo', 'nocoeffs',
        'scfenergies', 'scftargets', 'scfvalues',
        'vibanharms', 'vibdisps', 'vibfreqs', 'vibirs',
        'vibramans', 'vibsyms', 'scannames', 'scanenergies', 'scanparm',
        'scancoords', 'enthaply', 'freeenergy', 'temperature', 'entropy', 
        'optdone'
    ]

    # The expected types for all supported attributes in dectionary,
    # and their names which can be extracted as the keys.
    _attrtypes = {
        "aonames":        list,
        "aooverlaps":     numpy.ndarray,
        "atombasis":      list,
        "atomcharges":    dict,
        "atomcoords":     numpy.ndarray,
        "atommasses":     numpy.ndarray,
        "atomnos":        numpy.ndarray,
        "atomspins":      dict,
        "ccenergies":     numpy.ndarray,
        "charge":         int,
        "coreelectrons":  numpy.ndarray,
        "etenergies":     numpy.ndarray,
        "etoscs":         numpy.ndarray,
        "etrotats":       numpy.ndarray,
        "etsecs":         list,
        "etsyms":         list,
        "fonames":        list,
        "fooverlaps":     numpy.ndarray,
        "fragnames":      list,
        "frags":          list,
        'gbasis':         list,
        "geotargets":     numpy.ndarray,
        "geovalues":      numpy.ndarray,
        "grads":          numpy.ndarray,
        "hessian":        numpy.ndarray,
        "homos":          numpy.ndarray,
        "mocoeffs":       list,
        "moenergies":     list,
        "mosyms":         list,
        "mpenergies":     numpy.ndarray,
        "mult":           int,
        "natom":          int,
        "nbasis":         int,
        "nmo":            int,
        "nocoeffs":       numpy.ndarray,
        "scfenergies":    numpy.ndarray,
        "scftargets":     numpy.ndarray,
        "scfvalues":      list,
        "vibanharms":     numpy.ndarray,
        "vibdisps":       numpy.ndarray,
        "vibfreqs":       numpy.ndarray,
        "vibirs":         numpy.ndarray,
        "vibramans":      numpy.ndarray,
        "vibsyms":        list,
        "scannames":      list,
        "scanenergies":   list,
        "scanparm":       list,
        "scancoords":     numpy.ndarray,
        "enthaply":       float,
        "freeenergy":     float,
        "temperature":    float,
        "entropy":        float,
        "optdone":        bool
    }

    # Arrays are double precision by default, but these will be integer arrays.
    _intarrays = ['atomnos', 'coreelectrons', 'homos']

    # Attributes that should be lists of arrays (double precision).
    _listsofarrays = ['mocoeffs', 'moenergies', 'scfvalues']
    
    # Attributes that should be dictionaries of arrays (double precision).
    _dictsofarrays = ["atomcharges", "atomspins"]

    def __init__(self, attributes=None):
        """Initialize the cclibData object.
        
        Normally called in the parse() method of a Logfile subclass.
        
        Inputs:
            attributes - dictionary of attributes to load
        """

        if attributes:
            self.setattributes(attributes)
        
    def listify(self):
        """Converts all attributes that are arrays or lists/dicts of arrays to lists."""
        
        attrlist = [k for k in self._attrlist if hasattr(self, k)]
        for k in attrlist:
            v = self._attrtypes[k]
            if v == numpy.ndarray:
                setattr(self, k, getattr(self, k).tolist())
            elif v == list and k in self._listsofarrays:
                setattr(self, k, [x.tolist() for x in getattr(self, k)])
            elif v == dict and k in self._dictsofarrays:
                items = getattr(self, k).iteritems()
                pairs = [(key, val.tolist()) for key, val in items]
                setattr(self, k, dict(pairs))
    
    def arrayify(self):
        """Converts appropriate attributes to arrays or lists/dicts of arrays."""
        
        attrlist = [k for k in self._attrlist if hasattr(self, k)]
        for k in attrlist:
            v = self._attrtypes[k]
            precision = 'd'
            if k in self._intarrays:
                precision = 'i'
            if v == numpy.ndarray:
                setattr(self, k, numpy.array(getattr(self, k), precision))
            elif v == list and k in self._listsofarrays:
                setattr(self, k, [numpy.array(x, precision) for x in getattr(self, k)])
            elif v == dict and k in self._dictsofarrays:
                items = getattr(self, k).items()
                pairs = [(key, numpy.array(val, precision)) for key, val in items]
                setattr(self, k, dict(pairs))

    def getattributes(self, tolists=False):
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

    def setattributes(self, attributes):
        """Sets data attributes given in a dictionary.
        
        Inputs:
            attributes - dictionary of attributes to set
        Outputs:
            invalid - list of attributes names that were not set, which
                      means they are not specified in self._attrlist
        """
    
        if type(attributes) is not dict:
            raise TypeError("attributes must be in a dictionary")
    
        valid = [a for a in attributes if a in self._attrlist]
        invalid = [a for a in attributes if a not in self._attrlist]
    
        for attr in valid:
            setattr(self, attr, attributes[attr])
        self.arrayify()
        return invalid
