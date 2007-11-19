"""
cclib (http://cclib.sf.net) is (c) 2007, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""


import cPickle as pickle
import os
import sys

import numpy


# These modules are supported for read/write to JSON format by calling load_json.
json_support = [ ('simplejson','dumps', 'loads'),
                 ('cjson', 'encode', 'decode'),
                 ('json', 'write', 'read') ]


class ccData(object):
    """Class for objects containing data from cclib parsers and methods.

    Description of cclib attributes:
        aonames -- atomic orbital names (list)
        aooverlaps -- atomic orbital overlap matrix (array[2])
        atombasis -- indices of atomic orbitals on each atom (list of lists)
        atomcoords -- atom coordinates (array[3], angstroms)
        atomnos -- atomic numbers (array[1])
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
        vibdisps -- cartesian displacement vectors (array[3], delta angstrom)
        vibfreqs -- vibrational frequencies (array[1], 1/cm)
        vibirs -- IR intensities (array[1], km/mol)
        vibramans -- Raman intensities (array[1], A^4/Da)
        vibsyms -- symmetries of vibrations (list)
    (1) The term 'array' refers to a numpy array
    (2) The number of dimensions of an array is given in square brackets
    (3) Python indexes arrays/lists starting at zero, so if homos==[10], then
            the 11th molecular orbital is the HOMO
    """

    def __init__(self, attributes=None):
        """Initialize the cclibData object.
        
        Normally called in the parse() method of a Logfile subclass.
        
        Inputs:
            attributes - dictionary of attributes to load
        """

        # Names of all supported attributes.
        self._attrlist = ['aonames', 'aooverlaps', 'atombasis',
                          'atomcoords', 'atomnos',
                          'ccenergies', 'charge', 'coreelectrons',
                          'etenergies', 'etoscs', 'etrotats', 'etsecs', 'etsyms',
                          'fonames', 'fooverlaps', 'fragnames', 'frags',
                          'gbasis', 'geotargets', 'geovalues',
                          'hessian', 'homos',
                          'mocoeffs', 'moenergies', 'mosyms', 'mpenergies', 'mult',
                          'natom', 'nbasis', 'nmo', 'nocoeffs',
                          'scfenergies', 'scftargets', 'scfvalues',
                          'vibdisps', 'vibfreqs', 'vibirs', 'vibramans', 'vibsyms']

        # The expected types for all supported attributes.
        self._attrtypes = { "aonames":        list,
                            "aooverlaps":     numpy.ndarray,
                            "atombasis":      list,
                            "atomcoords":     numpy.ndarray,
                            "atomnos":        numpy.ndarray,
                            "charge":         int,
                            "coreelectrons":  numpy.ndarray,
                            "etenergies":     numpy.ndarray,
                            "etoscs":         numpy.ndarray,
                            "etrotats":       numpy.ndarray,
                            "etsecs":         list,
                            "etsyms":         list,
                            'gbasis':         list,
                            "geotargets":     numpy.ndarray,
                            "geovalues":      numpy.ndarray,
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
                            "vibdisps":       numpy.ndarray,
                            "vibfreqs":       numpy.ndarray,
                            "vibirs":         numpy.ndarray,
                            "vibramans":      numpy.ndarray,
                            "vibsyms":        list,
                          }

        # Arrays are double precision by default, but these will be integer arrays.
        self._intarrays = ['atomnos', 'coreelectrons', 'homos']

        # Attributes that should be lists of arrays (double precision).
        self._listsofarrays = ['mocoeffs', 'moenergies', 'scfvalues']
        
        if attributes:
            self.setattributes(attributes)
        
    def listify(self):
        """Converts all attributes that are arrays or lists of arrays to lists."""
        
        for k, v in self._attrtypes.iteritems():
            if hasattr(self, k):
                if v == numpy.ndarray:
                    setattr(self, k, getattr(self, k).tolist())
                elif v == list and k in self._listsofarrays:
                    setattr(self, k, [x.tolist() for x in getattr(self, k)])
    
    def arrayify(self):
        """Converts appropriate attributes to arrays or lists of arrays."""
        
        for k, v in self._attrtypes.iteritems():
            if hasattr(self, k):
                precision = 'd'
                if k in self._intarrays:
                    precision = 'i'
                if v == numpy.ndarray:
                    setattr(self, k, numpy.array(getattr(self, k), precision))
                elif v == list and k in self._listsofarrays:
                    setattr(self, k, [numpy.array(x, precision)
                                      for x in getattr(self, k)])

    def getattributes(self):
        """Returns a dictionary of existing data attributes."""
    
        attributes = {}
        for attr in self._attrlist:
            if hasattr(self, attr):
                attributes[attr] = getattr(self,attr)
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
            raise TypeError, "attributes must be in a dictionary"
    
        valid = [a for a in attributes if a in self._attrlist]
        invalid = [a for a in attributes if a not in self._attrlist]
    
        for attr in valid:
            setattr(self, attr, attributes[attr])
        self.arrayify()
        return invalid

    def writepickle(self, outfile=sys.stdout, protocol=0, listify=False, binary=False):
        """Returns pickled represenntation of data attributes in a dictionary.
        
        Inputs:
            outfile - output file name or file object, sys.stdout by default
            protocol - protocol (int) passed to pickle
            listify - flag (bool) to convert all attributes to lists
            binary - flag to write in binary mode (only if outfile is a file name)
        """
    
        if listify:
            self.listify()
            
        if type(outfile) is str:
            stream = open(outfile, "w"+"b"*binary)
            print >>stream, pickle.dumps(self.getattributes(), protocol=protocol)
            stream.close()
        elif type(outfile) is file:
            print >>outfile, pickle.dumps(self.getattributes(), protocol=protocol)
        else:
            raise TypeError, "outfile must be a file name or file object"
        
        if listify:
            self.arrayify()

    def writejson(self, outfile=sys.stdout, module=None, binary=False):
        """Returns JSON representation of data attributes in a dictionary.
        
        Inputs:
            outfile - output file name or file object, stdout by default
            module - optional name of JSON support module
            binary - flag to write in binary mode (only if outfile is a file name)
        """

        support = load_json(module)
        if not support:
            raise ImportError, "No JSON support module found"
        module, encoder, decoder = support

        # JSON does not understand numpy arrays, so listify them first.
        self.listify()

        if type(outfile) is str:
            stream = open(outfile, "w"+"b"*binary)
            print >>stream, encoder(self.getattributes())
            stream.close()
        elif type(outfile) is file:
            print >>outfile, encoder(self.getattributes())
        else:
            raise TypeError, "outfile must be a file name or file object"

        # Change all neccesary attributes back to arrays.
        self.arrayify()


def load_json(module=None):
    """Load a JSON support module."""
    
    imported = None
    for support in json_support:
        if not module or module == support[0]:
            try:
                exec "import %s" %support[0]
                module = eval(support[0])
                encode = getattr(module, support[1])
                decode = getattr(module, support[2])
                return module, encode, decode
            except ImportError:
                pass

def readpickle(indata=sys.stdin, binary=False):
    """Read pickle data and return a corresponding ccData object.
    
    Inputs:
        indata - raw pickle string, file name or file object, sys.stdin by default
        binary - flag to read in binary mode (only if outfile is a file name)
    """
    
    if type(indata) is str:
        if os.path.exists(indata):
            stream = open(indata,'r'+'b'*binary)
            attributes = pickle.loads(stream.read())
            stream.close()
        else:
            attributes = pickle.loads(indata)
    elif type(indata) is file:
        attributes = pickle.loads(indata.read())
    else:
        raise TypeError, "indata must be a string or file object"

    return ccData(attributes=attributes)

def readjson(indata=sys.stdin, module=None, binary=False):
    """Read JSON data and return a corresponding ccData object.
    
    Inputs:
        indata - raw JSON string, file name or file object, sys.stdin by default
        module - optional name of JSON support module
        binary - flag to read in binary mode (only if outfile is a file name)
    """
    
    support = load_json(module)
    if not support:
        raise ImportError, "No JSON support module found"
    module, encoder, decoder = support

    if type(indata) is str:
        if os.path.exists(indata):
            stream = open(indata,'r'+'b'*binary)
            attributes = decoder(stream.read())
            stream.close()
        else:
            attributes = decoder(indata)
    elif type(indata) is file:
        attributes = decoder(indata.read())
    else:
        raise TypeError, "indata must be a string or file object"

    return ccData(attributes=attributes)
