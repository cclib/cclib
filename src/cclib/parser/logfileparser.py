"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision$"

import sys
import logging
import inspect
import random

# If NumPy is not installed, try to import Numeric instead.
try:
    import numpy
except ImportError:
    import Numeric as numpy
    numpy.ndarray = numpy.arraytype

import utils


class Logfile(object):
    """Abstract class for logfile objects.

    Subclasses defined by cclib:
        ADF, GAMESS, GAMESSUK, Gaussian, Jaguar
    
    Attributes:
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
    (3) Python indexes arrays/lists starting at zero. So if homos==[10], then
        the 11th molecular orbital is the HOMO
    """

    def __init__(self, filename, progress=None, fupdate=0.05, cupdate=0.002, 
                                 loglevel=logging.INFO, logname="Log"):
        """Initialise the Logfile object.

        Typically called by subclasses in their own __init__ methods.
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
                          'natom', 'nbasis', 'nmo',
                          'scfenergies', 'scftargets', 'scfvalues',
                          'vibdisps', 'vibfreqs', 'vibirs', 'vibramans', 'vibsyms']

        # The expected types for all attributes.
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
                            "scfenergies":    numpy.ndarray,
                            "scftargets":     numpy.ndarray,
                            "scfvalues":      list,
                            "vibdisps":       numpy.ndarray,
                            "vibfreqs":       numpy.ndarray,
                            "vibirs":         numpy.ndarray,
                            "vibramans":      numpy.ndarray,
                            "vibsyms":        list,
                          }

        # Arrays are double precision by default, these will be integer arrays.
        self._tointarray = ['atomnos', 'coreelectrons', 'homos']

        # Attributes that should be lists of arrays.
        self._tolistofarrays = ['mocoeffs', 'moenergies', 'scfvalues']

        self.filename = filename

        # Progress indicator.
        self.progress = progress
        self.fupdate = fupdate
        self.cupdate = cupdate

        # Setup the logger.
        # Note that calling logging.getLogger() with one name always returns the same instance.
        # Presently in cclib, all parser instances of the same class use the same logger.
        # This means that care needs to be taken not to duplicate handlers.
        self.loglevel = loglevel
        self.logname  = logname
        self.logger = logging.getLogger('%s %s' % (self.logname,self.filename))
        self.logger.setLevel(self.loglevel)
        if len(self.logger.handlers) == 0:
                handler = logging.StreamHandler(sys.stdout)
                handler.setFormatter(logging.Formatter("[%(name)s %(levelname)s] %(message)s"))
                self.logger.addHandler(handler)

        # Periodic table of elements.
        self.table = utils.PeriodicTable()

        self.parsed = False

    def __setattr__(self, name, value):

        if hasattr(self, "logger"):
            # Only call logger.info if attribute is new and in list.
            if not hasattr(self, name) and name in self._attrlist:
                if type(value) in [numpy.ndarray, list]:
                    self.logger.info("Creating attribute %s[]" %name)
                else:
                    self.logger.info("Creating attribute %s: %s" %(name, str(value)))

        object.__setattr__(self, name, value)

    def parse(self, fupdate=None, cupdate=None):
        """Parse the logfile, using the assumed extract method of the child."""

        # Check that the sub-class has an extract attribute,
        #  that is callable with the proper number of arguemnts.
        if not hasattr(self, "extract"):
            raise AttributeError, "Class %s has no extract() method." %self.__class__.__name__
            return -1
        if not callable(self.extract):
            raise AttributeError, "Method %s._extract not callable." %self.__class__.__name__
            return -1
        if len(inspect.getargspec(self.extract)[0]) != 3:
            raise AttributeError, "Method %s._extract takes wrong number of arguments." %self.__class__.__name__
            return -1

        # Create list of attributes to keep after parsing.
        _nodelete = list(set(self.__dict__.keys() + self._attrlist))

        # Open the file object.
        inputfile = utils.openlogfile(self.filename)

        # Intialize self.progress.
        if self.progress:
            inputfile.seek(0,2)
            nstep = inputfile.tell()
            inputfile.seek(0)
            self.progress.initialize(nstep)
            self.progress.step = 0
            if fupdate:
                self.fupdate = fupdate
            if cupdate:
                self.cupdate = cupdate

        # Maybe the sub-class has something to do before parsing.
        if hasattr(self, "before_parsing"):
            self.before_parsing()

        # Loop over lines in the file object.
        for line in inputfile:

            self.updateprogress(inputfile, "Unsupported information", cupdate)

            # This call should check if the line begins a section of extracted data.
            # If it does, it parses some lines and sets the relevant attributes.
            self.extract(inputfile, line)

        # Close file object
        inputfile.close()

        # Make sure all attributes have correct type (including arrays).
        for attr in self._attrlist:
            if hasattr(self, attr):
                atype = self._attrtypes.get(attr, None)
                if atype and type(getattr(self, attr)) is not atype:
                    if atype is numpy.ndarray:
                        try:
                            precision = 'd'
                            if attr in self._tointarray:
                                precision = 'i'
                            setattr(self, attr, numpy.array(getattr(self, attr), precision))
                        except TypeError:
                            self.logger.info("Attribute %s cannot be converted to an array" %(attr))
                    else:
                        try:
                            setattr(self, attr, atype(getattr(self, attr)))
                        except ValueError:
                            self.logger.info("Attribute %s cannot be converted to type '%s'" %(attr, atype))

        # Make sure selected attrbutes are lists of arrays.
        for attr in self._tolistofarrays:
            if hasattr(self, attr):
                if not numpy.alltrue([type(x) is numpy.ndarray for x in getattr(self, attr)]):
                    try:
                        setattr(self, attr, [numpy.array(x, 'd') for x in getattr(self, attr)])
                    except ValueError:
                        self.logger.info("Elements of attribute %s cannot be converted to arrays." %attr)

        # If atomcoords were not parsed, but some input coordinates were ("inputcoords").
        # This is originally from the Gaussian parser, a regression fix.
        if not hasattr(self, "atomcoords") and hasattr(self, "inputcoords"):
            self.atomcoords = numpy.array(self.inputcoords, 'd')

        # Set nmo if not set already - to nbasis.
        if not hasattr(self, "nmo"):
            self.nmo = self.nbasis

        # Creating deafult coreelectrons array.
        if not hasattr(self, "coreelectrons"):
            self.coreelectrons = numpy.zeros(self.natom, "i")

        # Delete temporary attributes (set during parsing and not in attrlist).
        for attr in self.__dict__.keys():
            if not attr in _nodelete:
                self.__delattr__(attr)

        # Update self.progress as done.
        if self.progress:
            self.progress.update(nstep, "Done")

        self.parsed = True

    def updateprogress(self, inputfile, msg, xupdate=0.05):
        """Update progress."""

        if self.progress and random.random() < xupdate:
            newstep = inputfile.tell()
            if newstep != self.progress.step:
                self.progress.update(newstep, msg)
                self.progress.step = newstep

    def clean(self):
        """Delete all of the parsed attributes."""
        for attr in self._attrlist:
            if hasattr(self, attr):
                delattr(self, attr)
        self.parsed = False

    def normalisesym(self,symlabel):
        """Standardise the symmetry labels between parsers.

        This method should be overwritten by individual parsers, and should
        contain appropriate doctests. If is not overwritten, this is detected
        as an error by unit tests.
        """
        return "ERROR: This should be overwritten by this subclass"

    def float(self,number):
        """Convert a string to a float avoiding the problem with Ds.

        >>> t = Logfile("dummyfile")
        >>> t.float("123.2323E+02")
        12323.23
        >>> t.float("123.2323D+02")
        12323.23
        """
        number = number.replace("D","E")
        return float(number)

if __name__=="__main__":
    import doctest,logfileparser
    doctest.testmod(logfileparser,verbose=False)
