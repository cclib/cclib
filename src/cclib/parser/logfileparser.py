"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision$"

import logging, random, sys
import Numeric
import utils

class Logfile(object):
    """Abstract class for logfile objects.

    Subclasses:
        ADF, GAMESS, GAMESSUK, Gaussian, Jaguar
    
    Attributes:
        aonames -- atomic orbital names (list)
        aooverlaps -- atomic orbital overlap matrix (array[2])
        atomcoords -- atom coordinates (array[3], angstroms)
        atomnos -- atomic numbers (array[1])
        coreelectrons -- number of core electrons in an atom's pseudopotential (array[1])
        etenergies -- energy of electronic transitions (array[1], 1/cm)
        etoscs -- oscillator strength of electronic transition (array[1])
        etrotats -- rotatory strength of electronic transitions (array[1], ??)
        etsecs -- singly-excited configurations for electronic transitions (list of lists)
        etsyms -- symmetry of electronic transition (list)
        fonames -- fragment orbital names (list)
        fooverlaps -- fragment orbital overlap matrix (array[2])
        fragnames -- names of fragments (list)
        frags -- indices of atoms in a fragment (list of lists)
        gbasis -- coefficients and exponents of Gaussian basis functions (PyQuante format)
        geotargets -- targets for convergence of the geometry (array[1])
        geovalues -- current values for convergence of the geometry (array[1])
        homos -- molecular orbital index of HOMO(s) (array[1])
        mocoeffs -- molecular orbital coefficients (list of arrays[2])
        moenergies -- orbital energies (list of arrays[1], eV)
        mosyms -- orbital symmetries (list of lists)
        mpenergies -- molecular electronic energy after Moller-Plesset correction (array[2], eV)
        natom -- number of atoms (integer)
        nbasis -- number of basis functions (integer)
        nmo -- number of molecular orbitals (integer)
        scfenergies -- the electronic energy of the molecule (array[1], eV)
        scftargets -- targets for convergence of the SCF (array[2])
        scfvalues -- current values for convergence of the SCF (list of arrays[2])
        vibdisps -- Cartesian displacement vectors (array[3], delta angstrom)
        vibfreqs -- vibrational frequencies (array[1], 1/cm)
        vibirs -- IR intensity (array[1], km/mol)
        vibramans -- Raman intensity (array[1], A^4/Da)
        vibsyms -- symmetry of vibrations (list)
    (1) The term 'array' currently refers to a Numeric array
    (2) The number of dimensions of an array is given in square brackets
    (3) Python indexes arrays/lists starting at zero. So if homos==[10], then
        the 11th molecular orbital is the HOMO
    """

    
    def __init__(self,filename,progress=None,
                 loglevel=logging.INFO,logname="Log"):
        """Initialise the Logfile object.

        Typically called by subclasses in their own __init__ methods.
        """
        self.filename = filename
        self.progress = progress
        self.parsed = False
        self.loglevel = loglevel
        self.logname  = logname
        self.table = utils.PeriodicTable()

        self.attrlist = ['aonames', 'aooverlaps', 'atomcoords', 'atomnos',
                         'ccenergies', 'coreelectrons',
                         'etenergies', 'etoscs', 'etrotats', 'etsecs', 'etsyms',
                         'fonames', 'fooverlaps',
                         'geotargets', 'geovalues', 'homos', 'mocoeffs',
                         'moenergies', 'mosyms', 'mpenergies', 'natom', 'nbasis', 'nmo',
                         'scfenergies', 'scftargets', 'scfvalues',
                         'vibfreqs', 'vibirs', 'vibramans', 'vibsyms']

        self._tofloatarray = ['atomcoords', 'ccenergies',
                              'etenergies', 'etoscs',
                              'geotargets', 'geovalues', 'mpenergies',
                              'scfenergies', 'scftargets',
                              'vibdisps', 'vibfreqs', 'vibirs', 'vibramans']
        self._tointarray = ['atomnos', 'coreelectrons', 'homos']
        self._tolistofarrays = ['moenergies', 'scfvalues']

        # Set up the logger
        self.logger = logging.getLogger('%s %s' % (self.logname,self.filename))
        self.logger.setLevel(self.loglevel)
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter("[%(name)s %(levelname)s] %(message)s"))
        self.logger.addHandler(handler)

    def parse(self, fupdate=0.05, cupdate=0.02):
        """Parse the logfile, using the assumed extract method of the child."""

        # Open the file object
        inputfile = utils.openlogfile(self.filename)

        # Intialize self.progress
        if self.progress:
            inputfile.seek(0,2)
            nstep = inputfile.tell()
            inputfile.seek(0)
            self.progress.initialize(nstep)
            self.progress.step = 0

        try:
            # This method does the actual parsing of text,
            #  and should be defined by a subclass.
            self.extract(inputfile, fupdate=fupdate, cupdate=cupdate)
        except AttributeError:
            self.logger.info("Method parse() was called from generic LogFile class.")
            return

        # Close file object
        inputfile.close()

        # Make sure selected attributes are arrays
        for attr in self._tofloatarray:
            if hasattr(self, attr):
                if type(getattr(self, attr)) is not Numeric.arraytype:
                    setattr(self, attr, Numeric.array(getattr(self, attr), 'd'))
        for attr in self._tointarray:
            if hasattr(self, attr):
                if type(getattr(self, attr)) is not Numeric.arraytype:
                    setattr(self, attr, Numeric.array(getattr(self, attr), 'i'))


        # Make sure selected attrbutes are lists of arrays
        for attr in self._tolistofarrays:
            if hasattr(self, attr):
                if not Numeric.alltrue([type(x) is Numeric.arraytype for x in getattr(self, attr)]):
                    setattr(self, attr, [Numeric.array(x, 'd') for x in getattr(self, attr)])

        # Set nmo if not set already - to nbasis
        if not hasattr(self, "nmo"):
            self.nmo = self.nbasis

        # Creating deafult coreelectrons array
        if not hasattr(self, "coreelectrons"):
            self.logger.info("Creating attribute coreelectrons[]")
            self.coreelectrons = Numeric.zeros(self.natom, "i")

        # Update self.progress as done
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
        for attr in self.attrlist:
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
