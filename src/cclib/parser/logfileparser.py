"""
cclib is a parser for computational chemistry log files.

See http://cclib.sf.net for more information.

Copyright (C) 2006 Noel O'Boyle and Adam Tenderholt

 This program is free software; you can redistribute and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY, without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

Contributions (monetary as well as code :-) are encouraged.
"""
import logging, sys
import Numeric
import utils

class Logfile(object):
    """Abstract class for logfile objects.

    Subclasses:
        ADF, GAMESS, Gaussian, Jaguar
    
    Attributes:
        aonames -- "Ru_3p" (list)
        aooverlaps -- atomic orbital overlap matrix (array[2])
        atomnos -- atomic numbers (array)
        etenergies -- energy of electronic transitions (array[1], 1/cm)
        etoscs -- oscillator strength of electronic transition (array[1], ??)
        etrotats -- rotatory strength of electronic transitions (array[1], ??)
        etsecs -- singly-excited configurations comprising each electronic transition (??)
        etsyms -- symmetry of electronic transition (list)
        geotargets -- targets for convergence of the geometry (array[1])
        geovalues -- current values for convergence of the geometry (array[1], same units as geotargets)
        homos -- molecular orbital index of HOMO(s) (array[1])
        mocoeffs -- molecular orbital coefficients (array[3])
        moenergies -- orbital energies (array[2], eV)
        mosyms -- orbital symmetries (list[2])
        natom -- number of atoms (integer)
        nbasis -- number of basis functions (integer)
        nmo -- number of linearly-independent basis functions (integer)
        scfenergies -- the electronic energy of the molecule (array[1], a.u.)
        scftargets -- targets for convergence of the SCF (array[1])
        scfvalues -- current values for convergence of the SCF (array[2], same units as scftargets)
        vibfreqs -- vibrational frequencies (array, 1/cm)
        vibirs -- IR intensity (array, ??)
        vibramans -- Raman intensity (array, ??)
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

        # Set up the logger
        self.logger = logging.getLogger('%s %s' % (self.logname,self.filename))
        self.logger.setLevel(self.loglevel)
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter("[%(name)s %(levelname)s] %(message)s"))
        self.logger.addHandler(handler)

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
