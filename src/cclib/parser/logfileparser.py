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
import math,sys,logging,copy,re,os,time # How many of these are necessary?
import Numeric

def convertor(value,fromunits,tounits):
    """Convert from one set of units to another.

    >>> convertor(8,"eV","cm-1")
    64000
    """
    _convertor = {"eV_to_cm-1": lambda x: x*8065.6,
                  "nm_to_cm-1": lambda x: 1e7/x,
                  "cm-1_to_nm": lambda x: 1e7/x}

    return _convertor["%s_to_%s" % (fromunits,tounits)] (value)

class PeriodicTable(object):
    """Allows conversion between element name and atomic no.

    >>> t = PeriodicTable()
    >>> t.element[6]
    'C'
    >>> t.number['C']
    6
    """
    def __init__(self):
        self.element = [None,"H","He","Li","Be","B","C","N","O","F","Ne"]
        self.number = {}
        for i in range(1,len(self.element)):
            self.number[self.element[i]] = i

class Logfile(object):
    """Abstract class for logfile objects.

    Subclasses:
       G03
    
    Attributes:
        aonames -- list of "Ru_3p", etc.
        aooverlaps -- the atomic orbital overlap matrix
        atomnos -- atomic numbers
        etenergies -- energy of electronic transitions (i.e. UV-Vis, CD)
        etoscs -- oscillator strength of electronic transition
        etrotats -- rotatory strength(?) of electronic transitions (for CD)
        etsecs -- singly-excited configurations comprising each electronic transition
        etsyms -- symmetry of electronic transition
        geotargets -- targets for convergence of the geometry
        geovalues -- current values for convergence of the geometry
        homos -- location of HOMO(s)
        mocoeffs -- molecular orbital coefficients
        moenergies -- orbital energies
        mosyms -- orbital symmetries
        natom -- number of atoms
        nbasis -- number of basis functions
        nindep -- number of linearly-independent basis functions
        progress -- class (or None) for handling progress info
        scftargets -- targets for convergence of the SCF
        scfvalues -- current values for convergence of the SCF
        vibfreqs -- vibrational frequencies
        vibirs -- IR intensity
        vibramans -- Raman intensity
        vibsyms -- symmetry of vibrations    
    """
    def __init__(self,filename,progress=None):
        self.filename = filename
        self.progress = progress

    def float(self,number):
        """Convert a string to a float avoiding the problem with Ds."""
        number = number.replace("D","E")
        return float(number)

if __name__=="__main__":
    import doctest,parser
    doctest.testmod(parser,verbose=False)
