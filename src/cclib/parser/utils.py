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
import adfparser
import gamessparser
import gaussianparser
import jaguarparser

def guesstype(filename):
    """Guess the identity of a particular log file and return an instance of that.
    
    Notes: uses the first 20 lines of a log file to guess its identity

    Returns: one of ADF, GAMESS, Gaussian, Jaguar, or None (if it cannot figure it out).
    """
    filetype = None
    inputfile = open(filename,"r")
    for i,line in enumerate(inputfile):
        if i>20: break # not elegant, but less buggy than a while loop testing for EOF
        if line.find("Amsterdam Density Functional")>=0:
            filetype = adfparser.ADF
        elif line.find("GAMESS")>=0:
            filetype = gamessparser.GAMESS
        elif line.find("Gaussian")>=0:
            filetype = gaussianparser.Gaussian
        elif line.find("Jaguar")>=0:
            filetype = jaguarparser.Jaguar
    inputfile.close() # Need to close before creating an instance
    if filetype:
        filetype = apply(filetype,[filename]) # Create an instance of the chosen class
    return filetype

def convertor(value,fromunits,tounits):
    """Convert from one set of units to another.

    >>> print "%.1f" % convertor(8,"eV","cm-1")
    64524.8
    """
    _convertor = {"eV_to_cm-1": lambda x: x*8065.6,
                  "hartree_to_eV": lambda x: x*27.2114,
                  "nm_to_cm-1": lambda x: 1e7/x,
                  "cm-1_to_nm": lambda x: 1e7/x,
                  "au_to_Ang": lambda x: x*0.529177}

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
        self.element = [None,"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe", \
                        "Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo"]
        self.number = {}
        for i in range(1,len(self.element)):
            self.number[self.element[i]] = i

if __name__=="__main__":
    import doctest,utils
    doctest.testmod(utils,verbose=False)
