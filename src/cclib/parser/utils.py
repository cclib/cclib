# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

def convertor(value, fromunits, tounits):
    """Convert from one set of units to another.
    Source: NIST 2010 CODATA, GAMESS-US docs and others

    >>> print "%.1f" % convertor(8, "eV", "cm-1")
    64524.8
    """
    _convertor = {
        
        # Taken taken from NIST 2010 CODATA reference
        # Available at: http://physics.nist.gov/cuu/Constants/index.html
        "Angstrom_to_bohr": lambda x: x*1.889716,
        "bohr_to_Angstrom": lambda x: x*0.529177,
    
        "cm-1_to_eV":       lambda x: x/8065.54429,
        "cm-1_to_hartree":  lambda x: x/219474.6313708,
        "cm-1_to_kcal":     lambda x: x/349.7550112,
        "cm-1_to_kJmol-1":  lambda x: x/83.5934722814,
        "cm-1_to_nm":       lambda x: 1e7/x,
        
        "eV_to_cm-1":       lambda x: x*8065.54429,
        "eV_to_hartree":    lambda x: x/27.21138505,
        "eV_to_kcal":       lambda x: x*23.060548867,
        "eV_to_kJmol-1":    lambda x: x*96.4853364596,

        "hartree_to_cm-1":  lambda x: x*219474.6313708,
        "hartree_to_eV":    lambda x: x*27.21138505,
        "hartree_to_kcal":  lambda x: x*627.50947414,
        "hartree_to_kJmol-1":lambda x: x*2625.4996398,

        "kcal_to_cm-1":     lambda x: x*349.7550112,
        "kcal_to_eV":       lambda x: x/23.060548867,
        "kcal_to_hartree":  lambda x: x/627.50947414,
        "kcal_to_kJmol-1":  lambda x: x*4.184,

        "kJmol-1_to_cm-1":  lambda x: x*83.5934722814,
        "kJmol-1_to_eV":    lambda x: x/96.4853364596,
        "kJmol-1_to_hartree": lambda x: x/2625.49963978,
        "kJmol-1_to_kcal":  lambda x: x/4.184,
        "nm_to_cm-1":       lambda x: 1e7/x,

        # Taken from GAMESS docs, "Further information",
        # "Molecular Properties and Conversion Factors"
        "Debye^2/amu-Angstrom^2_to_km/mol": lambda x: x*42.255}

    return _convertor["%s_to_%s" % (fromunits, tounits)] (value)


class PeriodicTable(object):
    """Allows conversion between element name and atomic no.

    >>> t = PeriodicTable()
    >>> t.element[6]
    'C'
    >>> t.number['C']
    6
    >>> t.element[44]
    'Ru'
    >>> t.number['Au']
    79
    """
    
    def __init__(self):
        self.element = [None,
            'H', 'He',
            'Li', 'Be',
            'B', 'C', 'N', 'O', 'F', 'Ne',
            'Na', 'Mg',
            'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
            'K', 'Ca',
            'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
            'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
            'Rb', 'Sr',
            'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
            'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
            'Cs', 'Ba',
            'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
            'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
            'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
            'Fr', 'Ra',
            'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
            'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Uub']
        self.number = {}
        for i in range(1, len(self.element)):
            self.number[self.element[i]] = i


if __name__ == "__main__":
    import doctest, utils
    doctest.testmod(utils, verbose=False)
