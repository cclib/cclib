"""
This file is an implementation of the Charge Model 5 that can be used to calculate CM5 Charges, dipole moment and quadrupole moment

It is based on this Fortran implementation: https://github.com/hokru/cm5charges
-----------------------------------------------
Points to be noted:
-----------------------------------------------
1) Only ORCA output files are supported, and "Print[P_Hirshfeld] 1" must be entered in the ORCA input file to produce Hirshfeld charges.
2) Single bond atomic radii are considered.
3) For Carbon, only sp3 atomic radius is considered.

-----------------------------------------------
Example on how to Use this class:
-----------------------------------------------
from cclib.method import CM5
from cclib.parser import ORCA

parser = ORCA("water_hpa.out") #test file found in /cclib/data/ORCA/basicORCA4.2/
data = parser.parse()

cm5 = CM5(data, 1.20)
print(cm5.cm5_charges())
print(cm5.dipole_moment())
print(cm5.quadrupole_moment())

------------------------------------------------
Output:
------------------------------------------------
[-0.3880329  0.2227029  0.16703  ]
(array([-3.7422439 , -3.34907393,  0.80222668]), 5.0856910330516865)
[[ 2.9476879   0.39972848 -0.47709812]
 [ 2.648237   -0.2199436  -0.65167883]
 [ 1.183168    1.256816   -0.11061444]]
"""

import logging

import numpy as np
import periodictable.covalent_radius as pt

from cclib.method.calculationmethod import Method


class CM5(Method):
    def __init__(self, data, fscale=1.20, progress=None, loglevel=logging.INFO, logname="Log"):
        super().__init__(data, progress, loglevel, logname)

        self.required_attrs = ("natom", "atomcoords", "atomnos")
        self.fscale = fscale
        self.atomradius = np.empty(119)
        self.atomradius[0] = 0.20

        for line in pt.CorderoPyykko.split("\n"):
            fields = line.split()
            if line[0] == "#":
                continue
            else:
                fields = line.split()
                Z = int(fields[0])

                if len(fields) < 6:
                    for _ in range(len(fields), 6):
                        fields.append("0.0")

                if fields[2] == "-":
                    rC = 0.0
                    drC = 0.0
                else:
                    sfields = fields[2].split("(")
                    rC = float(sfields[0])
                    drC = float(sfields[1].split(")", 1)[0]) / 100 if len(sfields) == 2 else 0.0

                r1 = float(fields[3])
                r1_avg = (rC + r1) / 2 if r1 is not 0.0 else rC
                r2 = float(fields[4])
                r3 = float(fields[5])

            self.atomradius[Z] = r1_avg

    def cm5_charges(self):
        """Returns an array with cm5 charges"""
        nat = self.data.natom
        qcm5 = np.empty(nat)
        z = self.data.atomnos
        xyz = self.data.atomcoords[-1]
        hirshfeld_charges = self.data.atomcharges["hirshfeld"]

        alpha = 2.474  # Angstrom

        for i in range(0, nat - 1):
            s = 0
            for j in range(0, nat - 1):
                if i != j:
                    rij = np.linalg.norm(np.subtract(xyz[i], xyz[j]))  # = nuclear.get_distances()
                    bij = np.exp(
                        -alpha * (rij - self.atomradius[z[i]] - self.atomradius[z[j]])
                    )  # eq.2
                    s += tij(z[i], z[j]) * bij
            qcm5[i] = hirshfeld_charges[i] + s
        return qcm5

    def dipole_moment(self):
        """Returns the dipole moment computed from CM5 atomic charges."""
        nat = self.data.natom
        dipole = np.empty(nat)
        cm5_charges = self.cm5_charges()
        for i in range(0, nat):
            for j in range(0, 2):
                dipole[j] += cm5_charges[i] * self.data.atomcoords[-1, j, i] * self.fscale
        dipole *= 4.802889778
        s = np.sqrt(dipole[0] ** 2 + dipole[1] ** 2 + dipole[2] ** 2)
        return dipole, s

    def quadrupole_moment(self):
        """Returns the quadrupole moment computed from CM5 atomic charges."""
        bohr = 0.52917726
        nat = self.data.natom
        quad = np.empty([nat, nat])
        cm5_charges = self.cm5_charges()
        for k in range(0, nat - 1):
            dx = self.data.atomcoords[-1, k, 0] / bohr
            dy = self.data.atomcoords[-1, k, 1] / bohr
            dz = self.data.atomcoords[-1, k, 2] / bohr
            quad[0, 0] += dx * dx * cm5_charges[k]
            quad[1, 1] += dy * dy * cm5_charges[k]
            quad[2, 2] += dz * dz * cm5_charges[k]
            quad[0, 1] += dx * dy * cm5_charges[k]
            quad[0, 2] += dx * dz * cm5_charges[k]
            quad[1, 2] += dy * dz * cm5_charges[k]
        return quad


def tij(i, j):
    #              NULL H       He       Li   Be      B        C        N        O        F        Ne
    dz = np.array(
        [
            0.0,
            0.0056,
            -0.1543,
            0.0,
            0.0333,
            -0.1030,
            -0.0446,
            -0.1072,
            -0.0802,
            -0.0629,
            -0.1088,  #   Na      Mg   Al       Si       P        S        Cl       Ar
            0.0184,
            0.0,
            -0.0726,
            -0.0790,
            -0.0756,
            -0.0565,
            -0.0444,
            -0.07676,  #   K       Ca                                                Zn   Ga       Ge       As       Se       Br   Kr
            0.0130,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -0.0557,
            -0.0533,
            -0.0399,
            -0.0313,
            0.0,
            0.0,  #                                                                                   I        Xe
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -0.0220,
            0.0,
        ]
    )

    #               H-C     H-N     H-O     C-N     C-O     N-O
    dzz = np.array([0.0502, 0.1747, 0.1671, 0.0556, 0.0234, -0.0346])

    # catch special cases
    if i == 1:
        if j == 6:
            tij = dzz[1]
        elif j == 7:
            tij = dzz[2]
        elif j == 8:
            tij = dzz[3]
    elif i == 6:
        if j == 7:
            tij = dzz[4]
        elif j == 8:
            tij = dzz[5]
        elif j == 1:
            tij = -dzz[1]
    elif i == 7:
        if j == 8:
            tij = dzz[6]
        elif j == 1:
            tij = -dzz[2]
        elif j == 6:
            tij = -dzz[4]
    elif i == 8:
        if j == 1:
            tij = -dzz[3]
        elif j == 6:
            tij = -dzz[5]
        elif j == 7:
            tij = -dzz[6]
    else:
        tij = dz[i] - dz[j]

    return tij
