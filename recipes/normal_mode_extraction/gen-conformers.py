#!/usr/bin/env python

from __future__ import print_function

import sys, os
import logging
import numpy as np
from cclib.io import ccopen
from nms import nmsgenerator
from data import element_masses, element_symbols


def nms_ts(n):
    """Normal Mode Temperature and Samples

    Args:
        n (int): The number of non-hydrogen atoms.

    Returns:
        temperature (float): thermal energy to populate normal normal_modes
        samples     (int):   number of samples to consider (ignored)
    """
    stable =  {1: (2000.0, 500),
               2: (1500.0, 450),
               3: (1000.0, 425),
               4: ( 600.0, 400),
               5: ( 600.0, 200),
               6: ( 600.0,  30),
               7: ( 600.0,  20),
               8: ( 450.0,   5),
               9: ( 450.0,   3)}
    if n in stable: return stable[n]
    else: return (400, 2)

# read through a bunch of files on the command-line
for arg in sys.argv[1:]:
    # read the file with cclib
    file = ccopen(arg)
    file.logger.setLevel(logging.ERROR)
    molecule = file.parse()

    # the atom labels (here as atomic numbers)
    labels = molecule.atomnos
    try:
        # not all cclib parsers will populate atommasses
        # but this is more accurate, since it would include isotopes
        mass = np.asarray(molecule.atommasses).repeat(3)
    except AttributeError:
        mass = np.asarray([element_masses[i] for i in labels]).repeat(3)

    # get the last geometry from a geometry optimization
    # should end up as an ndarray (N, 3)
    xyz = np.asarray(molecule.atomcoords[-1]).reshape(-1,3)

    # vibrational frequencies (remove any zeros, e.g., translations, rotations)
    freq = np.asarray(molecule.vibfreqs)
    freq = freq[freq !=0.0] # should now be an array of size (3N-6,)

    # get flattened normal modes and reshape to (3N, 3N)
    # ensure they're in (CartesianCoordinate, ModeNumber)
    nmo = np.asarray(molecule.vibdisps).reshape(-1, xyz.shape[0]*3) # [3N, 3N]
    # orca filters out translations and rotations
    # they appear as modes with all-zero components
    # so please remove them
    nmo = nmo[np.array([not (i == 0.0).all() for i in nmo])]

    # compute mass weighted force constants for each of 3N-6 normal modes
    # shapes: [3N-6] * ([3N-6, 3N] * [3N,]).sum(axis=1) = [3N-6,]
    fcc = freq * freq * (nmo * nmo * mass).sum(axis=1) / 16.9744 / 10e4

    # get the temperature and number of samples for sampling
    t, s = nms_ts(np.count_nonzero(labels != 1))
    g = nmsgenerator(xyz, nmo, fcc, labels, t, 0)
    # generate 100 displaced geometries into new XYZ files
    filename, extension = os.path.splitext(arg)
    for i in range(100):
        energy, rxyz = g.get_random_structure(100)
        # print an xyz file
        xyzstr = '{}\n{}\n'.format(len(rxyz), energy)
        for l, x in zip(labels, rxyz):
            xyzstr += '{} {} {} {}\n'.format(element_symbols[l], *x)

        new_name = '{}.xyz'.format(i)
        with open(new_name, 'w') as f:
            f.write(xyzstr)
