# Copyright (c) 2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
"""A reader for TREXIO files."""

from cclib.io import filereader
from cclib.parser.data import ccData
from cclib.parser.utils import PeriodicTable

import trexio
from numpy import array as np_array


class Trexio(filereader.Reader):
    """A reader for TREXIO files."""

    def __init__(self, source, *args, **kwargs):
        super().__init__(source, *args, **kwargs)

        self.pt = PeriodicTable()

    def parse(self):

        foo_file = trexio.File(self.filename, mode="r", backend=self.backend)

        # read data from TREXIO file but ignore errors for missing properties
        try:
            natom = trexio.read_nucleus_num(foo_file)
            # atomnos can be either an array of integers or the number of protons in the atom nuclei
            # is it always the nucleus_charge of TREXIO ?
            atomnos = trexio.read_nucleus_charge(foo_file)
            atomcoords = trexio.read_nucleus_coord(foo_file)

            atom_labels = trexio.read_nucleus_label(foo_file)

            nbasis = trexio.read_basis_num(foo_file)  # noqa: F841
            nmo = trexio.read_mo_num(foo_file)  # noqa: F841

            mo_sym = trexio.read_mo_symmetry(foo_file)

            aooverlaps = trexio.read_ao_1e_int_overlap(foo_file)  # noqa: F841

            mo_coeff = trexio.read_mo_coefficient(foo_file)

        except trexio.Error as e:
            if e.error == trexio.TREXIO_ATTR_MISSING or e.error == trexio.TREXIO_DSET_MISSING:
                print(f"Missing property in TREXIO file: pass \n Caused by: {e.__cause__}")
                pass
            else:
                raise

        foo_file.close()

        # so far, we treat TREXIO data as spin restricted for simplicity
        all_atomcoords = np_array(1, atomcoords)
        all_mo_sym = [mo_sym]  # noqa: F841
        all_mo_coefficients = [mo_coeff]  # noqa: F841

        # get atommasses from the list of labels and/or nuclear charges
        if atom_labels or atomnos:
            from periodictable import elements

            if atom_labels:
                masses = [elements.symbol(el).mass for el in atom_labels]  # noqa: F841
            elif atomnos:
                local_labels = [self.pt.element[Z] for Z in atomnos]
                masses = [elements.symbol(el).mass for el in local_labels]  # noqa: F841

        # TODO: convert from TREXIO data to cclib:
        # aonames (list of str), atombasis (list of int),
        # coreelectrons (number of core elements in each atom's pseudopotential)
        # gbasis (list of lists of tuples) - info about Gaussian basis functions per atom
        # homos (a single value corresponding to HOMO: 1 value for restricted, 2 values for unrestricted, e.g. [4])

        attributes = {
            "natom": natom,
            "atomnos": atomnos,
            "atomcoords": all_atomcoords,
            "metadata": {"comments": "TREXIO_TODO"},
        }

        self.data = ccData(attributes)

        return self.data
