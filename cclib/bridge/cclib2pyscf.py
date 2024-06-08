# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in PySCF (https://github.com/pyscf/pyscf)."""

from cclib.parser.data import ccData
from cclib.parser.utils import PeriodicTable, convertor, find_package

import numpy as np

l_sym2num = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4}


class MissingAttributeError(Exception):
    pass


_found_pyscf = find_package("pyscf")
if _found_pyscf:
    from pyscf import gto


def _check_pyscf(found_pyscf):
    if not found_pyscf:
        raise ImportError("You must install `pyscf` to use this function")


def makepyscf(data, charge=0, mult=1):
    """Create a Pyscf Molecule."""
    _check_pyscf(_found_pyscf)
    required_attrs = {"atomcoords", "atomnos"}
    missing = [x for x in required_attrs if not hasattr(data, x)]
    if missing:
        missing = " ".join(missing)
        raise MissingAttributeError(
            f"Could not create pyscf molecule due to missing attribute: {missing}"
        )
    mol = gto.Mole(
        atom=[[f"{data.atomnos[i]}", data.atomcoords[-1][i]] for i in range(data.natom)],
        unit="Angstrom",
        charge=charge,
        multiplicity=mult,
    )
    inputattr = data.__dict__
    pt = PeriodicTable()
    if "gbasis" in inputattr:
        basis = {}  # object for internal PySCF format
        uatoms, uatoms_idx = np.unique(data.atomnos, return_index=True)  # find unique atoms
        for idx, i in enumerate(uatoms_idx):
            curr_atom_basis = data.gbasis[i]
            for jdx, j in enumerate(curr_atom_basis):
                curr_l = j[0]
                curr_e_prim = j[1]
                new_list = [l_sym2num[f"{curr_l}"]]
                new_list += curr_e_prim
                if f"{pt.element[uatoms[idx]]}" not in basis:
                    basis[f"{pt.element[uatoms[idx]]}"] = [new_list]
                else:
                    basis[f"{pt.element[uatoms[idx]]}"].append(new_list)
        mol.basis = basis
        mol.cart = True
    return mol


def makepyscf_mos(ccdata, mol):
    """
    Returns pyscf formatted MO properties from a cclib object.
    Parameters
    ---
    ccdata: cclib object
        cclib object from parsed output
    mol: pyscf Molecule object
       molecule object that must contain the mol.basis attribute

    Returns
    ----
    mo_coeff : n_spin x nmo x nao ndarray
        molecular coeffcients, unnormalized according to pyscf standards
    mo_occ : array
        molecular orbital occupation
    mo_syms : array
       molecular orbital symmetry labels
    mo_energies: array
        molecular orbital energies in units of Hartree
    """
    inputattrs = ccdata.__dict__
    mo_energies = convertor(np.asarray(ccdata.moenergies), "eV", "hartree")
    if "mocoeffs" in inputattrs:
        mol.build()
        s = mol.intor("int1e_ovlp")
        if np.shape(ccdata.mocoeffs)[0] == 1:
            mo_coeffs = np.einsum("i,ij->ij", np.sqrt(1 / s.diagonal()), ccdata.mocoeffs[0].T)
            mo_occ = np.zeros(ccdata.nmo)
            mo_occ[: ccdata.homos[0] + 1] = 2
            if hasattr(ccdata, "mosyms"):
                mo_syms = ccdata.mosyms
            else:
                mo_syms = np.full_like(ccdata.moenergies, "A", dtype=str)

        elif np.shape(ccdata.mocoeffs)[0] == 2:
            mo_coeff_a = np.einsum("i,ij->ij", np.sqrt(1 / s.diagonal()), ccdata.mocoeffs[0].T)
            mo_coeff_b = np.einsum("i,ij->ij", np.sqrt(1 / s.diagonal()), ccdata.mocoeffs[1].T)
            mo_occ = np.zeros((2, ccdata.nmo))
            mo_occ[0, : ccdata.homos[0] + 1] = 1
            mo_occ[1, : ccdata.homos[1] + 1] = 1
            mo_coeffs = np.array([mo_coeff_a, mo_coeff_b])
            if hasattr(ccdata, "mosyms"):
                mo_syms = ccdata.mosyms
            else:
                mo_syms = np.full_like(ccdata.moenergies, "A", dtype=str)
    return mo_coeffs, mo_occ, mo_syms, mo_energies


def makecclib(method, etmethod=None) -> ccData:
    """Create cclib attributes and return a ccData from a PySCF method object.

    The method object should naturally have already performed some sort of
    calculation.

    Inputs:
        atoms - an instance of ASE `Atoms`
        popname - population analysis to use for atomic partial charges and
            atomic spin densities. Molecular charge and multiplicity are
            evaluated from them.
    """
    _check_pyscf(_found_pyscf)
    attributes = {}

    mol = method.mol
    ptable = PeriodicTable()

    # Atoms.
    attributes["atomcoords"] = mol.atom_coords("Angstrom")
    attributes["atomnos"] = [ptable.number[element] for element in mol.elements]
    # attributes["atomcharges"] = mol.atom_charges() # is this the right type of atom charge?
    attributes["atommasses"] = mol.atom_mass_list(isotope_avg=True)

    attributes["charge"] = mol.charge
    attributes["multiplicity"] = mol.multiplicity
    attributes["coreelectrons"] = [
        mol.atom_nelec_core(i) for i in range(0, len(attributes["atomnos"]))
    ]

    # Excited states.
    if etmethod:
        attributes["etenergies"] = etmethod.e
        attributes["etoscs"] = etmethod.oscillator_strength(
            gauge="length"
        )  # or do we want velocity?
        # etmethod.analyse() prints real symmetries, so they must be available somehwere...
        attributes["etsyms"] = [
            "Singlet" if etmethod.singlet else "Triplet"
            for i in range(0, len(attributes["atomnos"]))
        ]

        # Orbital contributions.
        # In PySCF, occupied and virtual orbital indices are stored separately.
        attributes["etsecs"] = []
        nocc = np.count_nonzero(method.mo_occ == 2)
        for index in range(0, len(etmethod.e)):
            attributes["etsecs"].append([])

            # Assuming x are excitations, y are de-excitations.
            # y is ignored for now.
            x, y = etmethod.xy[index]

            # Flatten the x matrix.
            # The first index is the occupied orbital, the second is the virtual (both 0 indexed):
            o_indices, v_indices = np.where(x)
            for occupied, virtual in zip(o_indices, v_indices):
                attributes["etsecs"][index].append(
                    [(occupied, 0), (nocc + virtual, 0), x[occupied, virtual]]
                )

    return ccData(attributes)


del find_package
