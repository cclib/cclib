# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in PySCF (https://github.com/pyscf/pyscf)."""

import itertools

from cclib.parser.data import ccData
from cclib.parser.utils import PeriodicTable, convertor, find_package

import numpy as np

l_sym2num = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4}


class MissingAttributeError(Exception):
    pass


_found_pyscf = find_package("pyscf")
if _found_pyscf:
    import pyscf
    import pyscf.cc.ccsd
    import pyscf.mp.mp2
    import pyscf.prop.infrared.rhf
    import pyscf.scf.hf
    import pyscf.tdscf.rhf
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


def makecclib(method, ccsdt=None, opt_steps=None) -> ccData:
    """Create cclib attributes and return a ccData from a PySCF calculation.

    PySCF calculation results are stored in method objects, with each object representing a different part of the
    calculation. For example, you may have separate methods for the HF and MP parts of an MP2 calculation.
    This function will try to intelligently guess which methods you have available.
    If you already know which method objects correspond to which levels of theory, use cclibfrommethods() instead.

    Inputs:
        method - an instance of a PySCF method object that has already performed a calculation
        ccsdt - for CCSD(T), the (T) correction energy to the CC energy
        opt_steps - for an optimisation, a list of dictionaries containing the 'energy', 'gradients', and 'coords'
          at each opt step (including the last)
    """
    _check_pyscf(_found_pyscf)
    # What is our level of theory?
    scf, mp, cc, hess, freq, et = None, None, None, None, None, None

    # Assume the method is a base CC method (SCF, MP, CC etc.) unless we can prove otherwise.
    base_method = method

    if isinstance(method, pyscf.prop.infrared.rhf.Infrared):
        # Vibrational frequencies.
        base_method = method.base
        hess = method.mf_hess
        freq = method

    elif isinstance(method, pyscf.tdscf.rhf.TDBase):
        # Excited states.
        # TODO: What about multiple excited states?
        base_method = method._scf
        et = [method]

    if isinstance(base_method, pyscf.scf.hf.SCF) or isinstance(
        base_method, pyscf.scf.hf.KohnShamDFT
    ):
        scf = base_method

    elif isinstance(base_method, pyscf.mp.mp2.MP2):
        mp = base_method
        scf = base_method._scf

    elif isinstance(base_method, pyscf.cc.ccsd.CCSDBase):
        cc = base_method
        scf = base_method._scf

    else:
        # Panic.
        raise ValueError(
            f"Could not determine level of theory of base method '{type(base_method.__name__)}'"
        )

    return cclibfrommethods(
        scf=scf, mp=mp, cc=cc, ccsdt=ccsdt, hess=hess, freq=freq, opt=opt_steps, et=et
    )


def cclibfrommethods(
    # TODO: Types
    scf,
    mp=None,
    cc=None,
    ccsd_t=None,
    et=None,
    hess=None,
    freq=None,
    opt=None,
) -> ccData:
    """Create cclib attributes and return a ccData from a PySCF method object.

    Inputs:
        scf - an instance of a PySCF SCF method (RHF, UHF, RKS, UKS etc.).
        mp - an instance of a PySCF MPn method.
        cc - an instance of a PySCF CC method.
        ccsd_t - CCSD(T) correction to the CCSD energy.
        et - a list of instances of PySCF excited states methods (TD, TDA etc.).
              both singlet and triplet calculating methods can be passed simultaneously.
        hess - an instance of a PySCF hessian method.
        freq - an instance of an Infrared analysis class from pyscf.prop.infrared.*
        opt  - for an optimisation, a list of dictionaries containing the 'energy', 'gradients', and 'coords'
          at each opt step (including the last)

    """
    _check_pyscf(_found_pyscf)
    attributes = {}

    mol = scf.mol
    ptable = PeriodicTable()

    # TODO: A sanity check that all the supplied methods use the same mol object?

    # Metadata.
    attributes["metadata"] = {
        "package": "PySCF",
        "package_version": pyscf.__version__,
        # TODO: What if using a non-standard basis set?
        #'basis_set': mol.basis,
        "basis_set": ", ".join(set(mol.basis.values())),
        "methods": [],
    }

    # Atoms.
    if not opt:
        attributes["atomcoords"] = [mol.atom_coords("Angstrom")]

    else:
        attributes["atomcoords"] = [step["coords"] for step in opt]

    attributes["natoms"] = len(attributes["atomcoords"])
    attributes["atomnos"] = [ptable.number[element] for element in mol.elements]
    # attributes["atomcharges"] = mol.atom_charges() # is this the right type of atom charge?
    attributes["atommasses"] = mol.atom_mass_list(isotope_avg=True)

    attributes["charge"] = mol.charge
    attributes["mult"] = mol.multiplicity
    attributes["coreelectrons"] = [
        mol.atom_nelec_core(i) for i in range(0, len(attributes["atomnos"]))
    ]

    # Total energies.
    attributes["scfenergies"] = [scf.e_tot]
    attributes["metadata"]["success"] = scf.converged
    attributes["metadata"]["methods"].append("DFT" if hasattr(scf, "xc") else "HF")
    if hasattr(scf, "xc"):
        attributes["metadata"]["functional"] = scf.xc

    if mp:
        attributes["mpenergies"] = [mp.e_tot]
        attributes["metadata"]["success"] = mp.converged
        attributes["metadata"]["methods"].append("MP2")

    if cc:
        attributes["ccenergies"] = [cc.e_tot]
        attributes["metadata"]["success"] = cc.converged
        if cc.cc2:
            ccmethod = "CC2"

        elif ccsd_t:
            ccmethod = "CCSD(T)"

        else:
            ccmethod = "CCSD"

        attributes["metadata"]["methods"].append(ccmethod)

    # It's not immediately clear from the intermediate optimisation steps what level of theory the energy corresponds to,
    # we have to be smart based on what we asked for.
    if opt:
        if cc:
            attributes["ccenergies"] = [step["energy"] for step in opt]

        elif mp:
            attributes["mpenergies"] = [step["energy"] for step in opt]

        else:
            attributes["scfenergies"] = [step["energy"] for step in opt]

    # Orbitals.
    if scf.istype("UHF"):
        attributes["metadata"]["unrestricted"] = True

        nocc = [np.count_nonzero(scf.mo_occ[0] != 0), np.count_nonzero(scf.mo_occ[1] != 0)]
        attributes["moenergies"] = [
            convertor(scf.mo_energy[0], "hartree", "eV"),
            convertor(scf.mo_energy[1], "hartree", "eV"),
        ]
        attributes["homos"] = [nocc[0] - 1, nocc[1] - 1]
        attributes["mocoeffs"] = scf.mo_coeff

    else:
        attributes["metadata"]["unrestricted"] = False

        nocc = [np.count_nonzero(scf.mo_occ != 0)]
        attributes["moenergies"] = [convertor(scf.mo_energy, "hartree", "eV")]
        attributes["homos"] = [nocc[0] - 1]
        # Orbital coeffs.
        attributes["mocoeffs"] = [scf.mo_coeff]

    attributes["nmo"] = len(attributes["moenergies"][0])

    # Orbital symmetries.
    # if scf.mol.symmetry:
    if hasattr(scf, "get_orbsym"):
        # Symmetry labels are in scf.mol.irrep_name, symmetry ids are in scf.mol.irrep_id
        symmetry_mapping = {
            symm_id: symm_name for symm_id, symm_name in zip(scf.mol.irrep_id, scf.mol.irrep_name)
        }

        orbsyms = scf.get_orbsym() if scf.istype("UHF") else [scf.get_orbsym()]

        attributes["mosyms"] = []
        for alpha_beta in range(0, (len(attributes["moenergies"]))):
            attributes["mosyms"].append(
                [symmetry_mapping[symm_id] for symm_id in orbsyms[alpha_beta]]
            )

    # Excited states.
    if len(et) > 0:
        # PySCF tracks convergence for each state which is great.
        attributes["metadata"]["success"] = all(
            itertools.chain(*(etmethod.converged for etmethod in et))
        )

        # In cclib 1.x, 'energies' are actually expected to be cm-1. In 2.x, this will change to Hartree.
        attributes["etenergies"] = [
            convertor(hartree, "hartree", "wavenumber")
            for hartree in itertools.chain(*(etmethod.e for etmethod in et))
        ]
        attributes["etoscs"] = list(
            itertools.chain(*(etmethod.oscillator_strength(gauge="length") for etmethod in et))
        )  # or do we want velocity?
        # et.analyse() prints real symmetries, so they must be available somewhere...

        attributes["etsyms"] = []
        for etmethod in et:
            # From from pyscf/pyscf/tdscf/rhf.py
            # import pyscf.symm
            # orbsym = scf.get_orbsym()
            # x_sym = pyscf.symm.direct_prod(orbsym[mo_occ==2], orbsym[mo_occ==0], mol.groupname)

            # Need to be careful in interpreting the multiplicity of the excited state.
            if etmethod.singlet and mol.multiplicity != 2:
                et_mult = "Singlet"

            elif etmethod.singlet:
                et_mult = "Doublet"

            else:
                et_mult = "Triplet"

            attributes["etsyms"].extend([et_mult for i in range(len(etmethod.e))])

        # Orbital contributions.
        # In PySCF, occupied and virtual orbital indices are stored separately.
        attributes["etsecs"] = []

        for index in range(len(attributes["etenergies"])):
            attributes["etsecs"].append([])

            # Assuming x are excitations, y are de-excitations.
            # y is ignored for now.
            x, y = list(itertools.chain(*(etmethod.xy for etmethod in et)))[index]

            if not scf.istype("UHF"):
                # Flatten the x matrix.
                # The first index is the occupied orbital, the second is the virtual (both 0 indexed):
                o_indices, v_indices = np.where(x)
                for occupied, virtual in zip(o_indices, v_indices):
                    attributes["etsecs"][index].append(
                        [(occupied, 0), (nocc[0] + virtual, 0), x[occupied, virtual]]
                    )

            else:
                # Flatten the x matrix.
                # The first index is the occupied orbital, the second is the virtual (both 0 indexed):
                # Alpha -> alpha
                o_indices, v_indices = np.where(x[0])
                for occupied, virtual in zip(o_indices, v_indices):
                    attributes["etsecs"][index].append(
                        [(occupied, 0), (nocc[0] + virtual, 0), x[0][occupied, virtual]]
                    )

                # Beta -> Beta
                o_indices, v_indices = np.where(x[1])
                for occupied, virtual in zip(o_indices, v_indices):
                    attributes["etsecs"][index].append(
                        [(occupied, 1), (nocc[1] + virtual, 1), x[1][occupied, virtual]]
                    )

    # Hessian/frequencies
    if hess:
        # This doesn't seem to exist?
        # attributes["metadata"]["success"] = hess.converged
        pass
        # TODO: Don't know enough to work out what this is
        # cclib wants a rank 2 array, pyscf gives us a rank 4 array?
        # Dimensions appear to be natoms x natoms x 3 x 3
        # attributes["hessian"] = hess.de

    if freq:
        # Freq data, a dict with 'freq_error', 'freq_au', 'freq_wavenumber', 'norm_mode', 'reduced_mass',
        # 'vib_temperature', 'force_const_au', 'force_const_dyne'
        # TODO: This can include imaginary numbers, convert to 'negative' frequencies?
        attributes["vibfreqs"] = freq.vib_dict["freq_wavenumber"]
        # TODO: Check these units.
        attributes["vibfconsts"] = freq.vib_dict["force_const_dyne"]
        attributes["vibrmasses"] = freq.vib_dict["reduced_mass"]

        if hasattr(freq, "ir_inten"):
            # TODO: Units?
            attributes["vibirs"] = freq.ir_inten

    return ccData(attributes)


del find_package
