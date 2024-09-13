# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in PySCF (https://github.com/pyscf/pyscf)."""

import functools
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
    import pyscf.data.elements
    import pyscf.hessian.thermo
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


def makecclib(method, ccsd_t=None, scf_steps=[], opt_steps=[]) -> ccData:
    """Create cclib attributes and return a ccData from a PySCF calculation.

    PySCF calculation results are stored in method objects, with each object representing a different part of the
    calculation. For example, you may have separate methods for the HF and MP parts of an MP2 calculation.
    This function will try to intelligently guess which methods you have available.
    If you already know which method objects correspond to which levels of theory, use cclibfrommethods() instead.

    Inputs:
        method - an instance of a PySCF method object that has already performed a calculation
        ccsd_t - for CCSD(T), the (T) correction energy to the CC energy
        scf_steps - a list (per optimisation step) of lists (per SCF cycle) of dictionaries containing the results (and targets)
          of each SCF cycle. Each dict should contain at least 'e_tot', 'norm_gorb', 'conv_tol', and 'conv_tol_grad'.
        opt_steps - for an optimisation, a list of dictionaries containing the 'energy', 'gradients', and 'coords'
          at each opt step (including the last)
    """
    _check_pyscf(_found_pyscf)
    # What is our level of theory?
    scf, mp, cc, hess, freq, et = None, None, None, None, None, []

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
        scf=scf,
        mp=mp,
        cc=cc,
        ccsd_t=ccsd_t,
        hess=hess,
        freq=freq,
        scf_steps=scf_steps,
        opt_steps=opt_steps,
        et=et,
    )


def cclibfrommethods(
    # TODO: Types
    scf,
    mp=None,
    cc=None,
    ccsd_t=None,
    et=[],
    hess=None,
    freq=None,
    scf_steps=[],
    opt_steps=[],
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
        scf_steps - a list (per optimisation step) of lists (per SCF cycle) of dictionaries containing the results (and targets)
          of each SCF cycle. Each dict should contain at least 'e_tot', 'norm_gorb', 'conv_tol', and 'conv_tol_grad'.
        opt_steps  - for an optimisation, a list of dictionaries containing the 'energy', 'gradients', and 'coords'
          at each opt step (including the last)

    """
    _check_pyscf(_found_pyscf)
    attributes = {}

    mol = scf.mol
    ptable = PeriodicTable()

    # TODO: A sanity check that all the supplied methods use the same mol object?

    # Metadata.
    attributes["metadata"] = {
        # TODO: What if using a non-standard basis set?
        #'basis_set': mol.basis,
        "basis_set": mol.basis
        if isinstance(mol.basis, str)
        else ", ".join(set(mol.basis.values())),
        "input_file_contents": mol.tostring(),
        "legacy_package_version": pyscf.__version__,
        "methods": [],
        "package": "PySCF",
        "package_version": pyscf.__version__,
        "symmetry_detected": mol.topgroup.lower(),
        "symmetry_used": mol.groupname.lower(),
    }

    # Atoms.
    if not opt_steps:
        attributes["atomcoords"] = [mol.atom_coords("Angstrom")]

    else:
        attributes["atomcoords"] = [step["coords"] for step in opt_steps]

    attributes["natom"] = len(mol.atom_coords("Angstrom"))
    attributes["atomnos"] = [ptable.number[element] for element in mol.elements]
    # attributes["atomcharges"] = mol.atom_charges() # is this the right type of atom charge?

    # mol.atom_mass_list() does what it sounds like, and you can choose between average or exact mass with
    # atom_mass_list(isotope_avg = True | False). However, isotope_avg = False (what we want) strangely only returns
    # an integer mass. This is probably a PySCF bug.
    # Fortunately, PySCF ships a table of single isotope masses (COMMON_ISOTOPE_MASSES), so we'll just use that.
    attributes["atommasses"] = [
        pyscf.data.elements.COMMON_ISOTOPE_MASSES[atom_no] for atom_no in attributes["atomnos"]
    ]

    converter = functools.partial(convertor, fromunits="Angstrom", tounits="bohr")
    bohr_coords = [
        [list(map(converter, ang_coords)) for ang_coords in opt_step]
        for opt_step in attributes["atomcoords"]
    ]

    attributes["rotconsts"] = [
        pyscf.hessian.thermo.rotation_const(np.array(attributes["atommasses"]), opt_coords)
        for opt_coords
        # rotation_const expects atom positions in Bohr nor Angstrom.
        in bohr_coords
    ]

    attributes["charge"] = mol.charge
    attributes["mult"] = mol.multiplicity
    attributes["coreelectrons"] = [
        mol.atom_nelec_core(i) for i in range(0, len(attributes["atomnos"]))
    ]

    # Atomic orbitals and gbasis
    attributes["gbasis"] = []
    for atom_index in range(len(mol.atom_coords("Angstrom"))):
        # Shell types are found mol.bas_angular()
        # Basis exponents are found in mol.bas_exp()
        # Contraction coefficients are found in mol.bas_ctr_coeff()
        # Each function take an index (the orbital). The latter two return a list (of equal length)
        # for each contracted GTO.
        # However, bas_ctr_coeff() returns a list for each contracted orbital (which is normally of length 1).
        # Presumably this is to support SP type orbitals (with one exponent and multiple coefficients), but
        # it's not clear if this type of orbital is actually supported in PySCF?
        #
        # bas_angular returns the quantum number (index). We can convert this to a label using ANGULARMAP
        # (hopefully this is stable).
        #
        # The orbital indices for each atom can be found in mol.atom_shell_ids(), which takes
        # the atom index as argument.
        atom_basis = []
        for basis_index in mol.atom_shell_ids(atom_index):
            atom_basis.append(
                (
                    # Orbital label (S, P, D etc.)
                    pyscf.lib.parameters.ANGULAR[mol.bas_angular(basis_index)].upper(),
                    list(
                        zip(
                            # Exponent.
                            mol.bas_exp(basis_index),
                            # Coefficient, unpacked into a single item.
                            # TODO: What do if multiple coefficients are present?
                            [coeff_list[0] for coeff_list in mol.bas_ctr_coeff(basis_index)],
                        )
                    ),
                )
            )
        attributes["gbasis"].append(atom_basis)

    attributes["nbasis"] = mol.nao

    # mol.ao_labels() and aonames have almost the same format, just need to tweak a bit.
    attributes["aonames"] = [
        f"{atom_label}{atom_index + 1}_{shell.upper()}{xyz.upper()}"
        for atom_index, atom_label, shell, xyz in mol.ao_labels(False)
    ]

    # Build atombasis from ao_labels()
    ao_map = [atom_index for atom_index, atom_label, shell, xyz in mol.ao_labels(False)]
    attributes["atombasis"] = [
        [
            basis_index
            for basis_index, basis_atom_index in enumerate(ao_map)
            if basis_atom_index == atom_index
        ]
        for atom_index in range(attributes["natom"])
    ]

    # get_ovlp() looks like it has the correct format for aooverlaps already.
    attributes["aooverlaps"] = scf.get_ovlp()

    # Total energies.
    attributes["scfenergies"] = [convertor(scf.e_tot, "hartree", "eV")]
    attributes["metadata"]["success"] = scf.converged
    attributes["metadata"]["methods"].append("DFT" if hasattr(scf, "xc") else "HF")
    if hasattr(scf, "xc"):
        attributes["metadata"]["functional"] = scf.xc

    if mp:
        attributes["mpenergies"] = [convertor(mp.e_tot, "hartree", "eV")]
        attributes["metadata"]["success"] = mp.converged
        attributes["metadata"]["methods"].append("MP2")

    if cc:
        # We have to manually add in the CCSD(T) correction energy.
        attributes["ccenergies"] = [
            convertor((cc.e_tot + ccsd_t) if ccsd_t else cc.e_tot, "hartree", "eV")
        ]
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
    if opt_steps:
        if cc:
            attributes["ccenergies"] = [
                convertor(step["energy"], "hartree", "eV") for step in opt_steps
            ]

        elif mp:
            attributes["mpenergies"] = [
                convertor(step["energy"], "hartree", "eV") for step in opt_steps
            ]

        else:
            attributes["scfenergies"] = [
                convertor(step["energy"], "hartree", "eV") for step in opt_steps
            ]

    if scf_steps:
        # scf_steps contains all the info we need for scfvalues and scftargets, just need to unpack it.
        # TODO: scfvalues/scftargets should probably use a dicts to describe each convergence criteria,
        # but this is not specific to PySCF.
        attributes["scfvalues"] = []
        attributes["scftargets"] = []
        for opt_step in scf_steps:
            if len(opt_step) > 0:
                attributes["scftargets"].append(
                    (opt_step[0]["conv_tol"], opt_step[0]["conv_tol_grad"])
                )
                attributes["scfvalues"].append([])

            for scf_step in opt_step:
                attributes["scfvalues"][-1].append((scf_step["e_tot"], scf_step["norm_gorb"]))

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

    # Dipole moments.
    #
    scf_density_matrix = scf.make_rdm1(scf.mo_coeff, scf.mo_occ)
    # PySCF uses 0,0,0 as the origin, but this might not correspond to the centre of mass?
    origin = (0, 0, 0)
    attributes["moments"] = [origin, scf.dip_moment(scf.mol, scf_density_matrix, origin=origin)]

    # Quadrupole moments are new (introduced ~ August 2024) and not yet widely available.
    if hasattr(scf, "quad_moment"):
        attributes["moments"].append(
            scf.quad_moment(scf.mol, scf_density_matrix, origin=origin, unit="DebyeAngstrom")
        )

    # Mulliken.
    mulliken_pop, mulliken_charges = scf.mulliken_pop(
        scf.mol, scf_density_matrix, s=attributes["aooverlaps"]
    )
    # PySCF describes this as 'Mulliken population analysis, based on meta-Lowdin AOs', is this what we want?
    lowdin_pop, lowdin_charges = scf.mulliken_meta(
        scf.mol, scf_density_matrix, s=attributes["aooverlaps"]
    )
    attributes["atomcharges"] = {"mulliken": mulliken_charges, "lowdin": lowdin_charges}

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
