# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in PySCF (https://github.com/pyscf/pyscf)."""

import functools
import itertools
import logging
from typing import Any, Dict, List, Optional, Tuple
import numpy

from cclib.parser.data import ccData
from cclib.parser.utils import PeriodicTable, convertor, find_package

import numpy as np

l_sym2num = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4}


class MissingAttributeError(Exception):
    pass


_found_pyscf = find_package("pyscf")
if _found_pyscf:
    import pyscf.cc.ccsd
    import pyscf.data.elements
    import pyscf.gto
    import pyscf.hessian.thermo
    import pyscf.mp.mp2
    import pyscf.scf.hf
    import pyscf.tdscf.rhf

    # This is an optional install.
    try:
        import warnings

        warnings.filterwarnings(
            action="ignore", category=UserWarning, message="Since PySCF-2.3, B3LYP"
        )
        warnings.filterwarnings(
            action="ignore", category=UserWarning, message=r"Module [\w-]+ is under testing"
        )
        warnings.filterwarnings(
            action="ignore", category=UserWarning, message=r"Module [\w-]+ is not fully tested"
        )
        import pyscf.prop as pyscf_prop

    except ModuleNotFoundError:
        pyscf_prop = None


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
    mol = pyscf.gto.Mole(
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
    # mo_energies = convertor(np.asarray(ccdata.moenergies), "eV", "hartree")
    mo_energies = np.asarray(ccdata.moenergies)
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


def makecclib(
    # pyscf.lib.StreamObject is the ultimate base class I believe, but not sure
    # if we can give this as a type hint as we don't know if PySCF is even available...
    *methods: Any,
    ccsd_t: Optional[float] = None,
    scf_steps: List[List[Dict[str, float]]] = [],
    opt_steps: List[Dict[str, Any]] = [],
    opt_failed: bool = False,
    nmr_coupling: List[Tuple[Tuple[float, float, float],Tuple[float, float, float],Tuple[float, float, float]]]
) -> ccData:
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
        opt_failed - set to true when parsing an optimisation that failed to converge
    """
    _check_pyscf(_found_pyscf)
    # What is our level of theory?
    scf, mp, cc, hess, freq, et, nmr = None, None, None, None, None, [], None

    mols = {method.mol for method in methods if hasattr(method, "mol")}
    if len(mols) > 1:
        logger = logging.getLogger("cclib")
        logger.warning("not all PySCF methods are operating on the same molecule")

    for method in methods:
        # Assume the method is a base CC method (SCF, MP, CC etc.) unless we can prove otherwise.
        base_method = method

        if pyscf_prop is not None and isinstance(method, pyscf_prop.infrared.rhf.Infrared):
            # Vibrational frequencies.
            base_method = method.base
            hess = method.mf_hess
            freq = method

        elif isinstance(method, pyscf.tdscf.rhf.TDBase):
            # Excited states.
            # TODO: What about multiple excited states?
            base_method = method._scf
            et.append(method)

        if isinstance(base_method, (pyscf.scf.hf.SCF, pyscf.scf.hf.KohnShamDFT)):
            scf = base_method

        elif isinstance(base_method, pyscf.mp.mp2.MP2):
            mp = base_method
            scf = base_method._scf

        elif isinstance(base_method, pyscf.cc.ccsd.CCSDBase):
            cc = base_method
            scf = base_method._scf
        
        elif pyscf_prop is not None and isinstance(method, pyscf_prop.nmr.rhf.NMR):
            nmr = base_method

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
        opt_failed=opt_failed,
        nmr=nmr,
        nmr_coupling=nmr_coupling
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
    opt_failed=False,
    nmr=None,
    nmr_coupling=[]
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
        opt_failed - set to true when parsing an optimisation that failed to converge

    """
    _check_pyscf(_found_pyscf)
    attributes = {}

    mol = scf.mol
    ptable = PeriodicTable()

    # Metadata.
    attributes["metadata"] = {
        "basis_set": mol.basis
        if isinstance(mol.basis, str)
        else ", ".join(set(mol.basis.values())),
        "input_file_contents": mol.tostring(),
        "legacy_package_version": pyscf.__version__,
        "memory_available": mol.max_memory * 1000 * 1000,
        "num_cpu": pyscf.lib.num_threads(),
        "methods": [],
        "package": "PySCF",
        "package_version": pyscf.__version__,
        "symmetry_detected": mol.topgroup.lower(),
        "symmetry_used": mol.groupname.lower(),
    }

    # Solvent.
    if hasattr(scf, "with_solvent"):
        attributes["metadata"]["solvent_model"] = scf.with_solvent.method.replace("-", "")
        attributes["metadata"]["solvent_params"] = {"epsilon": scf.with_solvent.eps}

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

    coord_convertor = functools.partial(convertor, fromunits="Angstrom", tounits="bohr")
    bohr_coords = [
        [list(map(coord_convertor, ang_coords)) for ang_coords in opt_step]
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
        attributes["mpenergies"] = [[convertor(mp.e_tot, "hartree", "eV")]]
        attributes["metadata"]["methods"].append("MP2")

    if cc:
        # We have to manually add in the CCSD(T) correction energy.
        attributes["ccenergies"] = [
            (convertor(cc.e_tot + ccsd_t, "hartree", "eV"))
            if ccsd_t
            else convertor(cc.e_tot, "hartree", "eV")
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
            opt_e = attributes["ccenergies"]

        elif mp:
            attributes["mpenergies"] = [
                [convertor(step["energy"], "hartree", "eV")] for step in opt_steps
            ]
            opt_e = attributes["mpenergies"]

        else:
            attributes["scfenergies"] = [
                convertor(step["energy"], "hartree", "eV") for step in opt_steps
            ]
            opt_e = attributes["scfenergies"]

        attributes["optstatus"] = [ccData.OPT_UNKNOWN for _ in attributes["atomcoords"]]
        attributes["optstatus"][0] = ccData.OPT_NEW

        # Only energy for now.
        attributes["geovalues"] = [[energy] for energy in opt_e]

        # This is obviously not ideal, but due to the way PySCF handles optimisations (ie, it doesn't handle them directly)
        # it's difficult to come up with something better...
        if not opt_failed:
            attributes["optdone"] = [len(attributes["atomcoords"]) - 1]
            attributes["optstatus"][-1] = ccData.OPT_DONE

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
    attributes["moments"] = [
        origin,
        scf.dip_moment(scf.mol, scf_density_matrix, origin=origin, verbose=2),
    ]

    # Quadrupole moments are new (introduced ~ August 2024) and not yet widely available.
    # So far, PySCF only supports the traceless quadrupole (the raw version is calculated
    # internally but is not exposed), I propose we hold off on this until both become available.
    # if hasattr(scf, "quad_moment"):
    #     quad_moment_matrix = scf.quad_moment(
    #             scf.mol, scf_density_matrix, origin=origin, unit="DebyeAngstrom", verbose=2
    #         )
    #     attributes["moments"].append(
    #         [quad_moment_matrix[0,1], # XX
    #          quad_moment_matrix[0,1], # XY
    #          quad_moment_matrix[0,2], # XZ
    #          quad_moment_matrix[1,1], # YY
    #          quad_moment_matrix[1,2], # YZ
    #          quad_moment_matrix[2,2], # ZZ
    #         ]
    #     )

    # Mulliken.
    mulliken_pop, mulliken_charges = scf.mulliken_pop(
        scf.mol, scf_density_matrix, s=attributes["aooverlaps"], verbose=2
    )
    # PySCF describes this as 'Mulliken population analysis, based on meta-Lowdin AOs', is this what we want?
    lowdin_pop, lowdin_charges = scf.mulliken_meta(
        scf.mol, scf_density_matrix, s=attributes["aooverlaps"], verbose=2
    )
    attributes["atomcharges"] = {"mulliken": mulliken_charges, "lowdin": lowdin_charges}

    # Excited states.
    if len(et) > 0:
        # PySCF tracks convergence for each state which is great.
        attributes["metadata"]["success"] = all(
            itertools.chain(*(etmethod.converged for etmethod in et))
        )

        # In cclib 1.x, 'energies' are actually expected to be cm-1. In 2.x, this will change to Hartree.
        attributes["etenergies"] = list(
            itertools.chain(*(convertor(etmethod.e, "hartree", "wavenumber") for etmethod in et))
        )
        attributes["etoscs"] = list(
            itertools.chain(*(etmethod.oscillator_strength(gauge="length") for etmethod in et))
        )  # or do we want velocity?
        # et.analyse() prints real symmetries, so they must be available somewhere...

        attributes["etsyms"] = []
        for etmethod in et:
            if isinstance(etmethod, pyscf.tdscf.rks.TDA):
                attributes["metadata"]["excited_states_method"] = "TDA"

            elif isinstance(etmethod, pyscf.tdscf.rks.TDDFT):
                attributes["metadata"]["excited_states_method"] = "TD-DFT"

            elif isinstance(etmethod, pyscf.tdscf.rhf.TDA):
                attributes["metadata"]["excited_states_method"] = "CIS"

            elif isinstance(etmethod, pyscf.tdscf.rhf.TDHF):
                attributes["metadata"]["excited_states_method"] = "RPA"

            else:
                attributes["metadata"]["excited_states_method"] = type(etmethod).__name__

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

            attributes["etsyms"].extend([et_mult for _ in range(len(etmethod.e))])

        # Orbital contributions.
        # In PySCF, occupied and virtual orbital indices are stored separately.
        attributes["etsecs"] = []

        for index in range(len(attributes["etenergies"])):
            attributes["etsecs"].append([])

            # Assuming x are excitations, y are de-excitations.
            # y is ignored for now.
            x, y = list(itertools.chain(*(etmethod.xy for etmethod in et)))[index]

            def norm_xy(x, y):
                norm_factor = 1.0 / np.sqrt(np.linalg.norm(x) ** 2 - np.linalg.norm(y) ** 2)
                return (x * norm_factor, y * norm_factor)

            if not scf.istype("UHF"):
                # The coefficients of x (and presumably also y) are normalised to 0.5,
                # but we expect 1.0
                # Renormalize.
                x = norm_xy(x, y)[0]

                # Flatten the x matrix.
                # The first index is the occupied orbital, the second is the virtual (both 0 indexed):
                o_indices, v_indices = np.where(x)
                for occupied, virtual in zip(o_indices, v_indices):
                    attributes["etsecs"][index].append(
                        [(occupied, 0), (nocc[0] + virtual, 0), x[occupied, virtual]]
                    )

            else:
                x = [x[0], x[1]]
                x[0] = norm_xy(x[0], y[0])[0]
                x[1] = norm_xy(x[1], y[1])[0]

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
        pass
        # TODO: check this conversion
        # cclib wants a rank 2 array, pyscf gives us a rank 4 array?
        # Dimensions appear to be natoms x natoms x 3 x 3
        #
        # Units may also be mismatched.
        # attributes["hessian"] = hess.de.transpose(0,2,1,3).reshape(attributes["natom"], *3, attributes["natom"]*3)

    if freq:
        # Freq data, a dict with 'freq_error', 'freq_au', 'freq_wavenumber', 'norm_mode', 'reduced_mass',
        # 'vib_temperature', 'force_const_au', 'force_const_dyne'
        attributes["vibfreqs"] = freq.vib_dict["freq_wavenumber"]

        # Convert imaginary frequencies to negative ones.
        # Adapted from https://pyscf.org/_modules/pyscf/hessian/thermo.html
        if np.iscomplexobj(attributes["vibfreqs"]):
            attributes["vibfreqs"] = attributes["vibfreqs"].real - abs(attributes["vibfreqs"].imag)

        # TODO: Check these units.
        attributes["vibfconsts"] = freq.vib_dict["force_const_dyne"]
        attributes["vibrmasses"] = freq.vib_dict["reduced_mass"]

        # freq.de probably contains something useful, but not sure what.
        # The Infrared module is sadly less well documented that the other code.
        # apparently de is a "3-dim tensor: (atom number, atom coordinate components, dipole components)"
        # maybe one of these is vibdisps?

        if hasattr(freq, "ir_inten"):
            # TODO: Units?
            attributes["vibirs"] = freq.ir_inten
        
    # NMR
    if len(nmr_coupling):
        # nmr_coupling is a list, nmrtensors is a dict
        attributes['nmrtensors'] = {
            index: {
                "isotropic": numpy.mean(numpy.linalg.eigvals(value)),
                "total": value
            } for index, value in enumerate(nmr_coupling)
        }

    return ccData(attributes)


del find_package
