# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in ASE (https://wiki.fysik.dtu.dk/ase/)."""

import numpy as np

from cclib.parser.data import ccData
from cclib.parser.utils import find_package

_found_ase = find_package("ase")
if _found_ase:
    from ase import Atoms, units
    from ase.io.trajectory import Trajectory
    from ase.calculators.calculator import PropertyNotImplementedError


def _check_ase(found_ase):
    if not found_ase:
        raise ImportError("You must install `ase` to use this function")


def makease(
    atomcoords, atomnos, atomcharges=None, atomspins=None, atommasses=None
):
    """Create an ASE Atoms object from cclib attributes.

    ASE requires atomic partial charges and atomic spin densities rather than
    molecular charge and multiplicity, so we follow how other interfaces have
    done (e.g., MOPAC, Gaussian and XTB) and require atomcharges and atomspins,
    or leave undefined.

    Inputs:
        atomcoords - two dimensional array-like with atomic coordinates.
        atomnos - one dimensional array-like with atomic numbers.
        atomcharges - one dimensional array-like with atomic charges.
        atomspins - one dimensional array-like with atomic spin densities.
        atommasses - one dimensional array-like with atomic masses.
    """
    _check_ase(_found_ase)
    return Atoms(
        positions=atomcoords,
        numbers=atomnos,
        masses=atommasses,
        charges=atomcharges,
        magmoms=atomspins,
    )


def write_trajectory(filename, ccdata, popname="mulliken", index=None):
    """Write an ASE Trajectory object from a ccData object.

    We try to write the following properties: atomcoords, atomnos, atomcharges,
    atomspins, atommasses, scfenergies, grads, moments (dipole only) and
    freeenergy. No charge or mult is written since ASE calculates them from
    atomcharges and atomspins.

    When a program do a final single point calculation at the end of an
    optimization (e.g., ORCA), we'll have unequal number of grads and
    atomcoords. This means the last geometry in the written .traj will lack
    forces.

    Some unit conversions are done here, but things are converted back in
    read_trajectory/makecclib.

    Inputs:
        filename - path to traj file to be written.
        ccdata - an instance of ccData.
        popname - population analysis to use for atomic partial charges and
            atomic spin densities. Molecular charge and multiplicity are
            evaluated from them.
        index - sequence of integers indicating which atomcoords indices should
            be exported. By default, all are exported.
    """
    _check_ase(_found_ase)
    traj = Trajectory(filename, "w")
    for i, atomcoords in enumerate(ccdata.atomcoords):
        if index is not None and i not in index:
            continue

        atomspins = None
        if hasattr(ccdata, "atomspins"):
            atomspins = ccdata.atomspins[popname]
        atoms = makease(
            atomcoords,
            ccdata.atomnos,
            ccdata.atomcharges[popname],
            atomspins,
            ccdata.atommasses,
        )

        properties = {}
        if hasattr(ccdata, "scfenergies"):
            properties.update({"energy": ccdata.scfenergies[i]})
        if hasattr(ccdata, "grads"):
            try:
                properties.update(
                    {"forces": -ccdata.grads[i] * units.Hartree / units.Bohr}
                )
            except IndexError:
                pass

        if i == len(ccdata.atomcoords) - 1:  # last geometry
            if hasattr(ccdata, "moments"):
                properties.update({"dipole": ccdata.moments[1] * units.Bohr})
            if hasattr(ccdata, "free_energy"):
                properties.update(
                    {"free_energy": ccdata.freeenergy * units.Hartree}
                )

        traj.write(atoms, **properties)


def read_trajectory(filename):
    """Read an ASE Trajectory object and return a ccData object.

    The returned object has everything write_trajectory writes, plus natom,
    charge, mult and temperature.

    The following properties are taken from the last frame: atomnos,
    atomcharges, atomspins, atommasses, moments, freeenergy and temperature.
    charge, mult and natom also represent the last frame, since they depend on
    other propertes read from the last frame.

    Bear in mind that ASE calculates temperature from the kinetic energy, so
    anything "static" (which includes anything cclib parses) will have zero
    temperature.

    Inputs:
        filename - path to traj file to be read.
    """
    _check_ase(_found_ase)
    attributes = {"atomcoords": [], "scfenergies": [], "grads": []}

    for atoms in Trajectory(filename, "r"):
        ccdata = makecclib(atoms)
        attributes["atomcoords"].append(ccdata.atomcoords[-1])
        if hasattr(ccdata, "scfenergies"):
            attributes["scfenergies"].append(ccdata.scfenergies[-1])
        if hasattr(ccdata, "grads"):
            attributes["grads"].append(ccdata.grads[-1])

    # ccdata is now last frame
    attributes["atomnos"] = ccdata.atomnos
    attributes["atomcharges"] = ccdata.atomcharges
    attributes["atomspins"] = ccdata.atomspins
    attributes["atommasses"] = ccdata.atommasses

    if hasattr(ccdata, "moments"):
        attributes["moments"] = ccdata.moments
    if hasattr(ccdata, "freeenergy"):
        attributes["freeenergy"] = ccdata.freeenergy

    # remove if empty
    if not attributes["scfenergies"]:
        del attributes["scfenergies"]
    if not attributes["grads"]:
        del attributes["grads"]

    # extra stuff we can't write in write_trajectory
    attributes["temperature"] = ccdata.temperature
    attributes["charge"] = ccdata.charge
    attributes["mult"] = ccdata.mult
    attributes["natom"] = ccdata.natom

    return ccData(attributes)


def makecclib(atoms, popname="mulliken"):
    """Create cclib attributes and return a ccData from an ASE Atoms object.

    Available data (such as forces/gradients and potential energy/free
    energy) is assumed to be from SCF (see
    https://wiki.fysik.dtu.dk/ase/ase/atoms.html#adding-a-calculator).

    Bear in mind that ASE calculates temperature from the kinetic energy, so
    anything "static" (which includes anything cclib parses) will have zero
    temperature.

    Inputs:
        atoms - an instance of ASE `Atoms`
        popname - population analysis to use for atomic partial charges and
            atomic spin densities. Molecular charge and multiplicity are
            evaluated from them.
    """
    _check_ase(_found_ase)
    attributes = {
        "atomcoords": np.array([atoms.get_positions()]),
        "atomnos": atoms.get_atomic_numbers(),
        "atommasses": atoms.get_masses(),
        "natom": atoms.get_global_number_of_atoms(),
    }

    try:
        attributes["atomcharges"] = {popname: atoms.get_charges()}
    except (PropertyNotImplementedError, RuntimeError):
        attributes["atomcharges"] = {popname: atoms.get_initial_charges()}
    try:
        attributes["atomspins"] = {popname: atoms.get_magnetic_moments()}
    except (PropertyNotImplementedError, RuntimeError):
        attributes["atomspins"] = {
            popname: atoms.get_initial_magnetic_moments()
        }

    # the following is how ASE determines charge and multiplicity from initial
    # charges and initial magnetic moments in its Gaussian interface
    # (https://gitlab.com/ase/ase/-/blob/a26bda2160527ca7afc0135c69e4367a5bc5a264/ase/io/gaussian.py#L105)
    attributes["charge"] = attributes["atomcharges"][popname].sum()
    attributes["mult"] = attributes["atomspins"][popname].sum() + 1

    try:
        attributes["scfenergies"] = np.array([atoms.get_potential_energy()])
    except RuntimeError:
        pass
    try:
        attributes["grads"] = (
            -np.array([atoms.get_forces()]) * units.Bohr / units.Hartree
        )
    except RuntimeError:
        pass
    try:
        attributes["moments"] = [
            atoms.get_center_of_mass(),
            atoms.get_dipole_moment() / units.Bohr,
        ]
    except RuntimeError:
        pass
    try:
        attributes["freeenergy"] = (
            atoms.get_potential_energy(force_consistent=True) / units.Hartree
        )
    except RuntimeError:
        pass

    attributes["temperature"] = atoms.get_temperature()

    return ccData(attributes)


del find_package
