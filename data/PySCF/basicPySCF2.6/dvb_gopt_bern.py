# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from pyscf import dft, gto
from pyscf.geomopt.berny_solver import optimize


def calculate():
    # This is DVB.
    mol = gto.M(
        atom="""
    C          1.43386       -0.00000        0.00000
    C         -1.43386        0.00000        0.00000
    C          0.69976       -1.21791        0.00000
    C         -0.69976        1.21791        0.00000
    C         -0.70271       -1.21960        0.00000
    C          0.70271        1.21960        0.00000
    H          1.24445       -2.17211        0.00000
    H         -1.24445        2.17211        0.00000
    H         -1.24514       -2.17312        0.00000
    H          1.24514        2.17312        0.00000
    C         -2.93083        0.05118        0.00000
    C          2.93083       -0.05118        0.00000
    H         -3.35520        1.06709        0.00000
    H          3.35520       -1.06709        0.00000
    C         -3.76833       -1.00028        0.00000
    C          3.76833        1.00028        0.00000
    H         -3.41193       -2.03669        0.00000
    H          3.41193        2.03669        0.00000
    H         -4.85625       -0.86201        0.00000
    H          4.85625        0.86201        0.00000
    """,
        basis="STO-3G",
        symmetry=True,
    )

    scf_steps = []

    def store_scf(_locals):
        scf_steps.append(
            {
                "cycle": _locals["cycle"],
                "e_tot": _locals["e_tot"],
                "norm_gorb": _locals["norm_gorb"],
                "conv_tol": _locals["conv_tol"],
                "conv_tol_grad": _locals["conv_tol_grad"],
            }
        )

    method = dft.RKS(mol)
    method.callback = store_scf
    method.xc = "b3lyp"
    method.kernel()

    # Discard the initial SCF step.
    scf_steps = []

    opt_steps = []

    def store_opt(_locals):
        opt_steps.append(
            {
                # Convergence criteria (not used yet).
                "gradient_max": _locals["optimizer"]._state.params["gradientmax"],
                "gradient_rms": _locals["optimizer"]._state.params["gradientrms"],
                "step_max": _locals["optimizer"]._state.params["stepmax"],
                "step_rms": _locals["optimizer"]._state.params["steprms"],
                # Values at this opt step.
                "energy": _locals["energy"],
                "gradients": _locals["gradients"],
                "coords": _locals["mol"].atom_coords("Angstrom"),
            }
        )

    optimize(method, callback=store_opt, maxsteps=100)

    # Unflatten the SCF cycle list.
    scf_map = [index for index, scf_step in enumerate(scf_steps) if scf_step["cycle"] == 0] + [
        len(scf_steps)
    ]
    scf_steps = [scf_steps[scf_map[index] : end] for index, end in enumerate(scf_map[1:])]

    return {"methods": [method], "scf_steps": scf_steps, "opt_steps": opt_steps}


if __name__ == "__main__":
    calculate()
