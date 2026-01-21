# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from pyscf import dft, gto


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
        charge=1,
        spin=1,
    )

    scf_steps = []

    def store_intermediate(_locals):
        scf_steps.append(
            {
                "e_tot": _locals["e_tot"],
                "norm_gorb": _locals["norm_gorb"],
                "conv_tol": _locals["conv_tol"],
                "conv_tol_grad": _locals["conv_tol_grad"],
            }
        )

    method = dft.UKS(mol)
    method.callback = store_intermediate
    method.xc = "b3lyp"
    method.kernel()
    return {"methods": [method], "scf_steps": [scf_steps]}
