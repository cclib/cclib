# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from pyscf import gto, scf


def calculate():
    mol = gto.M(
        atom="""
    O         -0.00000       -0.11916        0.00000
    H         -0.79065        0.47664        0.00000
    H          0.79065        0.47664       -0.00000
    """,
        basis="STO-3G",
        symmetry=True,
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

    method = scf.HF(mol)
    method.callback = store_intermediate
    method.kernel()

    mp = method.MP2()
    mp.kernel()
    return {"methods": [mp], "scf_steps": [scf_steps]}
