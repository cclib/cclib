# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from pyscf import dft, gto, tddft


def calculate():
    # This is DVB.
    mol = gto.M(
        atom="""
    C            -1.41011800     0.26944500     0.00000000
    C            -0.92066200    -1.06478500     0.00000000
    C             0.45705800    -1.32538900     0.00000000
    C             1.41011800    -0.26944500     0.00000000
    C             0.92066200     1.06478500     0.00000000
    H            -1.62752200    -1.90447000     0.00000000
    H             0.81304100    -2.36475900     0.00000000
    H             1.62752200     1.90447000     0.00000000
    C             2.86996000    -0.60368800     0.00000000
    C             3.89207000     0.26944500     0.00000000
    H             3.74231400     1.35517500     0.00000000
    H             3.09015200    -1.68185900     0.00000000
    H             4.93302900    -0.07424700     0.00000000
    C            -2.86996000     0.60368800     0.00000000
    H            -3.09015200     1.68185900     0.00000000
    C            -3.89207000    -0.26944500     0.00000000
    H            -3.74231400    -1.35517500     0.00000000
    H            -4.93302900     0.07424700     0.00000000
    C            -0.45705800     1.32538900     0.00000000
    H            -0.81304100     2.36475900     0.00000000
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

    method = dft.RKS(mol)
    method.callback = store_intermediate
    method.xc = "b3lyp"
    method.kernel()

    # Now excited states.
    # TD-DFT proved extremely unstable, maybe a bug in PySCF?
    tdm = tddft.TDA(method)
    tdm.nstates = 5
    tdm.singlet = True
    tdm.kernel()

    assert all(tdm.converged)

    return {"methods": [tdm], "scf_steps": [scf_steps]}


if __name__ == "__main__":
    methods = calculate()
