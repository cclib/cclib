# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from pyscf import dft, gto
from pyscf.prop import infrared


def calculate():
    # This is DVB.
    mol = gto.M(
        atom="""
    C            -1.41006500     0.26933000     0.00000000
    C            -0.92062700    -1.06492700     0.00000000
    C             0.45713400    -1.32542900     0.00000000
    C             1.41006500    -0.26933000     0.00000000
    C             0.92062700     1.06492700     0.00000000
    H            -1.62749100    -1.90460200     0.00000000
    H             0.81321700    -2.36476300     0.00000000
    H             1.62749100     1.90460200     0.00000000
    C             2.86984500    -0.60347200     0.00000000
    C             3.89221900     0.26937300     0.00000000
    H             3.74294400     1.35516300     0.00000000
    H             3.08997800    -1.68167400     0.00000000
    H             4.93303400    -0.07474200     0.00000000
    C            -2.86984500     0.60347200     0.00000000
    H            -3.08997800     1.68167400     0.00000000
    C            -3.89221900    -0.26937300     0.00000000
    H            -3.74294400    -1.35516300     0.00000000
    H            -4.93303400     0.07474200     0.00000000
    C            -0.45713400     1.32542900     0.00000000
    H            -0.81321700     2.36476300     0.00000000
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

    irm = infrared.rhf.Infrared(method)
    irm.kernel()

    return {"methods": [irm], "scf_steps": [scf_steps]}
