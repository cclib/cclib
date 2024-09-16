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
    C         -1.41012        0.26944        0.00000
    C         -0.92066       -1.06479       -0.00000
    C          0.45706       -1.32539       -0.00000
    C          1.41012       -0.26944        0.00000
    C          0.92066        1.06479        0.00000
    H         -1.62752       -1.90447       -0.00000
    H          0.81304       -2.36476       -0.00000
    H          1.62752        1.90447        0.00000
    C          2.86996       -0.60369       -0.00000
    C          3.89207        0.26944       -0.00000
    H          3.74231        1.35518       -0.00000
    H          3.09015       -1.68186       -0.00000
    H          4.93303       -0.07425       -0.00000
    C         -2.86996        0.60369       -0.00000
    H         -3.09015        1.68186       -0.00000
    C         -3.89207       -0.26944       -0.00000
    H         -3.74231       -1.35518       -0.00000
    H         -4.93303        0.07425       -0.00000
    C         -0.45706        1.32539       -0.00000
    H         -0.81304        2.36476       -0.00000
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
