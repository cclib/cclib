# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from pyscf import dft, gto
from pyscf.prop import infrared


def calculate():
    # This is DVB.
    mol = gto.M(
        atom="""
    C         -1.41525        0.23022        0.00000
    C          1.41525       -0.23022        0.00000
    C         -0.49513        1.31446        0.00000
    C          0.49513       -1.31446        0.00000
    C          0.88941        1.09095        0.00000
    C         -0.88941       -1.09095        0.00000
    H         -0.87955        2.34373        0.00000
    H          0.87955       -2.34373        0.00000
    H          1.57790        1.94501        0.00000
    H         -1.57790       -1.94501        0.00000
    C          2.88458       -0.52109        0.00000
    C         -2.88458        0.52109        0.00000
    H          3.14034       -1.59196        0.00000
    H         -3.14034        1.59196        0.00000
    C          3.88004        0.38225        0.00000
    C         -3.88004       -0.38225        0.00000
    H          3.69468        1.46244        0.00000
    H         -3.69468       -1.46244        0.00000
    H          4.93164        0.07111        0.00000
    H         -4.93164       -0.07111        0.00000
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

    irm = infrared.rks.Infrared(method)
    irm.kernel()

    return {"methods": [irm], "scf_steps": [scf_steps]}
