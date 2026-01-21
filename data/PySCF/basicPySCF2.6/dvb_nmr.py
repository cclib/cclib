# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from pyscf import dft, gto
from pyscf.prop import nmr


def calculate():
    # This is DVB.
    mol = gto.M(
        atom="""
        C       -1.4152533224   0.2302217854   0.0000000000
        C        1.4152533224  -0.2302217854   0.0000000000
        C       -0.4951331558   1.3144608674   0.0000000000
        C        0.4951331558  -1.3144608674   0.0000000000
        C        0.8894090436   1.0909493743   0.0000000000
        C       -0.8894090436  -1.0909493743   0.0000000000
        H       -0.8795511985   2.3437343748   0.0000000000
        H        0.8795511985  -2.3437343748   0.0000000000
        H        1.5779041557   1.9450061275   0.0000000000
        H       -1.5779041557  -1.9450061275   0.0000000000
        C        2.8845844962  -0.5210893778   0.0000000000
        C       -2.8845844962   0.5210893778   0.0000000000
        H        3.1403356810  -1.5919605685   0.0000000000
        H       -3.1403356810   1.5919605685   0.0000000000
        C        3.8800428103   0.3822535424   0.0000000000
        C       -3.8800428103  -0.3822535424   0.0000000000
        H        3.6946765858   1.4624389570   0.0000000000
        H       -3.6946765858  -1.4624389570   0.0000000000
        H        4.9316453546   0.0711049543   0.0000000000
        H       -4.9316453546  -0.0711049543   0.0000000000
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

    nmrm = nmr.RKS(method)

    # Annoyingly, NMR coupling results are not saved in the NMR method object...
    nmrs = nmrm.kernel()

    return {"methods": [method, nmrm], "scf_steps": [scf_steps], "nmr_shielding": nmrs}
