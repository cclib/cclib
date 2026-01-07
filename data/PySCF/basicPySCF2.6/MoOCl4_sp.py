# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from pyscf import dft, gto


def calculate():
    mol = gto.M(
        atom="""
        Mo         0.00001        0.00000        0.29568
        O          0.00000        0.00000        1.98538
        Cl         2.42022       -0.00000       -0.41621
        Cl         0.00000       -2.42024       -0.41620
        Cl        -2.42025       -0.00000       -0.41620
        Cl         0.00000        2.42024       -0.41620
        """,
        basis={"Mo": "lanl2dz", "Cl": "lanl2dz", "O": "6-31G*"},
        ecp={"Mo": "lanl2dz", "Cl": "lanl2dz"},
        symmetry=True,
        charge=-2,
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
    return {"methods": [method], "scf_steps": [scf_steps]}


if __name__ == "__main__":
    calculate()
