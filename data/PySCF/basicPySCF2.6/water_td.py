# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from pyscf import dft, gto, tddft


def calculate():
    mol = gto.M(
        atom="""
    O         -0.00000       -0.11916        0.00000
    H         -0.79065        0.47664       -0.00000
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

    method = dft.RKS(mol)
    method.callback = store_intermediate
    method.xc = "b3lyp"
    method.kernel()

    # Now excited states.
    # TD-DFT proved extremely unstable, maybe a bug in PySCF?
    tdm = tddft.TD(method)
    tdm.nstates = 5
    tdm.singlet = True
    tdm.kernel()

    assert all(tdm.converged)

    return {"methods": [tdm], "scf_steps": [scf_steps]}


if __name__ == "__main__":
    methods = calculate()

    from cclib.bridge.cclib2pyscf import makecclib

    makecclib(*methods["methods"])
