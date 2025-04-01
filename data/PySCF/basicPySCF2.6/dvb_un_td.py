# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from pyscf import dft, gto, tddft


def calculate():
    # This is DVB.
    mol = gto.M(
        atom="""
    C             0.000000    0.000000   -0.643075
    O             0.000000    0.000000    0.482306
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

    method = dft.UKS(mol)
    method.callback = store_intermediate
    method.xc = "b3lyp"
    method.kernel()

    # Now excited states.
    singlets = tddft.TDA(method)
    singlets.nstates = 5
    singlets.singlet = True
    singlets.kernel()
    assert all(singlets.converged)

    triplets = tddft.TDA(method)
    triplets.nstates = 5
    triplets.singlet = False
    triplets.kernel()
    assert all(triplets.converged)

    return {"methods": [singlets, triplets], "scf_steps": [scf_steps]}


if __name__ == "__main__":
    methods = calculate()
