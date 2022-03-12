#!/usr/bin/env python

"""cclib_plot_geometry_convergence.py: Given a computational chemistry
logfile for a geometry optimization, plot how the convergence
parameters change over the course of the optimization.
"""


import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import argparse
import os.path

import cclib
from cclib.io import ccopen
from cclib.parser import utils


def main():
    """The main routine!"""

    parser = argparse.ArgumentParser()

    parser.add_argument("compchemfilename", nargs="+")
    parser.add_argument("--scaling-energy-change", type=float, default=10.0)

    args = parser.parse_args()
    compchemfilenames = args.compchemfilename

    for compchemfilename in compchemfilenames:

        stub = os.path.splitext(compchemfilename)[0]

        job = ccopen(compchemfilename)
        data = job.parse()

        fig, ax = plt.subplots()

        if type(job) == cclib.parser.qchemparser.QChem:

            scfenergies = [
                utils.convertor(scfenergy, "eV", "hartree") for scfenergy in data.scfenergies
            ]
            gradients = [geovalue[0] for geovalue in data.geovalues]
            displacements = [geovalue[1] for geovalue in data.geovalues]
            energy_changes = [
                (geovalue[2] * args.scaling_energy_change) for geovalue in data.geovalues
            ]

            # If this isn't true, something funny happened during the
            # parsing, so fail out.
            assert len(scfenergies) == len(gradients)

            steps = list(range(1, len(scfenergies) + 1))

            # ax.plot(steps, scfenergies, label='SCF energy')
            ax.plot(steps, gradients, label="max gradient")
            ax.plot(steps, displacements, label="max displacement")
            ax.plot(steps, energy_changes, label="energy change")

            ax.set_title(stub)
            ax.set_xlabel("optimization step #")
            ax.set_xlim((steps[0], steps[-1]))

        elif type(job) == cclib.parser.orcaparser.ORCA:

            pass

        else:
            pass

        ax.legend(loc="best", fancybox=True)

        fig.savefig(stub + ".pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
