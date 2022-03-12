#!/usr/bin/env python

"""qchem_make_freq_input_from_opt.py: Make an input file for a Q-Chem
frequency calculation based on the last possible geometry from a
Q-Chem geometry optimization.

The script assumes the output file being read from is called
'*opt(\d*)*.out', where 'opt' might be followed by a number. The
script will write an input file called '*freq*.in'.
"""


import os.path
import re
from collections import OrderedDict
from copy import deepcopy

from qchem_make_opt_input_from_opt import (
    form_molecule_section,
    form_molecule_section_from_fragments,
    make_file_iterator,
    parse_fragments_from_molecule,
    parse_user_input,
)

import cclib
from cclib.parser.utils import PeriodicTable


def template_input_freq(molecule, **rem_keywords):
    """The template for a input file that performs a frequency calculation."""
    avoid_these_keywords = ("jobtype",)
    rem_pieces = []
    for k in rem_keywords:
        if k not in avoid_these_keywords:
            rem_pieces.append(" {} = {}".format(k, rem_keywords[k]))
    rem_str = "\n".join(rem_pieces)
    # Return a string, which is the contents of the new input file,
    # with the '{}' fields appropriately replaced.
    return """$rem
 jobtype = freq
{rem_str}
$end

$molecule
{molecule}
$end
""".format(
        rem_str=rem_str, molecule=molecule
    )


def clean_up_rem(rem, do_fragment=False):
    """Make sure our $rem section is stringent enough for frequency
    calculations.
    """

    # These are the minimum values we'd like for frequency
    # calculations.
    min_scf_convergence = 9
    min_thresh = 12

    newrem = deepcopy(rem)

    if "thresh" in newrem:
        if int(newrem["thresh"]) < min_thresh:
            newrem["thresh"] = min_thresh
    if "scf_convergence" in newrem:
        if int(newrem["scf_convergence"]) < min_scf_convergence:
            newrem["scf_convergence"] = min_scf_convergence

    # If doing a no charge transfer fragment calculation, we can't use
    # analytic Hessians, so switch to numerical ones.
    if do_fragment:
        if "frgm_lpcorr" in newrem:
            if int(newrem["frgm_lpcorr"]) == 0:
                newrem["ideriv"] = 1
    # If we aren't doing a fragment calculation at all, certain
    # keywords need to be removed.
    else:
        keywords_to_remove = (
            "frgm_method",
            "frgm_lpcorr",
        )
        for k in keywords_to_remove:
            if k in newrem:
                del newrem[k]

    return newrem


def getargs():
    """Get command-line arguments."""

    import argparse

    # pylint: disable=C0103
    parser = argparse.ArgumentParser()

    parser.add_argument("outputfilename", nargs="+")

    parser.add_argument("--fragment", action="store_true")

    args = parser.parse_args()

    return args


def parse_rem_section(outputfilename):
    """Parse the $rem section in the repeated 'User input:' section of the
    output.
    """

    outputfile = make_file_iterator(outputfilename)

    rem = []

    line = ""
    while line.strip().lower() != "$rem":
        line = next(outputfile)
    line = next(outputfile)
    while "$end" not in line:
        sline = line.split()
        k = sline[0].replace("=", "").lower()
        v = sline[-1].replace("=", "").lower()
        rem.append((k, v))
        line = next(outputfile)

    return OrderedDict(rem)


if __name__ == "__main__":

    args = getargs()

    pt = PeriodicTable()
    # Format string template for the XYZ section.
    s = "{:3s} {:15.10f} {:15.10f} {:15.10f}"

    for outputfilename in args.outputfilename:

        job = cclib.io.ccopen(outputfilename)
        assert isinstance(job, cclib.parser.qchemparser.QChem)
        try:
            data = job.parse()
        # this is to deal with the Q-Chem parser not handling
        # incomplete SCF cycles properly
        except StopIteration:
            print("no output made: StopIteration in {}".format(outputfilename))
            continue

        # Determine the name of the file we're writing.
        assert outputfilename.endswith(".out")
        inputfilename = re.sub(r"opt\d*", "freq", outputfilename)
        inputfilename = inputfilename.replace(".out", ".in")
        inputfilename = os.path.basename(inputfilename)

        rem = parse_rem_section(outputfilename)

        # Make sure our $rem section is up to snuff for frequency
        # calculations.
        rem = clean_up_rem(rem, do_fragment=args.fragment)

        user_input = parse_user_input(outputfilename)

        # Form the atomic symbols and coordinates for each atom in
        # $molecule.
        element_list = [pt.element[Z] for Z in data.atomnos]
        last_geometry = data.atomcoords[-1]
        if args.fragment:
            charges, multiplicities, start_indices = parse_fragments_from_molecule(
                user_input["molecule"]
            )
            charges.insert(0, data.charge)
            multiplicities.insert(0, data.mult)
            molecule_section = form_molecule_section_from_fragments(
                element_list, last_geometry, charges, multiplicities, start_indices
            )
        else:
            molecule_section = form_molecule_section(
                element_list, last_geometry, data.charge, data.mult
            )
        molecule = "\n".join(molecule_section)

        with open(inputfilename, "w") as inputfile:
            inputfile.write(template_input_freq(molecule, **rem))

        print(inputfilename)
