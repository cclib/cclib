#!/usr/bin/env python

"""qchem_make_opt_input_from_opt.py: Make an input file for a Q-Chem
geometry optimization based on the last possible geometry from a
Q-Chem geometry optimization; this effectively 'restarts' the geometry
with a new filename.

The script assumes the output file being read from is called
'*opt(\d*).out', where 'opt' might be followed by a number. The script
will write an input file called '*opt(\d*)+1.in', with the previous
number incremented by one.
"""


import os.path
import re
from collections import OrderedDict

import cclib
from cclib.parser.utils import PeriodicTable


def make_file_iterator(filename):
    """Return an iterator over the contents of the given file name."""
    # pylint: disable=C0103
    with open(filename) as f:
        contents = f.read()
    return iter(contents.splitlines())


def getargs():
    """Get command-line arguments."""

    import argparse

    # pylint: disable=C0103
    parser = argparse.ArgumentParser()

    parser.add_argument("outputfilename", nargs="+")

    parser.add_argument("--fragment", action="store_true")

    args = parser.parse_args()

    return args


def parse_user_input(outputfilename):
    """Parse the $rem section in the repeated 'User input:' section of the
    output.

    The reason we do it this way rather than with shell tools is to
    handle any $section more easily and in a case-insensitive manner.
    """

    user_input = dict()

    outputfile = make_file_iterator(outputfilename)

    line = ""
    while "User input:" not in line:
        line = next(outputfile)
    line = next(outputfile)
    assert "----" in line
    line = next(outputfile)
    while "--------------------------------------------------------------" not in line:
        if line.strip() == "":
            pass
        elif line[0] == "$" and line.strip().lower() != "$end":
            section_header = line[1:].lower()
            user_input[section_header] = []
        elif line.strip().lower() == "$end":
            user_input[section_header] = "\n".join(user_input[section_header])
        else:
            user_input[section_header].append(line)
        line = next(outputfile)

    return user_input


def parse_fragments_from_molecule(molecule):
    """Given a $molecule section (without the $ lines), identify the
    charges and multiplicities of each fragment and the zero-based indices
    for the starting atom of each fragment.
    """

    charges = []
    multiplicities = []
    start_indices = []
    it = iter(molecule.splitlines())
    line = next(it)
    # sys_charge, sys_multiplicity = line.split()
    counter = 0
    # Gather the charges, spin multiplicities, and starting positions
    # of each fragment.
    for line in it:
        if "--" in line:
            line = next(it)
            charge, multiplicity = line.split()
            charges.append(charge)
            multiplicities.append(multiplicity)
            start_indices.append(counter)
        else:
            counter += 1
    assert len(charges) == len(multiplicities) == len(start_indices)

    return charges, multiplicities, start_indices


def form_molecule_section_from_fragments(
    elements, geometry, charges, multiplicities, start_indices
):
    """Form the Q-Chem $molecule section containing the charge,
    multiplicity, and atomic symbols and coordinates for multiple
    fragments.

    Returns a list that will need to be joined with newlines.
    """

    assert len(charges) == len(multiplicities) == (len(start_indices) + 1)

    s = "{:3s} {:15.10f} {:15.10f} {:15.10f}"
    # The first elements of the charge and multiplicity lists are for
    # the supersystem (whole molecule).
    molecule_section = ["{} {}".format(charges[0], multiplicities[0])]

    from itertools import count

    for (charge, multiplicity, idx_iter) in zip(charges[1:], multiplicities[1:], count(0)):
        molecule_section.append("--")
        molecule_section.append("{} {}".format(charge, multiplicity))
        idx_start = start_indices[idx_iter]
        try:
            idx_end = start_indices[idx_iter + 1]
        except IndexError:
            idx_end = len(elements)
        for element, coords in zip(elements[idx_start:idx_end], geometry[idx_start:idx_end]):
            molecule_section.append(s.format(element, *coords))

    return molecule_section


def form_molecule_section(elements, geometry, charge, multiplicity):
    """Form the Q-Chem $molecule section containing the charge,
    multiplicity, and atomic symbols and coordinates.

    Returns a list that will need to be joined with newlines.
    """

    s = "{:3s} {:15.10f} {:15.10f} {:15.10f}"
    molecule_section = ["{} {}".format(charge, multiplicity)]

    for (
        element,
        coords,
    ) in zip(elements, geometry):
        molecule_section.append(s.format(element, *coords))

    return molecule_section


if __name__ == "__main__":

    args = getargs()

    pt = PeriodicTable()

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
        numstr = re.search(r"opt(\d*)", outputfilename).groups()[0]
        if numstr == "":
            optnum = 2
        else:
            optnum = int(numstr) + 1
        inputfilename = re.sub(r"opt\d*", "opt{}".format(optnum), outputfilename)
        inputfilename = inputfilename.replace(".out", ".in")
        inputfilename = os.path.basename(inputfilename)

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
        user_input["molecule"] = "\n".join(molecule_section)

        with open(inputfilename, "w") as fh:
            for section_header in user_input:
                fh.write("${}\n".format(section_header))
                fh.write(user_input[section_header])
                fh.write("\n$end\n\n")

        print(inputfilename)
