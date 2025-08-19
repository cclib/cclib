# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for Molden writer."""

import os

import cclib
from cclib.io.filewriter import MissingAttributeError
from cclib.io.moldenwriter import MoldenReformatter, round_molden

import pytest

__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")
__testdir__ = __filedir__


class MOLDENTest:
    def test_missing_attribute_error(self) -> None:
        """Check if MissingAttributeError is raised as expected."""
        fpath = os.path.join(__datadir__, "data/GAMESS/basicGAMESS-US2018/dvb_un_sp.out")
        required_attrs = ["atomcoords", "atomnos", "natom"]
        for attr in required_attrs:
            data = cclib.io.ccread(fpath)
            delattr(data, attr)

            # Molden files cannot be written if required attrs are missing.
            with pytest.raises(MissingAttributeError):
                cclib.io.moldenwriter.MOLDEN(data)

    def test_atoms_section_size(self) -> None:
        """Check if size of Atoms section is equal to expected."""
        fpath = os.path.join(__datadir__, "data/GAMESS/basicGAMESS-US2018/dvb_un_sp.out")
        data = cclib.io.ccread(fpath)
        writer = cclib.io.moldenwriter.MOLDEN(data)
        # Check size of Atoms section.
        assert len(writer._coords_from_ccdata(-1)) == data.natom

    def test_atoms_section_size_with_ghost(self) -> None:
        """Check if size of Atoms section is equal to expected with a ghost atom present."""
        fpath = os.path.join(__datadir__, "data/GAMESS/basicGAMESS-US2018/dvb_un_sp_ghost.out")
        data = cclib.io.ccread(fpath)
        writer = cclib.io.moldenwriter.MOLDEN(data)
        writer.ghost = "GH"
        # Check for the ghost atom, and make sure it has the right label.
        assert "GH" in writer._coords_from_ccdata(-1)[0]
        # Check size of Atoms section.
        assert len(writer._coords_from_ccdata(-1)) == data.natom

    def test_gto_section_size(self) -> None:
        """Check if size of GTO section is equal to expected."""
        fpath = os.path.join(__datadir__, "data/GAMESS/basicGAMESS-US2018/dvb_un_sp.out")
        data = cclib.io.ccread(fpath)
        writer = cclib.io.moldenwriter.MOLDEN(data)
        # Check size of GTO section.
        size_gto_ccdata = 0
        for atom in data.gbasis:
            size_gto_ccdata += 1
            for prims in atom:
                size_gto_ccdata += len(prims[1]) + 1
        # Filter blank lines.
        size_gto_writer = len(list(filter(None, writer._gto_from_ccdata())))
        assert size_gto_writer == size_gto_ccdata

    def test_mo_section_size(self) -> None:
        """Check if size of MO section is equal to expected."""
        fpath = os.path.join(__datadir__, "data/GAMESS/basicGAMESS-US2018/dvb_un_sp.out")
        data = cclib.io.ccread(fpath)
        writer = cclib.io.moldenwriter.MOLDEN(data)
        # Check size of MO section.
        size_mo_ccdata = 0
        extra = 4 if hasattr(data, "mosyms") else 3
        for i in range(data.mult):
            size_mo_ccdata += len(data.moenergies[i]) * (len(data.mocoeffs[i][0]) + extra)
        # Filter blank lines.
        (mosyms, moenergies, mooccs, mocoeffs) = (
            writer._syms_energies_occs_coeffs_from_ccdata_for_moldenwriter()
        )
        size_mo_writer = len(
            list(filter(None, writer._mo_from_ccdata(mosyms, moenergies, mooccs, mocoeffs)))
        )
        assert size_mo_writer == size_mo_ccdata

    def test_no_section_size(self) -> None:
        """Check if size of NO section is equal to expected."""
        fpath = os.path.join(__datadir__, "data/GAMESS/basicGAMESS-US2018/water_cis_dets.out")
        data = cclib.io.ccread(fpath)
        writer = cclib.io.moldenwriter.MOLDEN(data)
        # Check size of NO section.
        size_no_ccdata = 0
        extra = 4
        size_no_ccdata += len(data.nooccnos) * (len(data.nocoeffs[0]) + extra)
        # Filter blank lines.
        (nosyms, noenergies, nooccs, nocoeffs) = (
            writer._syms_energies_occs_coeffs_from_ccdata_for_moldenwriter()
        )
        size_no_writer = len(
            list(filter(None, writer._mo_from_ccdata(nosyms, noenergies, nooccs, nocoeffs)))
        )
        assert size_no_writer == size_no_ccdata

    def test_round_molden(self) -> None:
        """Check if Molden Style number rounding works as expected."""
        # If the 6th digit after dot is greater than 5, but is not 7,
        # round the number upto 6th place.
        # Else truncate at 6th digit after dot.
        assert round_molden(1) == 1
        assert round_molden(-1) == -1
        assert round_molden(0.999995789) == 0.999995
        assert round_molden(-0.999995789) == -0.999995
        assert round_molden(0.999996789) == 0.999997
        assert round_molden(-0.999997789) == -0.999997
        assert round_molden(0.999997789) == 0.999997
        assert round_molden(-0.999998789) == -0.999999
        assert round_molden(-0.999999999) == -1.0

    def test_molden_cclib_diff(self) -> None:
        """Check if file written by cclib matched file written by Molden."""
        filenames = ["dvb_un_sp", "C_bigbasis", "water_mp2", "dvb_ir"]
        for fn in filenames:
            fpath = os.path.join(__datadir__, f"data/GAMESS/basicGAMESS-US2018/{fn}.out")
            data = cclib.io.ccread(fpath)
            cclib_out = cclib.io.moldenwriter.MOLDEN(data).generate_repr()
            # Reformat cclib's output to remove extra spaces.
            cclib_out_formatted = MoldenReformatter(cclib_out).reformat()
            fpath = os.path.join(__testdir__, f"data/molden5.7_{fn}.molden")
            with open(fpath) as handle:
                molden_out = handle.read()
            # Reformat Molden's output to remove extra spaces,
            # and fix number formatting.
            molden_out_formatted = MoldenReformatter(molden_out).reformat()
            # Assert if reformatted files from both writers are same.
            assert molden_out_formatted == cclib_out_formatted
