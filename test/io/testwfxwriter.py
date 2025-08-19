# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for wfx writer."""

import os

import cclib
from cclib.io.filewriter import MissingAttributeError
from cclib.io.wfxwriter import _list_format, _section

import numpy as np
import pytest

__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")


class WFXTest:
    def test_missing_attribute_error(self) -> None:
        """Check if MissingAttributeError is raised as expected."""
        fpath = os.path.join(__datadir__, "data/GAMESS/basicGAMESS-US2017/C_bigbasis.out")
        required_attrs = ("atomcoords", "atomnos", "gbasis", "charge", "homos", "mult", "mocoeffs")
        for attr in required_attrs:
            data = cclib.io.ccread(fpath)
            delattr(data, attr)

            # WFX files cannot be written if required attrs are missing.
            with pytest.raises(MissingAttributeError):
                cclib.io.wfxwriter.WFXWriter(data)

    def test_no_of_prims(self) -> None:
        """Check if number of primitives are calculated correctly."""
        num_orb = {"s": 1, "p": 3, "d": 6, "f": 10, "g": 15, "h": 21}
        gamessdir = os.path.join(__datadir__, "data/GAMESS/basicGAMESS-US2017")
        filenames = ["C_bigbasis.out", "dvb_un_sp.out"]
        filepaths = [os.path.join(gamessdir, fn) for fn in filenames]

        for fpath in filepaths:
            data = cclib.io.ccread(fpath)
            wfx = cclib.io.wfxwriter.WFXWriter(data)

            no_prims_writer = wfx._no_of_prims()
            no_prims_ccdata = 0
            for atom in data.gbasis:
                for prims in atom:
                    no_prims_ccdata += num_orb[prims[0].lower()] * len(prims[1])

            assert no_prims_writer == no_prims_ccdata

    def test_mo_normalization(self) -> None:
        """Check if MO section is printed correctly."""
        fpath = os.path.join(__datadir__, "data/GAMESS/basicGAMESS-US2017/C_bigbasis.out")
        data = cclib.io.ccread(fpath)
        wfx = cclib.io.wfxwriter.WFXWriter(data)

        normalized_mocoeffs = wfx._normalized_mocoeffs()
        if len(data.homos) > 1:
            assert len(normalized_mocoeffs) == wfx._no_electrons()
        else:
            assert len(normalized_mocoeffs) == wfx._no_of_mos()
        assert len(normalized_mocoeffs[0]) == wfx._no_of_prims()

    def test_mo_normalization_dat(self) -> None:
        """Check if MOs are normalized as expected."""
        fpath = os.path.join(__datadir__, "data/GAMESS/basicGAMESS-US2017/dvb_sp.out")
        data = cclib.io.ccread(fpath)
        wfx = cclib.io.wfxwriter.WFXWriter(data)

        normalized_mocoeffs_wfx = wfx._normalized_mocoeffs()
        datfile = os.path.join(__datadir__, "data/GAMESS/basicGAMESS-US2017/dvb_sp.dat")

        with open(datfile) as file:
            content = iter(file.readlines())

        normalized_mocoeffs_dat = []

        line = next(content)
        while "TOP OF INPUT FILE FOR BADER" not in line:
            line = next(content)

        while "END OF INPUT FILE FOR BADER" not in line:
            if "MO" in line and "OCC NO" in line and "ORB. ENERGY" in line:
                line = next(content)
                mo = []
                while "MO" not in line and "END DATA" not in line:
                    mo.extend([float(x) for x in line.rsplit()])
                    line = next(content)
                normalized_mocoeffs_dat.append(mo)
            else:
                line = next(content)

        for mo1, mo2 in zip(normalized_mocoeffs_dat, normalized_mocoeffs_wfx):
            for coeff1, coeff2 in zip(mo1, mo2):
                assert -0.09 <= abs(coeff1) - abs(coeff2) <= 0.09

    def test_programs(self) -> None:
        """Check other programs against reference data."""
        ref_file = "data/GAMESS/basicGAMESS-US2017/dvb_sp.out"
        programs = {
            "ORCA": "data/ORCA/basicORCA5.0/dvb_sp.out",
            "NWChem": "data/NWChem/basicNWChem6.5/dvb_sp_hf.out",
            "Psi4": "data/Psi4/basicPsi4-1.3.1/dvb_sp_rhf.out",
            "GAMESS-UK": "data/GAMESS-UK/basicGAMESS-UK8.0/dvb_sp_hf.out",
            "Firefly": "data/GAMESS/basicFirefly8.0/dvb_sp.out",
        }
        fpath_ref = os.path.join(__datadir__, ref_file)
        data_ref = cclib.io.ccread(fpath_ref)
        wfx_ref = cclib.io.wfxwriter.WFXWriter(data_ref)
        norm_mat_ref = wfx_ref._norm_mat()[0]
        mos_ref = wfx_ref._no_of_mos()

        for name in programs:
            fpath_prog = os.path.join(__datadir__, programs[name])
            data_prog = cclib.io.ccread(fpath_prog)
            wfx_prog = cclib.io.wfxwriter.WFXWriter(data_prog)

            norm_mat_prog = wfx_prog._norm_mat()[0]
            mos_prog = wfx_prog._no_of_mos()

            # Check if normalization matrix are matching.
            assert round(abs(max(np.array(norm_mat_prog) - np.array(norm_mat_ref))), 5) == 0
            # Check if number of orbitals to be printed are as expected.
            assert mos_prog == mos_ref

    def test_section_printing(self) -> None:
        """Check if wfx section printing works as expected."""
        float_section = _section("Test Section", 123.456)
        expected = ["<Test Section>", " 123.456", "</Test Section>"]
        assert float_section == expected

        list_section = _section("Test Section", ["1", "2"])
        expected = ["<Test Section>", "1", "2", "</Test Section>"]
        assert list_section == expected

        str_section = _section("Test Section", "Test Data")
        expected = ["<Test Section>", " Test Data", "</Test Section>"]
        assert str_section == expected

    def test_list_format(self) -> None:
        """Check if list formatting works as expected."""
        odd_list = _list_format([1, 2, 3], 2, "%8.1E")
        odd_expected = [" 1.0E+00 2.0E+00", " 3.0E+00"]
        assert odd_list == odd_expected
        even_list = _list_format([1, 2, 3, 4], 2, "%8.1E")
        even_expected = [" 1.0E+00 2.0E+00", " 3.0E+00 4.0E+00"]
        assert even_list == even_expected
