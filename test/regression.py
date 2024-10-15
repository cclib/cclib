# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A regression framework for parsing and testing logfiles.

The intention here is to make it easy to add new datafiles as bugs
are fixed and to write specific tests in the form of test functions.

In short, the file called regressionfiles.yaml contains a list of regression
logfiles, which is compared to the files found on disk. All these files
should be parsed correctly, and if there is an appropriately named function
defined, that function will be used as a test.

There is also a mechanism for running unit tests on old logfiles, which
have been moved here from the cclib repository when newer versions
became available. We still want those logfiles to parse and test correctly,
although sometimes special modification will be needed.

To run all the tests, run `python -m pytest test/regression.py` from the
top level directory in the cclib repository.

Running all regressions can take anywhere from 10-20s to several minutes
depending in your hardware. To aid debugging, you can limit which regressions
to parse and test using the standard pytest `-k` flag. For example, to limit
the test to a specific parser:
    python -m pytest test/regression.py -k 'Gaussian'
You can also limit a run to a single output file using its 'normalized' name
inside the data directory, like so:
    python -m pytest test/regression.py -k 'Gaussian_Gaussian03_borane_opt_log'
"""

# ruff: noqa: E402, F401
import datetime
import logging
import sys
from pathlib import Path

from cclib.io import ccread, moldenwriter
from cclib.parser import DALTON, Gaussian, ccData

import numpy
import pytest
from packaging.version import Version
from packaging.version import parse as parse_version

# This assume that the cclib-data repository is located at a specific location
# within the cclib repository. It would be better to figure out a more natural
# way to import the relevant tests from cclib here.
#
# This is safer than sys.path.append, and isn't sys.path.insert(0, ...) so
# virtualenvs work properly. See https://stackoverflow.com/q/10095037.
base_dir = (Path(__file__) / ".." / "..").resolve()
__regression_dir__ = base_dir / "data" / "regression"
__filedir__ = base_dir / "test"
sys.path.insert(1, str(__filedir__))
# TODO There are many seemingly unused test imports here so that pytest can
# parameterize them with the desired files from
# cclib-data/regressionfiles.yaml.  If one is removed, the regression tests
# that depend on the imported test class will simply not run rather than fail.
# A future solution is to modify the pytest conftest.py so that all requested
# tests in regressionfiles.yaml are automatically collected without needing
# additional imports.
from .data.testBasis import (
    GaussianBigBasisTest,
    GenericBasisTest,
    GenericBigBasisTest,
    MolcasBigBasisTest,
    MolproBigBasisTest,
    Psi4BigBasisTest,
    QChemBigBasisTest,
)
from .data.testBOMD import GenericBOMDTest
from .data.testCC import GenericCCTest
from .data.testCI import GAMESSCISTest, GaussianCISTest, GenericCISTest, QChemCISTest
from .data.testCore import ADFCoreTest, GenericCoreTest
from .data.testGeoOpt import ADFGeoOptTest, GenericGeoOptTest, OrcaGeoOptTest, Psi4GeoOptTest
from .data.testMP import (
    GaussianMP2Test,
    GaussianMP3Test,
    GaussianMP4SDQTest,
    GaussianMP4SDTQTest,
    GenericMP2Test,
    GenericMP3Test,
    GenericMP4SDQTest,
    GenericMP4SDTQTest,
    GenericMP5Test,
    QChemMP4SDQTest,
    QChemMP4SDTQTest,
)
from .data.testPolar import GenericPolarTest, ReferencePolarTest
from .data.testScan import GaussianRelaxedScanTest, GenericRelaxedScanTest
from .data.testSP import (
    ADFSPTest,
    GaussianSPTest,
    GenericHFSPTest,
    GenericSPTest,
    JaguarSPTest,
    MolcasSPTest,
    OrcaSPTest,
    PsiHFSPTest,
    PsiSPTest,
)
from .data.testSPun import GaussianSPunTest, GenericROSPTest, GenericSPunTest, JaguarSPunTest
from .data.testTD import (
    DALTONTDTest,
    GAMESSUSTDDFTTest,
    GaussianTDDFTTest,
    GenericTDDFTtrpTest,
    GenericTDTest,
    OrcaROCISTest,
    OrcaTDDFTTest,
    QChemTDDFTTest,
)
from .data.testTDun import GenericTDunTest
from .data.testvib import (
    ADFIRTest,
    FireflyIRTest,
    GamessIRTest,
    GaussianIRTest,
    GaussianRamanTest,
    GenericIRimgTest,
    GenericIRTest,
    GenericRamanTest,
    JaguarIRTest,
    OrcaIRTest,
    OrcaRamanTest,
    Psi4HFIRTest,
    QChemRamanTest,
)

# The following regression test functions were manually written, because they
# contain custom checks that were determined on a per-file basis. Care needs to be taken
# that the function name corresponds to the path of the logfile, with some characters
# changed according to normalisefilename().

# ADF #


def testADF_ADF2004_01_Fe_ox3_final_out(logfile) -> None:
    """Make sure HOMOS are correct."""
    assert logfile.data.homos[0] == 59 and logfile.data.homos[1] == 54

    assert logfile.data.metadata["legacy_package_version"] == "2004.01"
    assert logfile.data.metadata["package_version"] == "2004.01+200410211341"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testADF_ADF2013_01_dvb_gopt_b_unconverged_adfout(logfile) -> None:
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert hasattr(logfile.data, "optdone") and not logfile.data.optdone

    assert logfile.data.metadata["legacy_package_version"] == "2013.01"
    assert logfile.data.metadata["package_version"] == "2013.01+201309012319"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testADF_ADF2013_01_stopiter_dvb_sp_adfout(logfile):
    """This logfile has not SCF test lines so we have no way to check what happens."""
    # This is what we would have checked:
    # len(logfile.data.scfvalues[0]) == 10
    assert not hasattr(logfile.data, "scfvalues")

    assert logfile.data.metadata["package_version"] == "2013.01+201309012319"


def testADF_ADF2013_01_stopiter_dvb_sp_b_adfout(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    # Why is this not 3?
    assert len(logfile.data.scfvalues[0]) == 2

    assert logfile.data.metadata["package_version"] == "2013.01+201309012319"


def testADF_ADF2013_01_stopiter_dvb_sp_c_adfout(logfile):
    """This logfile has not SCF test lines so we have no way to check what happens."""
    # This is what we would have checked:
    # len(logfile.data.scfvalues[0]) == 6
    assert not hasattr(logfile.data, "scfvalues")

    assert logfile.data.metadata["package_version"] == "2013.01+201309012319"


def testADF_ADF2013_01_stopiter_dvb_sp_d_adfout(logfile):
    """This logfile has not SCF test lines so we have no way to check what happens."""
    # This is what we would have checked:
    # len(logfile.data.scfvalues[0]) == 7
    assert not hasattr(logfile.data, "scfvalues")

    assert logfile.data.metadata["package_version"] == "2013.01+201309012319"


def testADF_ADF2013_01_stopiter_dvb_un_sp_adfout(logfile):
    """This logfile has not SCF test lines so we have no way to check what happens."""
    # This is what we would have checked:
    # len(logfile.data.scfvalues[0]) == 7
    assert not hasattr(logfile.data, "scfvalues")

    assert logfile.data.metadata["package_version"] == "2013.01+201309012319"


def testADF_ADF2013_01_stopiter_dvb_un_sp_c_adfout(logfile):
    """This logfile has not SCF test lines so we have no way to check what happens."""
    # This is what we would have checked:
    # len(logfile.data.scfvalues[0]) == 10
    assert not hasattr(logfile.data, "scfvalues")

    assert logfile.data.metadata["package_version"] == "2013.01+201309012319"


def testADF_ADF2013_01_stopiter_MoOCl4_sp_adfout(logfile):
    """This logfile has not SCF test lines so we have no way to check what happens."""
    # This is what we would have checked:
    # len(logfile.data.scfvalues[0]) == 11
    assert not hasattr(logfile.data, "scfvalues")

    assert logfile.data.metadata["package_version"] == "2013.01+201309012319"


def testADF_ADF2014_01_DMO_ORD_orig_out(logfile):
    """In lieu of a unit test, make sure the polarizability (and
    potentially later the optical rotation) is properly parsed.
    """
    assert hasattr(logfile.data, "polarizabilities")
    assert len(logfile.data.polarizabilities) == 1
    assert logfile.data.polarizabilities[0].shape == (3, 3)

    # isotropic polarizability
    isotropic_calc = numpy.average(numpy.diag(logfile.data.polarizabilities[0]))
    isotropic_ref = 51.3359
    assert abs(isotropic_calc - isotropic_ref) < 1.0e-4

    assert logfile.data.metadata["legacy_package_version"] == "2014"
    assert logfile.data.metadata["package_version"] == "2014dev42059"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)
    assert logfile.data.metadata["package_version_date"] == "2014-06-11"
    assert logfile.data.metadata["package_version_description"] == "development version"


def testADF_ADF2016_166_tddft_0_31_new_out(logfile):
    """This file led to StopIteration (#430)."""
    assert logfile.data.metadata["legacy_package_version"] == "2016"
    assert logfile.data.metadata["package_version"] == "2016dev53619"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)
    assert logfile.data.metadata["package_version_date"] == "2016-07-21"
    assert "package_version_description" not in logfile.data.metadata


def testADF_ADF2016_fa2_adf_out(logfile):
    """This logfile, without symmetry, should get atombasis parsed."""
    assert hasattr(logfile.data, "atombasis")
    assert [b for ab in logfile.data.atombasis for b in ab] == list(range(logfile.data.nbasis))

    assert logfile.data.metadata["legacy_package_version"] == "2016"
    assert logfile.data.metadata["package_version"] == "2016dev50467"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)
    assert logfile.data.metadata["package_version_date"] == "2016-02-17"
    assert logfile.data.metadata["package_version_description"] == "branches/AndrewAtkins/ADF-Shar"


# DALTON #


def testDALTON_DALTON_2013_dvb_td_normalprint_out(logfile):
    r"""This original unit test prints a DFT-specific version of the excitation
    eigenvectors, which we do not parse.

    Here is an example of the general output (requiring `**RESPONSE/.PRINT 4`
    for older versions of DALTON), followed by "PBHT MO Overlap Diagnostic"
    which only appears for DFT calculations. Note that the reason we cannot
    parse this for etsyms is it doesn't contain the necessary
    coefficient. "K_IA" and "(r s) operator", which is $\kappa_{rs}$, the
    coefficient for excitation from the r -> s MO in the response vector, is
    not what most programs print; it is "(r s) scaled", which is $\kappa_{rs}
    * \sqrt{S_{rr} - S_{ss}}$. Because this isn't available from the PBHT
    output, we cannot parse it.

         Eigenvector for state no.  1

             Response orbital operator symmetry = 1
             (only scaled elements abs greater than   10.00 % of max abs value)

              Index(r,s)      r      s        (r s) operator      (s r) operator      (r s) scaled        (s r) scaled
              ----------    -----  -----      --------------      --------------      --------------      --------------
                 154        27(2)  28(2)        0.5645327267        0.0077924161        0.7983698385        0.0110201405
                 311        58(4)  59(4)       -0.4223079545        0.0137981027       -0.5972336367        0.0195134639

        ...

                                    PBHT MO Overlap Diagnostic
                                    --------------------------

              I    A    K_IA      K_AI   <|I|*|A|> <I^2*A^2>    Weight   Contrib

             27   28  0.564533  0.007792  0.790146  0.644560  0.309960  0.244913
             58   59 -0.422308  0.013798  0.784974  0.651925  0.190188  0.149293

    In the future, if `aooverlaps` and `mocoeffs` are available, it may be
    possible to calculate the necessary scaled coefficients for `etsecs`.
    """
    assert hasattr(logfile.data, "etenergies")
    assert not hasattr(logfile.data, "etsecs")
    assert hasattr(logfile.data, "etsyms")
    assert hasattr(logfile.data, "etoscs")

    assert logfile.data.metadata["legacy_package_version"] == "2013.4"
    assert (
        logfile.data.metadata["package_version"]
        == "2013.4+7abef2ada27562fe5e02849d6caeaa67c961732f"
    )
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testDALTON_DALTON_2015_dalton_atombasis_out(logfile):
    """This logfile didn't parse due to the absence of a line in the basis
    set section.
    """
    assert hasattr(logfile.data, "nbasis")
    assert logfile.data.nbasis == 37
    assert hasattr(logfile.data, "atombasis")

    assert logfile.data.metadata["legacy_package_version"] == "2015.0"
    assert (
        logfile.data.metadata["package_version"]
        == "2015.0+d34efb170c481236ad60c789dea90a4c857c6bab"
    )
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testDALTON_DALTON_2015_dalton_intgrl_out(logfile):
    """This logfile didn't parse due to the absence of a line in the basis
    set section.
    """
    assert hasattr(logfile.data, "nbasis")
    assert logfile.data.nbasis == 4
    assert hasattr(logfile.data, "atombasis")

    assert (
        logfile.data.metadata["package_version"]
        == "2015.0+d34efb170c481236ad60c789dea90a4c857c6bab"
    )


def testDALTON_DALTON_2015_dvb_td_normalprint_out(logfile):
    """This original unit test prints a DFT-specific version of the excitation
    eigenvectors, which we do not parse.
    """
    assert hasattr(logfile.data, "etenergies")
    assert not hasattr(logfile.data, "etsecs")
    assert hasattr(logfile.data, "etsyms")
    assert hasattr(logfile.data, "etoscs")

    assert (
        logfile.data.metadata["package_version"]
        == "2015.0+d34efb170c481236ad60c789dea90a4c857c6bab"
    )


def testDALTON_DALTON_2015_stopiter_dalton_dft_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert len(logfile.data.scfvalues[0]) == 8

    assert (
        logfile.data.metadata["package_version"]
        == "2015.0+d34efb170c481236ad60c789dea90a4c857c6bab"
    )


def testDALTON_DALTON_2015_stopiter_dalton_hf_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert len(logfile.data.scfvalues[0]) == 5

    assert (
        logfile.data.metadata["package_version"]
        == "2015.0+d34efb170c481236ad60c789dea90a4c857c6bab"
    )


def testDALTON_DALTON_2016_huge_neg_polar_freq_out(logfile):
    """This is an example of a multiple frequency-dependent polarizability
    calculation.
    """
    assert hasattr(logfile.data, "polarizabilities")
    assert len(logfile.data.polarizabilities) == 3
    assert abs(logfile.data.polarizabilities[2][0, 0] - 183.6308) < 1.0e-5

    assert logfile.data.metadata["legacy_package_version"] == "2016.2"
    assert (
        logfile.data.metadata["package_version"]
        == "2016.2+7db4647eac203e51aae7da3cbc289f55146b30e9"
    )
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testDALTON_DALTON_2016_huge_neg_polar_stat_out(logfile):
    """This logfile didn't parse due to lack of spacing between
    polarizability tensor elements.
    """
    assert hasattr(logfile.data, "polarizabilities")
    assert len(logfile.data.polarizabilities) == 1
    assert abs(logfile.data.polarizabilities[0][1, 1] + 7220.150408) < 1.0e-7

    assert (
        logfile.data.metadata["package_version"]
        == "2016.2+7db4647eac203e51aae7da3cbc289f55146b30e9"
    )


def testDALTON_DALTON_2016_Trp_polar_response_diplnx_out(logfile):
    """Check that only the xx component of polarizability is defined and
    all others are NaN even after parsing a previous file with full tensor.
    """
    # TODO replace with looking at file location from fixture?
    full_tens_path = __regression_dir__ / "DALTON" / "DALTON-2015" / "Trp_polar_response.out"
    DALTON(full_tens_path, loglevel=logging.ERROR).parse()
    assert hasattr(logfile.data, "polarizabilities")
    assert abs(logfile.data.polarizabilities[0][0, 0] - 95.11540019) < 1.0e-8
    assert numpy.count_nonzero(numpy.isnan(logfile.data.polarizabilities)) == 8

    assert (
        logfile.data.metadata["package_version"]
        == "2016.2+7db4647eac203e51aae7da3cbc289f55146b30e9"
    )


def testDALTON_DALTON_2018_dft_properties_nosym_H2O_cc_pVDZ_out(logfile):
    """The "simple" version string in newer development versions of DALTON wasn't
    being parsed properly.

    This file is in DALTON-2018, rather than DALTON-2019, because 2018.0 was
    just released.
    """
    assert logfile.data.metadata["legacy_package_version"] == "2019.alpha"
    assert logfile.data.metadata["package_version"] == "2019.alpha"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testDALTON_DALTON_2018_tdhf_2000_out(logfile):
    """Ensure etsecs are being parsed from a TDHF calculation without symmetry and
    a big print level.
    """
    assert hasattr(logfile.data, "etsecs")
    for attr in ("etenergies", "etsecs", "etsyms", "etoscs"):
        assert len(getattr(logfile.data, attr)) == 9
    assert logfile.data.etsecs[0][0] == [(1, 0), (2, 0), -0.9733558768]

    assert logfile.data.metadata["legacy_package_version"] == "2019.alpha"
    assert (
        logfile.data.metadata["package_version"]
        == "2019.alpha+25947a3d842ee2ebb42bff87a4dd64adbbd3ec5b"
    )
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testDALTON_DALTON_2018_tdhf_2000_sym_out(logfile):
    """Ensure etsecs are being parsed from a TDHF calculation with symmetry and a
    big print level.
    """
    assert hasattr(logfile.data, "etsecs")
    for attr in ("etenergies", "etsecs", "etsyms", "etoscs"):
        assert len(getattr(logfile.data, attr)) == 3
    assert logfile.data.etsecs[0][0] == [(1, 0), (2, 0), 0.9733562358]

    assert (
        logfile.data.metadata["package_version"]
        == "2019.alpha+25947a3d842ee2ebb42bff87a4dd64adbbd3ec5b"
    )


def testDALTON_DALTON_2018_tdhf_normal_out(logfile):
    """Ensure etsecs are being parsed from a TDHF calculation without symmetry and
    a normal print level.
    """
    assert hasattr(logfile.data, "etsecs")
    for attr in ("etenergies", "etsecs", "etsyms", "etoscs"):
        assert len(getattr(logfile.data, attr)) == 9
    assert logfile.data.etsecs[0][0] == [(1, 0), (2, 0), -0.9733558768]

    assert (
        logfile.data.metadata["package_version"]
        == "2019.alpha+25947a3d842ee2ebb42bff87a4dd64adbbd3ec5b"
    )


def testDALTON_DALTON_2018_tdhf_normal_sym_out(logfile):
    """Ensure etsecs are being parsed from a TDHF calculation with symmetry and a
    normal print level.
    """
    assert hasattr(logfile.data, "etsecs")
    for attr in ("etenergies", "etsecs", "etsyms", "etoscs"):
        assert len(getattr(logfile.data, attr)) == 3
    assert logfile.data.etsecs[0][0] == [(1, 0), (2, 0), 0.9733562358]

    assert (
        logfile.data.metadata["package_version"]
        == "2019.alpha+25947a3d842ee2ebb42bff87a4dd64adbbd3ec5b"
    )


def testDALTON_DALTON_2018_tdpbe_2000_out(logfile):
    """Ensure etsecs are being parsed from a TDDFT calculation without symmetry
    and a big print level.
    """
    assert hasattr(logfile.data, "etsecs")
    for attr in ("etenergies", "etsecs", "etsyms", "etoscs"):
        assert len(getattr(logfile.data, attr)) == 9
    assert logfile.data.etsecs[0][0] == [(1, 0), (2, 0), 0.9992665559]

    assert (
        logfile.data.metadata["package_version"]
        == "2019.alpha+25947a3d842ee2ebb42bff87a4dd64adbbd3ec5b"
    )


def testDALTON_DALTON_2018_tdpbe_2000_sym_out(logfile):
    """Ensure etsecs are being parsed from a TDDFT calculation with symmetry and a
    big print level.
    """
    assert hasattr(logfile.data, "etsecs")
    for attr in ("etenergies", "etsecs", "etsyms", "etoscs"):
        assert len(getattr(logfile.data, attr)) == 3
    assert logfile.data.etsecs[0][0] == [(1, 0), (2, 0), 0.9992672154]

    assert (
        logfile.data.metadata["package_version"]
        == "2019.alpha+25947a3d842ee2ebb42bff87a4dd64adbbd3ec5b"
    )


def testDALTON_DALTON_2018_tdpbe_normal_out(logfile):
    """Ensure etsecs are being parsed from a TDDFT calculation without symmetry
    and a normal print level.
    """
    assert hasattr(logfile.data, "etsecs")
    for attr in ("etenergies", "etsecs", "etsyms", "etoscs"):
        assert len(getattr(logfile.data, attr)) == 9
    assert logfile.data.etsecs[0][0] == [(1, 0), (2, 0), 0.9992665559]

    assert (
        logfile.data.metadata["package_version"]
        == "2019.alpha+25947a3d842ee2ebb42bff87a4dd64adbbd3ec5b"
    )


def testDALTON_DALTON_2018_tdpbe_normal_sym_out(logfile):
    """Ensure etsecs are being parsed from a TDDFT calculation with symmetry and a
    normal print level.
    """
    assert hasattr(logfile.data, "etsecs")
    for attr in ("etenergies", "etsecs", "etsyms", "etoscs"):
        assert len(getattr(logfile.data, attr)) == 3
    assert logfile.data.etsecs[0][0] == [(1, 0), (2, 0), 0.9992672154]

    assert (
        logfile.data.metadata["package_version"]
        == "2019.alpha+25947a3d842ee2ebb42bff87a4dd64adbbd3ec5b"
    )


# Formatted checkpoint #


def testFChk_Gaussian03_dvb_gopt_unconverged_fchk(logfile):
    metadata = logfile.data.metadata
    assert metadata["package"] == "FChk[Gaussian]"
    # Impossible to determined based upon current parsed data
    assert "success" not in metadata


def testFChk_Gaussian16_dvb_gopt_unconverged_fchk(logfile):
    metadata = logfile.data.metadata
    assert metadata["package"] == "FChk[Gaussian]"
    # >= g16 has "Job Status"
    assert "success" in metadata
    assert not metadata["success"]


def testFChk_QChem5_3_dvb_gopt_unconverged_in_fchk(logfile):
    metadata = logfile.data.metadata
    assert metadata["package"] == "FChk[QChem]"
    # Determined because mocoeffs are missing
    assert "success" in metadata
    assert not metadata["success"]


def testFChk_QChem5_3_dvb_sp_unconverged_in_fchk(logfile):
    metadata = logfile.data.metadata
    assert metadata["package"] == "FChk[QChem]"
    # Determined because mocoeffs are missing
    assert "success" in metadata
    assert not metadata["success"]


# Firefly #


def testGAMESS_Firefly8_0_dvb_gopt_a_unconverged_out(logfile):
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert hasattr(logfile.data, "optdone") and not logfile.data.optdone

    assert logfile.data.metadata["legacy_package_version"] == "8.0.1"
    assert logfile.data.metadata["package_version"] == "8.0.1+8540"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGAMESS_Firefly8_0_h2o_log(logfile):
    """Check that molecular orbitals are parsed correctly (cclib/cclib#208)."""
    assert logfile.data.mocoeffs[0][0][0] == -0.994216

    assert logfile.data.metadata["legacy_package_version"] == "8.0.0"
    assert logfile.data.metadata["package_version"] == "8.0.0+7651"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGAMESS_Firefly8_0_stopiter_firefly_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert len(logfile.data.scfvalues[0]) == 6

    assert logfile.data.metadata["package_version"] == "8.0.1+8540"


def testGAMESS_Firefly8_1_benzene_am1_log(logfile):
    """Molecular orbitals were not parsed (cclib/cclib#228)."""
    assert hasattr(logfile.data, "mocoeffs")

    assert logfile.data.metadata["legacy_package_version"] == "8.1.0"
    assert logfile.data.metadata["package_version"] == "8.1.0+9035"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGAMESS_Firefly8_1_naphtalene_t_0_out(logfile):
    """Molecular orbitals were not parsed (cclib/cclib#228)."""
    assert hasattr(logfile.data, "mocoeffs")

    assert logfile.data.metadata["legacy_package_version"] == "8.1.1"
    assert logfile.data.metadata["package_version"] == "8.1.1+9295"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGAMESS_Firefly8_1_naphtalene_t_0_SP_out(logfile):
    """Molecular orbitals were not parsed (cclib/cclib#228)."""
    assert hasattr(logfile.data, "mocoeffs")

    assert logfile.data.metadata["package_version"] == "8.1.1+9295"


# GAMESS #


def testGAMESS_GAMESS_US2008_N2_UMP2_out(logfile):
    """Check that the new format for GAMESS MP2 is parsed."""
    assert hasattr(logfile.data, "mpenergies")
    assert len(logfile.data.mpenergies) == 1
    assert abs(logfile.data.mpenergies[0] - -109.3647999161) < 1.0e-10

    assert logfile.data.metadata["legacy_package_version"] == "2008R1"
    assert logfile.data.metadata["package_version"] == "2008.r1"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGAMESS_GAMESS_US2008_N2_ROMP2_out(logfile):
    """Check that the new format for GAMESS MP2 is parsed."""
    assert hasattr(logfile.data, "mpenergies")
    assert len(logfile.data.mpenergies) == 1
    assert abs(logfile.data.mpenergies[0] - -109.3647999184) < 1.0e-10

    assert logfile.data.metadata["package_version"] == "2008.r1"


def testGAMESS_GAMESS_US2009_open_shell_ccsd_test_log(logfile):
    """Parse ccenergies from open shell CCSD calculations."""
    assert hasattr(logfile.data, "ccenergies")
    assert len(logfile.data.ccenergies) == 1
    assert abs(logfile.data.ccenergies[0] - -128.6777922565) < 1.0e-10

    assert logfile.data.metadata["legacy_package_version"] == "2009R3"
    assert logfile.data.metadata["package_version"] == "2009.r3"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGAMESS_GAMESS_US2009_paulo_h2o_mp2_out(logfile):
    """Check that the new format for GAMESS MP2 is parsed."""
    assert hasattr(logfile.data, "mpenergies")
    assert len(logfile.data.mpenergies) == 1
    assert abs(logfile.data.mpenergies[0] - -76.1492222841) < 1.0e-10

    assert logfile.data.metadata["package_version"] == "2009.r3"


def testGAMESS_GAMESS_US2012_dvb_gopt_a_unconverged_out(logfile):
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert hasattr(logfile.data, "optdone") and not logfile.data.optdone

    assert logfile.data.metadata["legacy_package_version"] == "2012R2"
    assert logfile.data.metadata["package_version"] == "2012.r2"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGAMESS_GAMESS_US2012_stopiter_gamess_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert len(logfile.data.scfvalues[0]) == 10

    assert logfile.data.metadata["package_version"] == "2012.r1"


def testGAMESS_GAMESS_US2013_N_UHF_out(logfile):
    """An UHF job that has an LZ value analysis between the alpha and beta orbitals."""
    assert len(logfile.data.moenergies) == 2

    assert logfile.data.metadata["legacy_package_version"] == "2013R1"
    assert logfile.data.metadata["package_version"] == "2013.r1"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGAMESS_GAMESS_US2014_CdtetraM1B3LYP_log(logfile):
    """This logfile had coefficients for only 80 molecular orbitals."""
    assert len(logfile.data.mocoeffs) == 2
    assert numpy.count_nonzero(logfile.data.mocoeffs[0][79 - 1 :, :]) == 258
    assert numpy.count_nonzero(logfile.data.mocoeffs[0][80 - 1 : 0 :]) == 0
    assert logfile.data.mocoeffs[0].all() == logfile.data.mocoeffs[1].all()

    assert logfile.data.metadata["legacy_package_version"] == "2014R1"
    assert logfile.data.metadata["package_version"] == "2014.r1"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGAMESS_GAMESS_US2018_exam45_log(logfile):
    """This logfile has EOM-CC electronic transitions (not currently supported)."""
    assert not hasattr(logfile.data, "etenergies")

    assert logfile.data.metadata["legacy_package_version"] == "2018R2"
    assert logfile.data.metadata["package_version"] == "2018.r2"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGAMESS_GAMESS_US2018_exam46_log(logfile):
    """
    This logfile has >100 scf iterations, which used to cause
    a parsing error.
    """
    assert len(logfile.data.scfvalues[0]) == 113
    assert logfile.data.metadata["legacy_package_version"] == "2018R3"
    assert logfile.data.metadata["package_version"] == "2018.r3"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGAMESS_WinGAMESS_dvb_td_trplet_2007_03_24_r1_out(logfile):
    """Do some basic checks for this old unit test that was failing.

    The unit tests are not run automatically on this old unit logfile,
    because we know the output has etsecs whose sum is way off.
    So, perform a subset of the basic assertions for GenericTDTesttrp.
    """
    number = 5
    assert len(logfile.data.etenergies) == number
    idx_lambdamax = [i for i, x in enumerate(logfile.data.etoscs) if x == max(logfile.data.etoscs)][
        0
    ]
    assert (
        abs(logfile.data.etenergies[idx_lambdamax] - (-381.9320539243 - -382.0432999970)) < 1.0e-5
    )
    assert len(logfile.data.etoscs) == number
    assert abs(max(logfile.data.etoscs) - 0.0) < 0.01
    assert len(logfile.data.etsecs) == number

    assert logfile.data.metadata["legacy_package_version"] == "2007R1"
    assert logfile.data.metadata["package_version"] == "2007.r1"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testnoparseGAMESS_WinGAMESS_H2O_def2SVPD_triplet_2019_06_30_R1_out(filename):
    """Check if the molden writer can handle an unrestricted case"""
    data = ccread(__filedir__ / filename)
    writer = moldenwriter.MOLDEN(data)
    (mosyms, moenergies, mooccs, mocoeffs) = (
        writer._syms_energies_occs_coeffs_from_ccdata_for_moldenwriter()
    )
    molden_lines = writer._mo_from_ccdata(mosyms, moenergies, mooccs, mocoeffs)
    # Check size of Atoms section.
    assert len(molden_lines) == (data.nbasis + 4) * (data.nmo * 2)
    # check docc orbital
    beta_idx = (data.nbasis + 4) * data.nmo
    assert "Beta" in molden_lines[beta_idx + 2]
    assert "Occup=   1.000000" in molden_lines[beta_idx + 3]
    assert "0.989063" in molden_lines[beta_idx + 4]


# GAMESS-UK #


def testGAMESS_UK_GAMESS_UK8_0_dvb_gopt_hf_unconverged_out(logfile):
    assert hasattr(logfile.data, "optdone") and not logfile.data.optdone

    assert logfile.data.metadata["legacy_package_version"] == "8.0"
    assert logfile.data.metadata["package_version"] == "8.0+6248"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGAMESS_UK_GAMESS_UK8_0_stopiter_gamessuk_dft_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert len(logfile.data.scfvalues[0]) == 7

    assert logfile.data.metadata["package_version"] == "8.0+6248"


def testGAMESS_UK_GAMESS_UK8_0_stopiter_gamessuk_hf_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert len(logfile.data.scfvalues[0]) == 5

    assert logfile.data.metadata["package_version"] == "8.0+6248"


# Gaussian #


def testGaussian_Gaussian98_C_bigmult_log(logfile):
    """
    This file failed first becuase it had a double digit multiplicity.
    Then it failed because it had no alpha virtual orbitals.
    """
    assert logfile.data.charge == -3
    assert logfile.data.mult == 10
    assert logfile.data.homos[0] == 8
    assert logfile.data.homos[1] == -1  # No occupied beta orbitals

    assert logfile.data.metadata["legacy_package_version"] == "98revisionA.11.3"
    assert logfile.data.metadata["package_version"] == "1998+A.11.3"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGaussian_Gaussian98_NIST_CCCBDB_1himidaz_m21b0_out(logfile):
    """A G3 computation is a sequence of jobs."""

    # All steps deal with the same molecule, so we extract the coordinates
    # from all steps.
    assert len(logfile.data.atomcoords) == 10

    # Different G3 steps do perturbation to different orders, and so
    # we expect only the last MP2 energy to be extracted.
    assert len(logfile.data.mpenergies) == 1

    assert logfile.data.metadata["legacy_package_version"] == "98revisionA.7"
    assert logfile.data.metadata["package_version"] == "1998+A.7"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGaussian_Gaussian98_NIST_CCCBDB_1himidaz_m23b6_out(logfile):
    """A job that was killed before it ended, should have several basic attributes parsed."""
    assert hasattr(logfile.data, "charge")
    assert hasattr(logfile.data, "metadata")
    assert hasattr(logfile.data, "mult")

    assert logfile.data.metadata["package_version"] == "1998+A.7"


def testGaussian_Gaussian98_test_Cu2_log(logfile):
    """An example of the number of basis set function changing."""
    assert logfile.data.nbasis == 38

    assert (
        logfile.data.metadata["cpu_time"]
        == logfile.data.metadata["wall_time"]
        == [datetime.timedelta(seconds=25, microseconds=800000)]
    )
    assert logfile.data.metadata["legacy_package_version"] == "98revisionA.11.4"
    assert logfile.data.metadata["package_version"] == "1998+A.11.4"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGaussian_Gaussian98_test_H2_log(logfile):
    """
    The atomic charges from a natural population analysis were
    not parsed correctly, and they should be zero for dihydrogen.
    """
    assert logfile.data.atomcharges["natural"][0] == 0.0
    assert logfile.data.atomcharges["natural"][1] == 0.0

    assert logfile.data.metadata["package_version"] == "1998+A.11.4"


def testGaussian_Gaussian98_water_zmatrix_nosym_log(logfile):
    """This file is missing natom.

    This file had no atomcoords as it did not contain either an
    "Input orientation" or "Standard orientation section".
    As a result it failed to parse. Fixed in r400.
    """
    assert len(logfile.data.atomcoords) == 1
    assert logfile.data.natom == 3

    assert logfile.data.metadata["package_version"] == "1998+A.11.3"


def testGaussian_Gaussian03_AM1_SP_out(logfile):
    """Previously, caused scfvalue parsing to fail."""
    assert len(logfile.data.scfvalues[0]) == 13

    assert logfile.data.metadata["legacy_package_version"] == "03revisionE.01"
    assert logfile.data.metadata["package_version"] == "2003+E.01"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGaussian_Gaussian03_OPT_td_out(logfile):
    """Working fine - adding to ensure that CD is parsed correctly."""
    assert len(logfile.data.etrotats) == 10
    assert logfile.data.etrotats[0] == -0.4568

    assert logfile.data.metadata["package_version"] == "2003+B.05"


def testGaussian_Gaussian03_anthracene_log(logfile):
    """This file exposed a bug in extracting the vibsyms."""
    assert len(logfile.data.vibsyms) == len(logfile.data.vibfreqs)

    assert logfile.data.metadata["legacy_package_version"] == "03revisionC.02"
    assert logfile.data.metadata["package_version"] == "2003+C.02"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGaussian_Gaussian03_borane_opt_log(logfile):
    """An example of changing molecular orbital count."""
    assert logfile.data.optstatus[-1] == logfile.data.OPT_DONE
    assert logfile.data.nmo == 609

    assert logfile.data.metadata["package_version"] == "2003+E.01"


def testGaussian_Gaussian03_chn1_log(logfile):
    """
    This file failed to parse, due to the use of 'pop=regular'.
    We have decided that mocoeffs should not be defined for such calculations.
    """
    assert not hasattr(logfile.data, "mocoeffs")

    assert logfile.data.metadata["legacy_package_version"] == "03revisionB.04"
    assert logfile.data.metadata["package_version"] == "2003+B.04"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGaussian_Gaussian03_cyclopropenyl_rhf_g03_cut_log(logfile):
    """
    Not using symmetry at all (option nosymm) means standard orientation
    is not printed. In this case inputcoords are copied by the parser,
    which up till now stored the last coordinates.
    """
    assert len(logfile.data.atomcoords) == len(logfile.data.geovalues)

    assert logfile.data.metadata["package_version"] == "2003+C.02"


def testGaussian_Gaussian03_DCV4T_C60_log(logfile):
    """This is a test for a very large Gaussian file with > 99 atoms.

    The log file is too big, so we are just including the start.
    Previously, parsing failed in the pseudopotential section.
    """
    assert len(logfile.data.coreelectrons) == 102
    assert logfile.data.coreelectrons[101] == 2

    assert logfile.data.metadata["legacy_package_version"] == "03revisionD.02"
    assert logfile.data.metadata["package_version"] == "2003+D.02"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGaussian_Gaussian03_dvb_gopt_symmfollow_log(logfile):
    """Non-standard treatment of symmetry.

    In this case the Standard orientation is also printed non-standard,
    which caused only the first coordinates to be read previously.
    """
    assert len(logfile.data.atomcoords) == len(logfile.data.geovalues)

    # This calc only uses one CPU, so wall_time == cpu_time.
    assert (
        logfile.data.metadata["cpu_time"]
        == logfile.data.metadata["wall_time"]
        == [datetime.timedelta(seconds=99)]
    )
    assert logfile.data.metadata["legacy_package_version"] == "03revisionC.01"
    assert logfile.data.metadata["package_version"] == "2003+C.01"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGaussian_Gaussian03_mendes_out(logfile):
    """Previously, failed to extract coreelectrons."""
    centers = [9, 10, 11, 27]
    for i, x in enumerate(logfile.data.coreelectrons):
        if i in centers:
            assert x == 10
        else:
            assert x == 0

    assert logfile.data.metadata["package_version"] == "2003+C.02"


def testGaussian_Gaussian03_Mo4OSibdt2_opt_log(logfile):
    """
    This file had no atomcoords as it did not contain any
    "Input orientation" sections, only "Standard orientation".
    """
    assert logfile.data.optstatus[-1] == logfile.data.OPT_DONE
    assert hasattr(logfile.data, "atomcoords")

    assert logfile.data.metadata["package_version"] == "2003+C.02"


def testGaussian_Gaussian03_orbgs_log(logfile):
    """Check that the pseudopotential is being parsed correctly."""
    assert hasattr(logfile.data, "coreelectrons"), "Missing coreelectrons"
    assert logfile.data.coreelectrons[0] == 28
    assert logfile.data.coreelectrons[15] == 10
    assert logfile.data.coreelectrons[20] == 10
    assert logfile.data.coreelectrons[23] == 10

    assert logfile.data.metadata["package_version"] == "2003+C.02"


def testGaussian_Gaussian09_100_g09(logfile):
    """Check that the final system is the one parsed (cclib/cclib#243)."""
    assert logfile.data.natom == 54
    assert logfile.data.homos == [104]

    assert logfile.data.metadata["legacy_package_version"] == "09revisionB.01"
    assert logfile.data.metadata["package_version"] == "2009+B.01"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGaussian_Gaussian09_25DMF_HRANH_log(logfile):
    """Check that the anharmonicities are being parsed correctly."""
    assert hasattr(logfile.data, "vibanharms"), "Missing vibanharms"
    anharms = logfile.data.vibanharms
    N = len(logfile.data.vibfreqs)
    assert 39 == N == anharms.shape[0] == anharms.shape[1]
    assert abs(anharms[0][0] + 43.341) < 0.01
    assert abs(anharms[N - 1][N - 1] + 36.481) < 0.01

    assert logfile.data.metadata["package_version"] == "2009+B.01"


def testGaussian_Gaussian09_2D_PES_all_converged_log(logfile):
    """Check that optstatus has no UNCOVERGED values."""
    assert ccData.OPT_UNCONVERGED not in logfile.data.optstatus

    assert logfile.data.metadata["legacy_package_version"] == "09revisionD.01"
    assert logfile.data.metadata["package_version"] == "2009+D.01"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)

    # The energies printed in the scan summary are misformated.
    assert numpy.all(numpy.isnan(logfile.data.scanenergies))


def testGaussian_Gaussian09_2D_PES_one_unconverged_log(logfile):
    """Check that optstatus contains UNCOVERGED values."""
    assert ccData.OPT_UNCONVERGED in logfile.data.optstatus

    assert logfile.data.metadata["package_version"] == "2009+D.01"


def testGaussian_Gaussian09_534_out(logfile):
    """Previously, caused etenergies parsing to fail."""
    assert logfile.data.etsyms[0] == "Singlet-?Sym"
    assert abs(logfile.data.etenergies[0] - 0.09532039604871197) < 1.0e-5

    assert logfile.data.metadata["legacy_package_version"] == "09revisionA.02"
    assert logfile.data.metadata["package_version"] == "2009+A.02"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGaussian_Gaussian09_BSL_opt_freq_DFT_out(logfile):
    """Failed for converting to CJSON when moments weren't parsed for
    Gaussian.
    """
    assert hasattr(logfile.data, "moments")
    # dipole Y
    assert logfile.data.moments[1][1] == 0.5009
    # hexadecapole ZZZZ
    assert logfile.data.moments[4][-1] == -77.9600

    assert logfile.data.metadata["package_version"] == "2009+D.01"


def testGaussian_Gaussian09_dvb_gopt_unconverged_log(logfile):
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert hasattr(logfile.data, "optdone") and not logfile.data.optdone
    assert logfile.data.optstatus[-1] == logfile.data.OPT_UNCONVERGED

    assert (
        logfile.data.metadata["cpu_time"]
        == logfile.data.metadata["wall_time"]
        == [datetime.timedelta(seconds=27, microseconds=700000)]
    )
    assert logfile.data.metadata["package_version"] == "2009+D.01"


def testGaussian_Gaussian09_dvb_lowdin_log(logfile):
    """Check if both Mulliken and Lowdin charges are parsed."""
    assert "mulliken" in logfile.data.atomcharges
    assert "lowdin" in logfile.data.atomcharges

    assert logfile.data.metadata["package_version"] == "2009+A.02"


def testGaussian_Gaussian09_Dahlgren_TS_log(logfile):
    """Failed to parse ccenergies for a variety of reasons"""
    assert hasattr(logfile.data, "ccenergies")
    assert abs(logfile.data.ccenergies[0] - (-434.37573219)) < 1.0e-6

    assert logfile.data.metadata["package_version"] == "2009+A.02"


def testGaussian_Gaussian09_irc_point_log(logfile):
    """Failed to parse vibfreqs except for 10, 11"""
    assert hasattr(logfile.data, "vibfreqs")
    assert len(logfile.data.vibfreqs) == 11

    assert logfile.data.metadata["package_version"] == "2009+D.01"


def testGaussian_Gaussian09_issue_460_log(logfile):
    """Lots of malformed lines when parsing for scfvalues:

    RMSDP=3.79D-04 MaxDP=4.02D-02              OVMax= 4.31D-02
    RMSDP=1.43D-06 MaxDP=5.44D-04 DE=-6.21D-07 OVMax= 5.76D-04
    RMSDP=2.06D-05 MaxDP=3.84D-03 DE= 4.82D-04 O E= -2574.14897924075     Delta-E=        0.000439804468 Rises=F Damp=F
    RMSDP=8.64D-09 MaxDP=2.65D-06 DE=-1.67D-10 OVMax= 3. E= -2574.14837678675     Delta-E=       -0.000000179038 Rises=F Damp=F
    RMSDP= E= -2574.14931865182     Delta-E=       -0.000000019540 Rises=F Damp=F
    RMSDP=9.34D- E= -2574.14837612206     Delta-E=       -0.000000620705 Rises=F Damp=F
    RMSDP=7.18D-05 Max E= -2574.14797761904     Delta-E=       -0.000000000397 Rises=F Damp=F
    RMSDP=1.85D-06 MaxD E= -2574.14770506975     Delta-E=       -0.042173156160 Rises=F Damp=F
    RMSDP=1.69D-06 MaxDP= E= -2574.14801776548     Delta-E=        0.000023521317 Rises=F Damp=F
    RMSDP=3.80D-08 MaxDP=1 E= -2574.14856570920     Delta-E=       -0.000002960194 Rises=F Damp=F
    RMSDP=4.47D-09 MaxDP=1.40 E= -2574.14915435699     Delta-E=       -0.000255709558 Rises=F Damp=F
    RMSDP=5.54D-08 MaxDP=1.55D-05 DE=-2.55D-0 E= -2574.14854319757     Delta-E=       -0.000929740010 Rises=F Damp=F
    RMSDP=7.20D-09 MaxDP=1.75D-06 DE=- (Enter /QFsoft/applic/GAUSSIAN/g09d.01_pgi11.9-ISTANBUL/g09/l703.exe)
    RMSDP=5.24D-09 MaxDP=1.47D-06 DE=-1.82D-11 OVMax= 2.15 (Enter /QFsoft/applic/GAUSSIAN/g09d.01_pgi11.9-ISTANBUL/g09/l703.exe)
    RMSDP=1.71D-04 MaxDP=1.54D-02    Iteration    2 A^-1*A deviation from unit magnitude is 1.11D-15 for    266.
    """
    assert hasattr(logfile.data, "scfvalues")
    assert logfile.data.scfvalues[0][0, 0] == 3.37e-03
    assert numpy.isnan(logfile.data.scfvalues[0][0, 2])

    assert logfile.data.metadata["package_version"] == "2009+D.01"


def testGaussian_Gaussian09_OPT_td_g09_out(logfile):
    """Couldn't find etrotats as G09 has different output than G03."""
    assert len(logfile.data.etrotats) == 10
    assert logfile.data.etrotats[0] == -0.4568

    assert logfile.data.metadata["package_version"] == "2009+A.02"


def testGaussian_Gaussian09_OPT_oniom_log(logfile):
    """AO basis extraction broke with ONIOM"""

    assert logfile.data.metadata["package_version"] == "2009+D.01"


def testGaussian_Gaussian09_oniom_IR_intensity_log(logfile):
    """Problem parsing IR intensity from mode 192"""
    assert hasattr(logfile.data, "vibirs")
    assert len(logfile.data.vibirs) == 216

    assert logfile.data.metadata["package_version"] == "2009+C.01"


def testGaussian_Gaussian09_Ru2bpyen2_H2_freq3_log(logfile):
    """Here atomnos wans't added to the gaussian parser before."""
    assert len(logfile.data.atomnos) == 69

    assert logfile.data.metadata["package_version"] == "2009+A.02"


def testGaussian_Gaussian09_benzene_HPfreq_log(logfile):
    """Check that higher precision vib displacements obtained with freq=hpmodes) are parsed correctly."""
    assert abs(logfile.data.vibdisps[0, 0, 2] - (-0.04497)) < 0.00001

    assert logfile.data.metadata["package_version"] == "2009+C.01"


def testGaussian_Gaussian09_benzene_freq_log(logfile):
    """Check that default precision vib displacements are parsed correctly."""
    assert abs(logfile.data.vibdisps[0, 0, 2] - (-0.04)) < 0.00001

    assert logfile.data.metadata["package_version"] == "2009+C.01"


def testGaussian_Gaussian09_relaxed_PES_testH2_log(logfile):
    """Check that all optimizations converge in a single step."""
    atomcoords = logfile.data.atomcoords
    optstatus = logfile.data.optstatus
    assert len(optstatus) == len(atomcoords)

    assert all(s == ccData.OPT_DONE + ccData.OPT_NEW for s in optstatus)

    assert logfile.data.metadata["package_version"] == "2009+D.01"


def testGaussian_Gaussian09_relaxed_PES_testCO2_log(logfile):
    """A relaxed PES scan with some uncoverged and some converged runs."""
    atomcoords = logfile.data.atomcoords
    optstatus = logfile.data.optstatus
    assert len(optstatus) == len(atomcoords)

    new_points = numpy.where(optstatus & ccData.OPT_NEW)[0]

    # The first new point is just the beginning of the scan.
    assert new_points[0] == 0

    # The next two new points are at the end of unconverged runs.
    assert optstatus[new_points[1] - 1] == ccData.OPT_UNCONVERGED
    assert all(
        optstatus[i] == ccData.OPT_UNKNOWN for i in range(new_points[0] + 1, new_points[1] - 1)
    )
    assert optstatus[new_points[2] - 1] == ccData.OPT_UNCONVERGED
    assert all(
        optstatus[i] == ccData.OPT_UNKNOWN for i in range(new_points[1] + 1, new_points[2] - 1)
    )

    # The next new point is after a convergence.
    assert optstatus[new_points[3] - 1] == ccData.OPT_DONE
    assert all(
        optstatus[i] == ccData.OPT_UNKNOWN for i in range(new_points[2] + 1, new_points[3] - 1)
    )

    # All subsequent point are both new and converged, since they seem
    # to have converged in a single step.
    assert all(s == ccData.OPT_DONE + ccData.OPT_NEW for s in optstatus[new_points[3] :])

    assert logfile.data.metadata["package_version"] == "2009+D.01"


def testGaussian_Gaussian09_stopiter_gaussian_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert len(logfile.data.scfvalues[0]) == 4

    assert logfile.data.metadata["package_version"] == "2009+D.01"


def testGaussian_Gaussian09_benzene_excited_states_optimization_issue889_log(logfile):
    """Check that only converged geometry excited states properties are reported."""
    assert logfile.data.etdips.shape == (20, 3)
    assert len(logfile.data.etenergies) == 20
    assert logfile.data.etmagdips.shape == (20, 3)
    assert len(logfile.data.etoscs) == 20
    assert len(logfile.data.etrotats) == 20
    assert len(logfile.data.etsecs) == 20
    assert logfile.data.etveldips.shape == (20, 3)


def testGaussian_Gaussian09_issue1150_log(logfile):
    """Symmetry parsing for Gaussian09 was broken"""
    assert logfile.metadata["symmetry_detected"] == "c1"


def testGaussian_Gaussian16_co_pbe1pbe_631ppGss_log(logfile):
    """A calculation where the platform couldn't be parsed."""
    # *********************************************
    # Gaussian 16:  Apple M1-G16RevC.02  7-Dec-2021
    #                   6-Mar-2024
    # *********************************************
    assert logfile.metadata["package_version"] == "2016+C.02"
    assert logfile.metadata["platform"] == "Apple M1"


def testGaussian_Gaussian16_H3_natcharge_log(logfile):
    """A calculation with natural charges calculated. Test issue 1055 where
    only the beta set of charges was parsed rather than the spin independent"""

    assert isinstance(logfile.data.atomcharges, dict)
    assert "mulliken" in logfile.data.atomcharges
    assert "natural" in logfile.data.atomcharges
    assert numpy.all(logfile.data.atomcharges["natural"] == [0.0721, -0.1442, 0.0721])


def testGaussian_Gaussian16_naturalspinorbitals_parsing_log(logfile):
    """A UHF calculation with natural spin orbitals."""

    assert isinstance(logfile.data.nocoeffs, numpy.ndarray)
    assert isinstance(logfile.data.nooccnos, numpy.ndarray)
    assert isinstance(logfile.data.aonames, list)
    assert isinstance(logfile.data.atombasis, list)

    assert numpy.shape(logfile.data.nocoeffs) == (2, logfile.data.nmo, logfile.data.nmo)
    assert len(logfile.data.nooccnos[0]) == logfile.data.nmo
    assert len(logfile.data.nooccnos[1]) == logfile.data.nmo
    assert len(logfile.data.aonames) == logfile.data.nbasis
    assert len(numpy.ravel(logfile.data.atombasis)) == logfile.data.nbasis

    assert logfile.data.nooccnos[0][14] == 0.00506
    assert logfile.data.nooccnos[1][14] == 0.00318
    assert logfile.data.nocoeffs[0][14, 12] == 0.00618
    assert logfile.data.nocoeffs[1][14, 9] == 0.79289
    assert logfile.data.aonames[41] == "O2_9D 0"
    assert logfile.data.atombasis[1][0] == 23

    assert logfile.data.metadata["cpu_time"] == [
        datetime.timedelta(seconds=74, microseconds=400000)
    ]
    assert logfile.data.metadata["wall_time"] == [
        datetime.timedelta(seconds=3, microseconds=500000)
    ]
    assert logfile.data.metadata["legacy_package_version"] == "16revisionA.03"
    assert logfile.data.metadata["package_version"] == "2016+A.03"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testGaussian_Gaussian16_issue851_log(logfile):
    """Surface scan from cclib/cclib#851 where attributes were not lists."""

    assert isinstance(logfile.data.scannames, list)
    assert isinstance(logfile.data.scanparm, list)
    assert isinstance(logfile.data.scanenergies, list)


def testGaussian_Gaussian16_issue962_log(logfile):
    """For issue 962, this shouldn't have scftargets but should parse fully"""

    assert not hasattr(logfile.data, "scftargets")


def testGaussian_Gaussian16_C01_CC_log(logfile):
    """For issue 1110, check parsing of ccenergies in newer Gaussian version"""

    assert hasattr(logfile.data, "ccenergies")


def testGaussian_Gaussian16_Ethane_mp5_log(logfile):
    """For issue 1163, check we can parse a log file that has MPn in its description."""

    # This issue is about failing to parse if certain strings are present in the Gaussian log file description section.
    # Check we can still parse MP energies up to MP5
    assert hasattr(logfile.data, "mpenergies")
    assert len(logfile.data.mpenergies) == 1
    assert len(logfile.data.mpenergies[0]) == 4


# Jaguar #

# It would be good to have an unconverged geometry optimization so that
# we can test that optdone is set properly.
# def testJaguarX.X_dvb_gopt_unconverged:
#    assert hasattr(logfile.data, 'optdone') and not logfile.data.optdone


def testJaguar_Jaguar8_3_stopiter_jaguar_dft_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert len(logfile.data.scfvalues[0]) == 4

    assert logfile.data.metadata["legacy_package_version"] == "8.3"
    assert logfile.data.metadata["package_version"] == "8.3+13"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testJaguar_Jaguar8_3_stopiter_jaguar_hf_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert len(logfile.data.scfvalues[0]) == 3

    assert logfile.data.metadata["package_version"] == "8.3+13"


# Molcas #


def testMolcas_Molcas18_test_standard_000_out(logfile):
    """Don't support parsing MOs for multiple symmetry species."""
    assert not hasattr(logfile.data, "moenergies")
    assert not hasattr(logfile.data, "mocoeffs")

    assert logfile.data.metadata["legacy_package_version"] == "18.09"
    assert (
        logfile.data.metadata["package_version"]
        == "18.09+52-ge15dc38.81d3fb3dc6a5c5df6b3791ef1ef3790f"
    )
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testMolcas_Molcas18_test_standard_001_out(logfile):
    """This logfile has two calculations, and we currently only want to parse the first."""
    assert logfile.data.natom == 8

    # There are also four symmetry species, and orbital count should cover all of them.
    assert logfile.data.nbasis == 30
    assert logfile.data.nmo == 30

    assert (
        logfile.data.metadata["package_version"]
        == "18.09+52-ge15dc38.81d3fb3dc6a5c5df6b3791ef1ef3790f"
    )


def testMolcas_Molcas18_test_standard_003_out(logfile):
    """This logfile has extra charged monopoles (not part of the molecule)."""
    assert logfile.data.charge == 0

    assert (
        logfile.data.metadata["package_version"]
        == "18.09+52-ge15dc38.81d3fb3dc6a5c5df6b3791ef1ef3790f"
    )


def testMolcas_Molcas18_test_standard_005_out(logfile):
    """Final geometry in optimization has fewer atoms due to symmetry, and so is ignored."""
    assert len(logfile.data.atomcoords) == 2

    assert (
        logfile.data.metadata["package_version"]
        == "18.09+52-ge15dc38.81d3fb3dc6a5c5df6b3791ef1ef3790f"
    )


def testMolcas_Molcas18_test_stevenv_001_out(logfile):
    """Don't support parsing MOs for RAS (active space)."""
    assert not hasattr(logfile.data, "moenergies")
    assert not hasattr(logfile.data, "mocoeffs")

    assert (
        logfile.data.metadata["package_version"]
        == "18.09+52-ge15dc38.81d3fb3dc6a5c5df6b3791ef1ef3790f"
    )


def testMolcas_Molcas18_test_stevenv_desym_out(logfile):
    """This logfile has iterations interrupted by a Fermi aufbau procedure."""
    assert len(logfile.data.scfvalues) == 1
    assert len(logfile.data.scfvalues[0]) == 26

    assert (
        logfile.data.metadata["package_version"]
        == "18.09+52-ge15dc38.81d3fb3dc6a5c5df6b3791ef1ef3790f"
    )


# Molpro #


def testMolpro_Molpro2008_ch2o_molpro_casscf_out(logfile):
    """A CASSCF job with symmetry and natural orbitals."""

    # The last two atoms are equivalent, so the last ends up having no
    # functions asigned. This is not obvious, because the functions are
    # distributed between the last two atoms in the block where gbasis
    # is parsed, but it seems all are assigned to the penultimate atom later.
    assert logfile.data.atombasis[-1] == []
    assert len(logfile.data.aonames) == logfile.data.nbasis

    # The MO coefficients are printed in several block, each corresponding
    # to one irrep, so make sure we have reconstructed the coefficients correctly.
    assert len(logfile.data.moenergies) == 1
    assert logfile.data.moenergies[0].shape == (logfile.data.nmo,)
    assert len(logfile.data.mocoeffs) == 1
    assert logfile.data.mocoeffs[0].shape == (logfile.data.nmo, logfile.data.nbasis)

    # These coefficients should be zero due to symmetry.
    assert logfile.data.mocoeffs[0][-2][0] == 0.0
    assert logfile.data.mocoeffs[0][0][-2] == 0.0

    assert isinstance(logfile.data.nocoeffs, numpy.ndarray)
    assert isinstance(logfile.data.nooccnos, numpy.ndarray)
    assert logfile.data.nocoeffs.shape == logfile.data.mocoeffs[0].shape
    assert len(logfile.data.nooccnos) == logfile.data.nmo
    assert logfile.data.nooccnos[27] == 1.95640

    assert logfile.data.metadata["legacy_package_version"] == "2012.1"
    assert logfile.data.metadata["package_version"] == "2012.1"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testMolpro_Molpro2012_CHONHSH_HF_STO_3G_out(logfile):
    """Formatting of the basis function is slightly different than expected."""
    assert len(logfile.data.gbasis) == 7
    assert len(logfile.data.gbasis[0]) == 3  # C
    assert len(logfile.data.gbasis[1]) == 3  # N
    assert len(logfile.data.gbasis[2]) == 3  # O
    assert len(logfile.data.gbasis[3]) == 5  # S
    assert len(logfile.data.gbasis[4]) == 1  # H
    assert len(logfile.data.gbasis[5]) == 1  # H
    assert len(logfile.data.gbasis[6]) == 1  # H

    assert logfile.data.metadata["legacy_package_version"] == "2012.1"
    assert (
        logfile.data.metadata["package_version"]
        == "2012.1.23+f8cfea266908527a8826bdcd5983aaf62e47d3bf"
    )
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testMolpro_Molpro2012_dvb_gopt_unconverged_out(logfile):
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert hasattr(logfile.data, "optdone") and not logfile.data.optdone

    assert logfile.data.metadata["legacy_package_version"] == "2012.1"
    assert (
        logfile.data.metadata["package_version"]
        == "2012.1.12+e112a8ab93d81616c1987a1f1ef3707d874b6803"
    )
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testMolpro_Molpro2012_stopiter_molpro_dft_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert len(logfile.data.scfvalues[0]) == 6

    assert logfile.data.metadata["legacy_package_version"] == "2012.1"
    assert (
        logfile.data.metadata["package_version"]
        == "2012.1+c18f7d37f9f045f75d4f3096db241dde02ddca0a"
    )
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testMolpro_Molpro2012_stopiter_molpro_hf_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert len(logfile.data.scfvalues[0]) == 6

    assert (
        logfile.data.metadata["package_version"]
        == "2012.1+c18f7d37f9f045f75d4f3096db241dde02ddca0a"
    )


# MOPAC #


def testMOPAC_MOPAC2016_9S3_uuu_Cs_cation_freq_PM7_out(logfile):
    """There was a syntax error in the frequency parsing."""
    assert hasattr(logfile.data, "vibfreqs")

    assert logfile.data.metadata["legacy_package_version"] == "2016"
    assert logfile.data.metadata["package_version"] == "16.175"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


# NWChem #


def testNWChem_NWChem6_0_dvb_gopt_hf_unconverged_out(logfile):
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert hasattr(logfile.data, "optdone") and not logfile.data.optdone

    assert logfile.data.metadata["legacy_package_version"] == "6.0"
    assert logfile.data.metadata["package_version"] == "6.0"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testNWChem_NWChem6_0_dvb_sp_hf_moments_only_quadrupole_out(logfile):
    """Quadrupole moments are printed/parsed, but not lower moments (no shape)."""
    assert hasattr(logfile.data, "moments") and len(logfile.data.moments) == 3
    assert len(logfile.data.moments[0]) == 3
    assert not logfile.data.moments[1].shape
    assert len(logfile.data.moments[2]) == 6

    assert logfile.data.metadata["package_version"] == "6.0"


def testNWChem_NWChem6_0_dvb_sp_hf_moments_only_octupole_out(logfile):
    """Quadrupole moments are printed/parsed, but not lower moments (no shape)."""
    assert hasattr(logfile.data, "moments") and len(logfile.data.moments) == 4
    assert len(logfile.data.moments[0]) == 3
    assert not logfile.data.moments[1].shape
    assert not logfile.data.moments[2].shape
    assert len(logfile.data.moments[3]) == 10

    assert logfile.data.metadata["package_version"] == "6.0"


def testNWChem_NWChem6_0_hydrogen_atom_ROHF_cc_pVDZ_out(logfile):
    """A lone hydrogen atom is a common edge case; it has no beta
    electrons.
    """
    assert logfile.data.charge == 0
    assert logfile.data.natom == 1
    assert logfile.data.nbasis == 5
    assert logfile.data.nmo == 5
    assert len(logfile.data.moenergies) == 1
    assert logfile.data.moenergies[0].shape == (5,)
    assert logfile.data.homos.shape == (2,)
    assert logfile.data.homos[0] == 0
    assert logfile.data.homos[1] == -1

    assert logfile.data.metadata["package_version"] == "6.0"


def testNWChem_NWChem6_0_hydrogen_atom_UHF_cc_pVDZ_out(logfile):
    """A lone hydrogen atom is a common edge case; it has no beta
    electrons.

    Additionally, this calculations has no title, which caused some
    issues with skip_lines().
    """
    assert logfile.data.charge == 0
    assert logfile.data.natom == 1
    assert logfile.data.nbasis == 5
    assert logfile.data.nmo == 5
    assert len(logfile.data.moenergies) == 2
    assert logfile.data.moenergies[0].shape == (5,)
    assert logfile.data.moenergies[1].shape == (5,)
    assert logfile.data.homos.shape == (2,)
    assert logfile.data.homos[0] == 0
    assert logfile.data.homos[1] == -1

    assert logfile.data.metadata["package_version"] == "6.0"


def testNWChem_NWChem6_5_stopiter_nwchem_dft_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert len(logfile.data.scfvalues[0]) == 3

    assert logfile.data.metadata["legacy_package_version"] == "6.5"
    assert logfile.data.metadata["package_version"] == "6.5+26243"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testNWChem_NWChem6_5_stopiter_nwchem_hf_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert len(logfile.data.scfvalues[0]) == 2

    assert logfile.data.metadata["package_version"] == "6.5+26243"


# def testNWChem_NWChem6_8_1057_out(logfile):
#    """Multistep job caused a premature end of parsing."""
#    assert not hasattr(logfile.data, "atomnos")
#    assert not hasattr(logfile.data, "gbasis")
#
#    assert logfile.data.metadata["legacy_package_version"] == "6.8.1"
#    assert logfile.data.metadata["package_version"] == "6.8.1+g2272a644e"
#    assert isinstance(
#        parse_version(logfile.data.metadata["package_version"]), Version
#    )


def testNWChem_NWChem6_8_526_out(logfile):
    """If `print low` is present in the input, SCF iterations are not
    printed.
    """
    assert not hasattr(logfile.data, "scftargets")
    assert not hasattr(logfile.data, "scfvalues")

    assert logfile.data.metadata["legacy_package_version"] == "6.8.1"
    assert logfile.data.metadata["package_version"] == "6.8.1+g08bf49b"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


# ORCA #


def testORCA_ORCA2_8_co_cosmo_out(logfile):
    """This is related to bug 3184890.

    The scfenergies were not being parsed correctly for this geometry
    optimization run, for two reasons.
    First, the printing of SCF total energies is different inside
    geometry optimization steps than for single point calculations,
    which also affects unit tests.
    However, this logfile uses a setting that causes an SCF run to
    terminate prematurely when a set maximum number of cycles is reached.
    In this case, the last energy reported should probably be used,
    and the number of values in scfenergies preserved.
    """
    assert hasattr(logfile.data, "scfenergies") and len(logfile.data.scfenergies) == 4

    assert logfile.data.metadata["legacy_package_version"] == "2.8"
    assert logfile.data.metadata["package_version"] == "2.8+2287"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testORCA_ORCA2_9_job_out(logfile):
    """First output file and request to parse atomic spin densities.

    Make sure that the sum of such densities is one in this case (or reasonaby close),
    but remember that this attribute is a dictionary, so we must iterate.
    """
    assert all([abs(sum(v) - 1.0) < 0.0001 for k, v in logfile.data.atomspins.items()])

    assert logfile.data.metadata["legacy_package_version"] == "2.9.0"
    assert logfile.data.metadata["package_version"] == "2.9.0"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testORCA_ORCA2_9_qmspeedtest_hf_out(logfile):
    """Check precision of SCF energies (cclib/cclib#210)."""
    energy = logfile.data.scfenergies[-1]
    expected = -644.675706036271
    assert abs(energy - expected) < 1.0e-8

    assert logfile.data.metadata["legacy_package_version"] == "2.9.1"
    assert logfile.data.metadata["package_version"] == "2.9.1"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testORCA_ORCA3_0_chelpg_out(logfile):
    """ORCA file with chelpg charges"""
    assert "chelpg" in logfile.data.atomcharges
    charges = logfile.data.atomcharges["chelpg"]
    assert len(charges) == 9
    assert charges[0] == 0.363939
    assert charges[1] == 0.025695


def testORCA_ORCA3_0_dvb_gopt_unconverged_out(logfile):
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert hasattr(logfile.data, "optdone") and not logfile.data.optdone

    assert logfile.data.metadata["legacy_package_version"] == "3.0.1"
    assert logfile.data.metadata["package_version"] == "3.0.1"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testORCA_ORCA3_0_polar_rhf_cg_out(logfile):
    """Alternative CP-SCF solver for the polarizability wasn't being detected."""
    assert hasattr(logfile.data, "polarizabilities")

    assert logfile.data.metadata["legacy_package_version"] == "3.0.3"
    assert logfile.data.metadata["package_version"] == "3.0.3"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testORCA_ORCA3_0_polar_rhf_diis_out(logfile):
    """Alternative CP-SCF solver for the polarizability wasn't being detected."""
    assert hasattr(logfile.data, "polarizabilities")

    assert logfile.data.metadata["package_version"] == "3.0.3"


def testORCA_ORCA3_0_stopiter_orca_scf_compact_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert len(logfile.data.scfvalues[0]) == 1

    assert logfile.data.metadata["package_version"] == "3.0.1"


def testORCA_ORCA3_0_stopiter_orca_scf_large_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert len(logfile.data.scfvalues[0]) == 9

    assert logfile.data.metadata["package_version"] == "2.9.1"


def testORCA_ORCA4_0_1_ttt_td_out(logfile):
    """RPA is slightly different from TDA, see #373."""
    assert hasattr(logfile.data, "etsyms")
    assert len(logfile.data.etsecs) == 24
    assert len(logfile.data.etsecs[0]) == 1
    assert numpy.isnan(logfile.data.etsecs[0][0][2])
    assert len(logfile.data.etrotats) == 24
    assert logfile.data.etrotats[13] == -0.03974

    assert logfile.data.metadata["legacy_package_version"] == "4.0.0"
    assert logfile.data.metadata["package_version"] == "4.0.0"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testORCA_ORCA4_0_hydrogen_fluoride_numfreq_out(logfile):
    """Frequencies from linear molecules weren't parsed correctly (#426)."""
    numpy.testing.assert_equal(logfile.data.vibfreqs, [4473.96])


def testORCA_ORCA4_0_hydrogen_fluoride_usesym_anfreq_out(logfile):
    """Frequencies from linear molecules weren't parsed correctly (#426)."""
    numpy.testing.assert_equal(logfile.data.vibfreqs, [4473.89])


def testORCA_ORCA4_0_invalid_literal_for_float_out(logfile):
    """MO coefficients are glued together, see #629."""
    assert hasattr(logfile.data, "mocoeffs")
    assert logfile.data.mocoeffs[0].shape == (logfile.data.nmo, logfile.data.nbasis)

    # Test the coefficients from this line where things are glued together:
    # 15C   6s       -154.480939-111.069870-171.460819-79.052025241.536860-92.159399
    assert logfile.data.mocoeffs[0][102][378] == -154.480939
    assert logfile.data.mocoeffs[0][103][378] == -111.069870
    assert logfile.data.mocoeffs[0][104][378] == -171.460819
    assert logfile.data.mocoeffs[0][105][378] == -79.052025
    assert logfile.data.mocoeffs[0][106][378] == 241.536860
    assert logfile.data.mocoeffs[0][107][378] == -92.159399

    assert logfile.data.metadata["legacy_package_version"] == "4.0.1.2"
    assert logfile.data.metadata["package_version"] == "4.0.1.2"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testORCA_ORCA4_0_IrCl6_sp_out(logfile):
    """Tests ECP and weird SCF printing."""
    assert hasattr(logfile.data, "scfvalues")
    assert len(logfile.data.scfvalues) == 1
    vals_first = [0.000000000000, 28.31276975, 0.71923638]
    vals_last = [0.000037800796, 0.00412549, 0.00014041]
    numpy.testing.assert_almost_equal(logfile.data.scfvalues[0][0], vals_first)
    numpy.testing.assert_almost_equal(logfile.data.scfvalues[0][-1], vals_last)


def testORCA_ORCA4_0_comment_or_blank_line_out(logfile):
    """Coordinates with blank lines or comments weren't parsed correctly (#747)."""
    assert hasattr(logfile.data, "atomcoords")
    assert logfile.data.atomcoords.shape == (1, 8, 3)

    assert logfile.data.metadata["legacy_package_version"] == "4.0.1.2"
    assert logfile.data.metadata["package_version"] == "4.0.1.2"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testORCA_ORCA4_1_725_out(logfile):
    """This file uses embedding potentials, which requires `>` after atom names in
    the input file and that confuses different parts of the parser.

    In #725 we decided to not include these potentials in the parsed results.
    """
    assert logfile.data.natom == 7
    numpy.testing.assert_equal(
        logfile.data.atomnos, numpy.array([20, 17, 17, 17, 17, 17, 17], dtype=int)
    )
    assert len(logfile.data.atomcharges["mulliken"]) == 7
    assert len(logfile.data.atomcharges["lowdin"]) == 7

    assert logfile.data.metadata["legacy_package_version"] == "4.1.x"
    assert logfile.data.metadata["package_version"] == "4.1dev+13440"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testORCA_ORCA4_1_orca_from_issue_736_out(logfile):
    """ORCA file with no whitespace between SCF iteration columns."""
    assert len(logfile.data.scfvalues) == 23
    # The first iteration in the problematic block:
    # ITER       Energy         Delta-E        Max-DP      RMS-DP      [F,P]     Damp
    #           ***  Starting incremental Fock matrix formation  ***
    # 0   -257.0554667435   0.000000000000537.42184135  4.76025534  0.4401076 0.8500
    assert abs(logfile.data.scfvalues[14][0][1] - 537) < 1.0, logfile.data.scfvalues[14][0]


def testORCA_ORCA4_1_porphine_out(logfile):
    """ORCA optimization with multiple TD-DFT gradients and absorption spectra."""
    assert len(logfile.data.etenergies) == 1


def testORCA_ORCA4_1_single_atom_freq_out(logfile):
    """ORCA frequency with single atom."""
    assert len(logfile.data.vibdisps) == 0
    assert len(logfile.data.vibfreqs) == 0
    assert len(logfile.data.vibirs) == 0
    # These values are different from what ORCA prints as the total enthalpy,
    # because for single atoms that includes a spurious correction. We build the
    # enthalpy ourselves from electronic and translational energies (see #817 for details).
    numpy.testing.assert_almost_equal(logfile.data.enthalpy, -460.14376, 5)
    numpy.testing.assert_almost_equal(logfile.data.entropy, 6.056e-5, 8)
    numpy.testing.assert_almost_equal(logfile.data.freeenergy, -460.16182, 6)


def testORCA_ORCA4_2_947_out(logfile):
    """A constrained geometry optimization which prints the extra line

    WARNING: THERE ARE 5 CONSTRAINED CARTESIAN COORDINATES

    just before the gradient.
    """
    assert len(logfile.data.atomcoords) == 7
    assert len(logfile.data.grads) == 6


def testORCA_ORCA4_2_MP2_gradient_out(logfile):
    """ORCA numerical frequency calculation with gradients."""
    assert logfile.data.metadata["package_version"] == "4.2.0"
    assert hasattr(logfile.data, "grads")
    assert logfile.data.grads.shape == (1, 3, 3)
    # atom 2, y-coordinate.
    idx = (0, 1, 1)
    assert logfile.data.grads[idx] == -0.00040549


def testORCA_ORCA4_2_ligando_30_SRM1_S_ZINDO_out(logfile):
    """ORCA says that ZINDO uses the def2-SVP basis set before echoing the
    input file despite actually using STO-3G fit to Slater functions (#1187).
    """
    assert logfile.data.metadata["basis_set"] == "STO-3G"
    assert logfile.data.metadata["methods"] == ["ZINDO/S"]
    assert hasattr(logfile.data, "etsyms")


def testORCA_ORCA4_2_long_input_out(logfile):
    """Long ORCA input file (#804)."""
    assert logfile.data.metadata["package_version"] == "4.2.0"
    assert hasattr(logfile.data, "atomcoords")
    assert logfile.data.atomcoords.shape == (100, 12, 3)


def testORCA_ORCA4_2_water_dlpno_ccsd_out(logfile):
    """DLPNO-CCSD files have extra lines between E(0) and E(TOT) than normal CCSD
    outputs:

        ----------------------
        COUPLED CLUSTER ENERGY
        ----------------------

        E(0)                                       ...    -74.963574242
        E(CORR)(strong-pairs)                      ...     -0.049905771
        E(CORR)(weak-pairs)                        ...      0.000000000
        E(CORR)(corrected)                         ...     -0.049905771
        E(TOT)                                     ...    -75.013480013
        Singles Norm <S|S>**1/2                    ...      0.013957180
        T1 diagnostic                              ...      0.004934608
    """
    assert hasattr(logfile.data, "ccenergies")


def testORCA_ORCA4_2_longer_input_out(logfile):
    """Longer ORCA input file (#1034)."""
    assert (
        logfile.data.metadata["input_file_contents"][-47:-4]
        == "H   1.066878310   1.542378768  -0.602599044"
    )


def testORCA_ORCA4_2_casscf_out(logfile):
    """ORCA casscf input file (#1044)."""
    assert numpy.isclose(logfile.data.etenergies[0], 0.128812)


def testORCA_ORCA5_0_ADBNA_Me_Mes_MesCz_log(logfile):
    """Check we can parse etsyms in difficult cases."""
    assert hasattr(logfile.data, "etsyms")


def testORCA_ORCA5_0_Benzene_opt_etsyms_log(logfile):
    """Check we can parse etsyms in opt + excited states calc."""
    assert hasattr(logfile.data, "etsyms")


# PSI 3 #


def testPsi3_Psi3_4_water_psi3_log(logfile):
    """An RHF for water with D orbitals and C2v symmetry.

    Here we can check that the D orbitals are considered by checking atombasis and nbasis.
    """
    assert logfile.data.nbasis == 25
    assert [len(ab) for ab in logfile.data.atombasis] == [15, 5, 5]

    # FIXME not present? wasn't failing earlier?
    # assert logfile.data.metadata["legacy_package_version"] == "3.4"
    assert logfile.data.metadata["package_version"] == "3.4alpha"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


# PSI 4 #


def testPsi4_Psi4_beta5_dvb_gopt_hf_unconverged_out(logfile):
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert logfile.data.metadata["legacy_package_version"] == "beta5"
    assert logfile.data.metadata["package_version"] == "0!0.beta5"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)
    assert hasattr(logfile.data, "optdone") and not logfile.data.optdone


def testPsi4_Psi4_beta5_sample_cc54_0_01_0_1_0_1_out(logfile):
    """TODO"""
    assert logfile.data.metadata["legacy_package_version"] == "beta2+"
    assert (
        logfile.data.metadata["package_version"]
        == "0!0.beta2.dev+fa5960b375b8ca2a5e4000a48cb95e7f218c579a"
    )
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testPsi4_Psi4_beta5_stopiter_psi_dft_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert logfile.data.metadata["package_version"] == "0!0.beta5"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)
    assert len(logfile.data.scfvalues[0]) == 7


def testPsi4_Psi4_beta5_stopiter_psi_hf_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert logfile.data.metadata["package_version"] == "0!0.beta5"
    assert len(logfile.data.scfvalues[0]) == 6


def testPsi4_Psi4_0_5_sample_scf5_out(logfile):
    assert logfile.data.metadata["legacy_package_version"] == "0.5"
    assert logfile.data.metadata["package_version"] == "1!0.5.dev+master-dbe9080"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testPsi4_Psi4_0_5_water_fdgrad_out(logfile):
    """Ensure that finite difference gradients are parsed."""
    assert logfile.data.metadata["legacy_package_version"] == "1.2a1.dev429"
    assert logfile.data.metadata["package_version"] == "1!1.2a1.dev429+fixsym-7838fc1-dirty"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)
    assert hasattr(logfile.data, "grads")
    assert logfile.data.grads.shape == (1, 3, 3)
    assert abs(logfile.data.grads[0, 0, 2] - 0.05498126903657) < 1.0e-12
    # In C2v symmetry, there are 5 unique displacements for the
    # nuclear gradient, and this is at the MP2 level.
    assert logfile.data.mpenergies.shape == (5, 1)


def testPsi4_Psi4_1_2_ch4_hf_opt_freq_out(logfile):
    """Ensure that molecular orbitals and normal modes are parsed in Psi4 1.2"""
    assert logfile.data.metadata["legacy_package_version"] == "1.2.1"
    assert logfile.data.metadata["package_version"] == "1!1.2.1.dev+HEAD-406f4de"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)
    assert hasattr(logfile.data, "mocoeffs")
    assert hasattr(logfile.data, "vibdisps")
    assert hasattr(logfile.data, "vibfreqs")


# Q-Chem #


def testQChem_QChem4_2_CH3___Na__RS_out(logfile):
    """An unrestricted fragment job with BSSE correction.

    Contains only the Roothaan step energies for the CP correction.

    The fragment SCF sections are printed.

    This is to ensure only the supersystem is parsed.
    """

    assert logfile.data.metadata["legacy_package_version"] == "4.2.2"
    assert logfile.data.metadata["package_version"] == "4.2.2"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)

    assert logfile.data.charge == 1
    assert logfile.data.mult == 2
    assert len(logfile.data.moenergies) == 2
    assert len(logfile.data.atomcoords[0]) == 5
    assert len(logfile.data.atomnos) == 5

    # Fragments: A, B, RS_CP(A), RS_CP(B), Full
    assert len(logfile.data.scfenergies) == 1
    scfenergy = -201.9388745658
    assert abs(logfile.data.scfenergies[0] - scfenergy) < 1.0e-10

    assert logfile.data.nbasis == logfile.data.nmo == 40
    assert len(logfile.data.moenergies[0]) == 40
    assert len(logfile.data.moenergies[1]) == 40
    assert isinstance(logfile.data.moenergies, list)
    assert isinstance(logfile.data.moenergies[0], numpy.ndarray)
    assert isinstance(logfile.data.moenergies[1], numpy.ndarray)


def testQChem_QChem4_2_CH3___Na__RS_SCF_out(logfile):
    """An unrestricted fragment job with BSSE correction.

    Contains both the Roothaan step and full SCF energies for the CP correction.

    The fragment SCF sections are printed.

    This is to ensure only the supersystem is printed.
    """

    assert logfile.data.metadata["legacy_package_version"] == "4.1.0.1"
    assert logfile.data.metadata["package_version"] == "4.1.0.1"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)

    assert logfile.data.charge == 1
    assert logfile.data.mult == 2
    assert len(logfile.data.moenergies) == 2
    assert len(logfile.data.atomcoords[0]) == 5
    assert len(logfile.data.atomnos) == 5

    # Fragments: A, B, RS_CP(A), RS_CP(B), SCF_CP(A), SCF_CP(B), Full
    assert len(logfile.data.scfenergies) == 1
    scfenergy = -201.9396979324
    assert abs(logfile.data.scfenergies[0] - scfenergy) < 1.0e-10

    assert logfile.data.nbasis == logfile.data.nmo == 40
    assert len(logfile.data.moenergies[0]) == 40
    assert len(logfile.data.moenergies[1]) == 40
    assert isinstance(logfile.data.moenergies, list)
    assert isinstance(logfile.data.moenergies[0], numpy.ndarray)
    assert isinstance(logfile.data.moenergies[1], numpy.ndarray)


def testQChem_QChem4_2_CH4___Na__out(logfile):
    """A restricted fragment job with no BSSE correction.

    The fragment SCF sections are printed.

    This is to ensure only the supersystem is parsed.
    """

    assert logfile.data.metadata["legacy_package_version"] == "4.2.0"
    assert logfile.data.metadata["package_version"] == "4.2.0"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)

    assert logfile.data.charge == 1
    assert logfile.data.mult == 1
    assert len(logfile.data.moenergies) == 1
    assert len(logfile.data.atomcoords[0]) == 6
    assert len(logfile.data.atomnos) == 6

    # Fragments: A, B, Full
    assert len(logfile.data.scfenergies) == 1
    scfenergy = -202.6119443654
    assert abs(logfile.data.scfenergies[0] - scfenergy) < 1.0e-10

    assert logfile.data.nbasis == logfile.data.nmo == 42
    assert len(logfile.data.moenergies[0]) == 42
    assert isinstance(logfile.data.moenergies, list)
    assert isinstance(logfile.data.moenergies[0], numpy.ndarray)


def testQChem_QChem4_2_CH3___Na__RS_SCF_noprint_out(logfile):
    """An unrestricted fragment job with BSSE correction.

    Contains both the Roothaan step and full SCF energies for the CP correction.

    The fragment SCF sections are not printed.

    This is to ensure only the supersystem is parsed.
    """

    assert logfile.data.metadata["legacy_package_version"] == "4.3.0"
    assert logfile.data.metadata["package_version"] == "4.3.0"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)

    assert logfile.data.charge == 1
    assert logfile.data.mult == 2
    assert len(logfile.data.moenergies) == 2
    assert len(logfile.data.atomcoords[0]) == 5
    assert len(logfile.data.atomnos) == 5

    assert len(logfile.data.scfenergies) == 1
    scfenergy = -201.9396979324
    assert abs(logfile.data.scfenergies[0] - scfenergy) < 1.0e-10

    assert logfile.data.nbasis == logfile.data.nmo == 40
    assert len(logfile.data.moenergies[0]) == 40
    assert len(logfile.data.moenergies[1]) == 40
    assert isinstance(logfile.data.moenergies, list)
    assert isinstance(logfile.data.moenergies[0], numpy.ndarray)
    assert isinstance(logfile.data.moenergies[1], numpy.ndarray)


def testQChem_QChem4_2_CH3___Na__RS_noprint_out(logfile):
    """An unrestricted fragment job with BSSE correction.

    Contains only the Roothaan step energies for the CP correction.

    The fragment SCF sections are not printed.

    This is to ensure only the supersystem is parsed.
    """

    assert logfile.data.metadata["package_version"] == "4.3.0"

    assert logfile.data.charge == 1
    assert logfile.data.mult == 2
    assert len(logfile.data.moenergies) == 2
    assert len(logfile.data.atomcoords[0]) == 5
    assert len(logfile.data.atomnos) == 5

    assert len(logfile.data.scfenergies) == 1
    scfenergy = -201.9388582085
    assert abs(logfile.data.scfenergies[0] - scfenergy) < 1.0e-10

    assert logfile.data.nbasis == logfile.data.nmo == 40
    assert len(logfile.data.moenergies[0]) == 40
    assert len(logfile.data.moenergies[1]) == 40
    assert isinstance(logfile.data.moenergies, list)
    assert isinstance(logfile.data.moenergies[0], numpy.ndarray)
    assert isinstance(logfile.data.moenergies[1], numpy.ndarray)


def testQChem_QChem4_2_CH4___Na__noprint_out(logfile):
    """A restricted fragment job with no BSSE correction.

    The fragment SCF sections are not printed.

    This is to ensure only the supersystem is parsed.
    """

    assert logfile.data.metadata["package_version"] == "4.3.0"

    assert logfile.data.charge == 1
    assert logfile.data.mult == 1
    assert len(logfile.data.moenergies) == 1
    assert len(logfile.data.atomcoords[0]) == 6
    assert len(logfile.data.atomnos) == 6

    assert len(logfile.data.scfenergies) == 1
    scfenergy = -202.6119443654
    assert abs(logfile.data.scfenergies[0] - scfenergy) < 1.0e-10

    assert logfile.data.nbasis == logfile.data.nmo == 42
    assert len(logfile.data.moenergies[0]) == 42
    assert isinstance(logfile.data.moenergies, list)
    assert isinstance(logfile.data.moenergies[0], numpy.ndarray)


def testQChem_QChem4_2_CO2_out(logfile):
    """A job containing a specific number of orbitals requested for
    printing.
    """

    assert logfile.data.metadata["package_version"] == "4.2.2"

    nbasis = 45
    nmo = 45
    nalpha = 11
    assert logfile.data.nbasis == nbasis
    assert logfile.data.nmo == nmo
    assert len(logfile.data.mocoeffs) == 1
    assert logfile.data.mocoeffs[0].shape == (nmo, nbasis)
    assert logfile.data.mocoeffs[0][0, 0] == -0.0001434
    assert logfile.data.mocoeffs[0][nalpha + 5 - 1, nbasis - 1] == -0.0000661
    assert len(logfile.data.moenergies) == 1
    assert len(logfile.data.moenergies[0]) == nmo


def testQChem_QChem4_2_CO2_cation_UHF_out(logfile):
    """A job containing a specific number of orbitals requested for
    printing."""

    assert logfile.data.metadata["package_version"] == "4.2.2"

    nbasis = 45
    nmo = 45
    nalpha = 11
    nbeta = 10
    assert logfile.data.nbasis == nbasis
    assert logfile.data.nmo == nmo
    assert len(logfile.data.mocoeffs) == 2
    assert logfile.data.mocoeffs[0].shape == (nmo, nbasis)
    assert logfile.data.mocoeffs[1].shape == (nmo, nbasis)
    assert logfile.data.mocoeffs[0][0, 0] == -0.0001549
    assert logfile.data.mocoeffs[0][nalpha + 5 - 1, nbasis - 1] == -0.0000985
    assert logfile.data.mocoeffs[1][0, 0] == -0.0001612
    assert logfile.data.mocoeffs[1][nbeta + 5 - 1, nbasis - 1] == -0.0027710
    assert len(logfile.data.moenergies) == 2
    assert len(logfile.data.moenergies[0]) == nmo
    assert len(logfile.data.moenergies[1]) == nmo


def testQChem_QChem4_2_CO2_cation_ROHF_bigprint_allvirt_out(logfile):
    """A job containing a specific number of orbitals requested for
    printing."""

    assert logfile.data.metadata["package_version"] == "4.2.2"

    nbasis = 45
    nmo = 45
    nalpha = 11
    nbeta = 10
    assert logfile.data.nbasis == nbasis
    assert logfile.data.nmo == nmo
    assert len(logfile.data.mocoeffs) == 2
    assert logfile.data.mocoeffs[0].shape == (nmo, nbasis)
    assert logfile.data.mocoeffs[1].shape == (nmo, nbasis)
    assert logfile.data.mocoeffs[0][0, 0] == -0.0001543
    assert logfile.data.mocoeffs[0][nalpha + 5 - 3, nbasis - 1] == -0.0132848
    assert logfile.data.mocoeffs[1][2, 0] == 0.9927881
    assert logfile.data.mocoeffs[1][nbeta + 5 - 1, nbasis - 1] == 0.0018019
    assert len(logfile.data.moenergies) == 2
    assert len(logfile.data.moenergies[0]) == nmo
    assert len(logfile.data.moenergies[1]) == nmo


def testQChem_QChem4_2_CO2_linear_dependence_printall_out(logfile):
    """A job with linear dependency and all MOs printed."""

    assert logfile.data.metadata["package_version"] == "4.2.2"

    nbasis = 138
    nmo = 106
    assert logfile.data.nbasis == nbasis
    assert logfile.data.nmo == nmo
    assert len(logfile.data.mocoeffs) == 1
    assert logfile.data.mocoeffs[0].shape == (nmo, nbasis)
    assert logfile.data.mocoeffs[0].T[59, 15] == -0.28758
    assert logfile.data.mocoeffs[0].T[59, 16] == -0.00000


def testQChem_QChem4_2_CO2_linear_dependence_printall_final_out(logfile):
    """A job with linear dependency and all MOs printed.

    The increased precision is due to the presence of `scf_final_print
    = 3` giving a separate block with more decimal places.
    """

    assert logfile.data.metadata["package_version"] == "4.2.2"

    nbasis = 138
    nmo = 106
    assert logfile.data.nbasis == nbasis
    assert logfile.data.nmo == nmo
    assert len(logfile.data.mocoeffs) == 1
    assert logfile.data.mocoeffs[0].shape == (nmo, nbasis)
    assert logfile.data.mocoeffs[0].T[59, 15] == -0.2875844
    # Even though all MO coefficients are printed in the less precise
    # block, they aren't parsed.
    # assert logfile.data.mocoeffs[0].T[59, 16] == -0.00000
    assert numpy.isnan(logfile.data.mocoeffs[0].T[59, 16])


def testQChem_QChem4_2_CO2_linear_dependence_printdefault_out(logfile):
    """A job with linear dependency and the default number of MOs printed
    (all occupieds and 5 virtuals).
    """

    assert logfile.data.metadata["package_version"] == "4.2.2"

    nbasis = 138
    nmo = 106
    assert logfile.data.nbasis == nbasis
    assert logfile.data.nmo == nmo
    assert len(logfile.data.mocoeffs) == 1
    assert logfile.data.mocoeffs[0].shape == (nmo, nbasis)
    assert logfile.data.mocoeffs[0].T[59, 15] == -0.28758
    assert numpy.isnan(logfile.data.mocoeffs[0].T[59, 16])


def testQChem_QChem4_2_dvb_gopt_unconverged_out(logfile):
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert logfile.data.metadata["package_version"] == "4.2.0"
    assert hasattr(logfile.data, "optdone") and not logfile.data.optdone


def testQChem_QChem4_2_dvb_sp_multipole_10_out(logfile):
    """Multipole moments up to the 10-th order.

    Since this example has various formats for the moment ranks, we can test
    the parser by making sure the first moment (pure X) is as expected.
    """
    assert logfile.data.metadata["package_version"] == "4.2.0"
    assert hasattr(logfile.data, "moments") and len(logfile.data.moments) == 11
    tol = 1.0e-6
    assert logfile.data.moments[1][0] < tol
    assert abs(logfile.data.moments[2][0] - -50.9647) < tol
    assert abs(logfile.data.moments[3][0] - 0.0007) < tol
    assert abs(logfile.data.moments[4][0] - -1811.1540) < tol
    assert abs(logfile.data.moments[5][0] - 0.0159) < tol
    assert abs(logfile.data.moments[6][0] - -57575.0744) < tol
    assert abs(logfile.data.moments[7][0] - 0.3915) < tol
    assert numpy.isnan(logfile.data.moments[8][0])
    assert abs(logfile.data.moments[9][0] - 10.1638) < tol
    assert numpy.isnan(logfile.data.moments[10][0])


def testQChem_QChem4_2_MoOCl4_sp_noprint_builtin_mixed_all_Cl_out(logfile):
    """ECP on all Cl atoms, but iprint is off, so coreelectrons must be
    guessed.
    """
    assert logfile.data.metadata["package_version"] == "4.2.0"
    assert logfile.data.charge == -2
    assert logfile.data.mult == 1
    assert hasattr(logfile.data, "coreelectrons")
    coreelectrons = numpy.array([0, 0, 10, 10, 10, 10], dtype=int)
    assert numpy.all(coreelectrons == logfile.data.coreelectrons)


def testQChem_QChem4_2_MoOCl4_sp_noprint_builtin_mixed_both_out(logfile):
    """ECP on Mo and all Cl atoms, but iprint is off, so coreelectrons
    can't be guessed.

    Uses `ecp = gen`.
    """
    assert logfile.data.metadata["package_version"] == "4.2.0"
    assert logfile.data.charge == -2
    assert logfile.data.mult == 1
    assert not hasattr(logfile.data, "coreelectrons")


def testQChem_QChem4_2_MoOCl4_sp_noprint_builtin_mixed_single_Mo_out(logfile):
    """ECP on Mo, but iprint is off, so coreelectrons must be guessed."""
    assert logfile.data.metadata["package_version"] == "4.2.0"
    assert logfile.data.charge == -2
    assert logfile.data.mult == 1
    assert hasattr(logfile.data, "coreelectrons")
    coreelectrons = numpy.array([28, 0, 0, 0, 0, 0], dtype=int)
    assert numpy.all(coreelectrons == logfile.data.coreelectrons)


def testQChem_QChem4_2_MoOCl4_sp_noprint_builtin_out(logfile):
    """ECP on Mo and all Cl atoms, but iprint is off, so coreelectrons
    can't be guessed.

    Uses `ecp = <builtin>`.
    """
    assert logfile.data.metadata["package_version"] == "4.2.0"
    assert logfile.data.charge == -2
    assert logfile.data.mult == 1
    assert not hasattr(logfile.data, "coreelectrons")


def testQChem_QChem4_2_MoOCl4_sp_noprint_user_Mo_builtin_all_Cl_out(logfile):
    """ECP on Mo and all Cl atoms, but iprint is off; the coreelectrons
    count is given for Mo, and Cl can be guessed.
    """
    assert logfile.data.metadata["package_version"] == "4.2.0"
    assert logfile.data.charge == -2
    assert logfile.data.mult == 1
    assert hasattr(logfile.data, "coreelectrons")
    coreelectrons = numpy.array([28, 0, 10, 10, 10, 10], dtype=int)
    assert numpy.all(coreelectrons == logfile.data.coreelectrons)


def testQChem_QChem4_2_MoOCl4_sp_print_builtin_mixed_single_Mo_single_Cl_out(logfile):
    """ECP on Mo and all Cl atoms; iprint is on, so coreelectrons can be
    calculated.

    This was intended to only have an ECP on a single Cl, but Q-Chem
    silently puts it on all.
    """
    assert logfile.data.metadata["package_version"] == "4.2.0"
    assert logfile.data.charge == -2
    assert logfile.data.mult == 1
    assert hasattr(logfile.data, "coreelectrons")
    coreelectrons = numpy.array([28, 0, 10, 10, 10, 10], dtype=int)
    assert numpy.all(coreelectrons == logfile.data.coreelectrons)


def testQChem_QChem4_2_print_frgm_false_opt_out(logfile):
    """Fragment calculation: geometry optimization.

    Fragment sections are not printed.
    """

    assert logfile.data.metadata["package_version"] == "4.3.0"

    assert logfile.data.charge == -1
    assert logfile.data.mult == 1

    assert len(logfile.data.scfenergies) == 11
    assert len(logfile.data.grads) == 11


def testQChem_QChem4_2_print_frgm_true_opt_out(logfile):
    """Fragment calculation: geometry optimization.

    Fragment sections are printed.
    """

    assert logfile.data.metadata["package_version"] == "4.3.0"

    assert logfile.data.charge == -1
    assert logfile.data.mult == 1

    assert len(logfile.data.scfenergies) == 11
    assert len(logfile.data.grads) == 11


def testQChem_QChem4_2_print_frgm_false_sp_out(logfile):
    """Fragment calculation: single point energy.

    Fragment sections are not printed.
    """

    assert logfile.data.metadata["package_version"] == "4.3.0"

    assert logfile.data.charge == -1
    assert logfile.data.mult == 1

    assert len(logfile.data.scfenergies) == 1


def testQChem_QChem4_2_print_frgm_true_sp_out(logfile):
    """Fragment calculation: single point energy.

    Fragment sections are printed.
    """

    assert logfile.data.metadata["package_version"] == "4.3.0"

    assert logfile.data.charge == -1
    assert logfile.data.mult == 1

    assert len(logfile.data.scfenergies) == 1


def testQChem_QChem4_2_print_frgm_true_sp_ccsdt_out(logfile):
    """Fragment calculation: single point energy, CCSD(T).

    Fragment sections are printed.
    """

    assert logfile.data.metadata["package_version"] == "4.3.0"

    assert len(logfile.data.mpenergies[0]) == 1
    assert len(logfile.data.ccenergies) == 1


def testQChem_QChem4_2_qchem_tddft_rpa_out(logfile):
    """An RPA/TD-DFT job.

    Here Q-Chem prints both the TDA and RPA results. These differ somewhat, since
    TDA allows only X vectors (occupied-virtual transitions) whereas RPA also
    allows Y vectors (virtual-occupied deexcitations), and the formatting in these
    two cases is subtly different (see cclib/cclib#154 for details).

    Currently cclib will store the second set of transitions (RPA), but this
    could change in the future if we support multistep jobs.
    """

    assert logfile.data.metadata["package_version"] == "4.2.0"

    assert len(logfile.data.etsecs) == 10
    assert len(logfile.data.etsecs[0]) == 13

    # Check a few vectors manually, since we know the output. X vectors are transitions
    # from occupied to virtual orbitals, whereas Y vectors the other way around, so cclib
    # should be switching the indices. Here is the corresponding fragment in the logfile:
    #     Excited state 1: excitation energy (eV) = 3.1318
    #     Total energy for state 1: -382.185270280389
    #     Multiplicity: Triplet
    #     Trans. Mom.: 0.0000 X 0.0000 Y 0.0000 Z
    #     Strength : 0.0000
    #     X: D( 12) --> V( 13) amplitude = 0.0162
    #     X: D( 28) --> V( 5) amplitude = 0.1039
    #     Y: D( 28) --> V( 5) amplitude = 0.0605
    assert logfile.data.etsecs[0][0] == [(11, 0), (47, 0), 0.0162]
    assert logfile.data.etsecs[0][1] == [(27, 0), (39, 0), 0.1039]
    assert logfile.data.etsecs[0][2] == [(39, 0), (27, 0), 0.0605]


def testQChem_QChem4_2_read_molecule_out(logfile):
    """A two-calculation output with the charge/multiplicity not specified
    in the user section."""

    assert logfile.data.metadata["package_version"] == "4.3.0"

    # These correspond to the second calculation.
    assert logfile.data.charge == 1
    assert logfile.data.mult == 2
    assert len(logfile.data.moenergies) == 2

    # However, we currently take data from both, since they aren't
    # exactly fragment calculations.
    assert len(logfile.data.scfenergies) == 2


def testQChem_QChem4_2_stopiter_qchem_out(logfile):
    """Check to ensure that an incomplete SCF is handled correctly."""
    assert logfile.data.metadata["legacy_package_version"] == "4.0.0.1"
    assert logfile.data.metadata["package_version"] == "4.0.0.1"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)
    assert len(logfile.data.scfvalues[0]) == 7


def testQChem_QChem4_3_R_propylene_oxide_force_ccsd_out(logfile):
    """Check to see that the CCSD gradient (not the HF gradient) is being
    parsed.
    """
    assert logfile.data.metadata["package_version"] == "4.3.0"
    assert hasattr(logfile.data, "grads")
    assert logfile.data.grads.shape == (1, logfile.data.natom, 3)
    # atom 9, y-coordinate.
    idx = (0, 8, 1)
    assert logfile.data.grads[idx] == 0.00584973


def testQChem_QChem4_3_R_propylene_oxide_force_hf_numerical_energies_out(logfile):
    """Check to see that the HF numerical gradient (from energies) is
    being parsed.
    """
    assert logfile.data.metadata["package_version"] == "4.3.0"
    # This isn't implemented yet.
    assert not hasattr(logfile.data, "grads")


def testQChem_QChem4_3_R_propylene_oxide_force_mp2_out(logfile):
    """Check to see that the MP2 gradient (not the HF gradient) is
    being parsed.
    """
    assert logfile.data.metadata["package_version"] == "4.3.0"
    assert hasattr(logfile.data, "grads")
    assert logfile.data.grads.shape == (1, logfile.data.natom, 3)
    # atom 9, y-coordinate.
    idx = (0, 8, 1)
    assert logfile.data.grads[idx] == 0.00436177


def testQChem_QChem4_3_R_propylene_oxide_force_rimp2_out(logfile):
    """Check to see that the RI-MP2 gradient (not the HF gradient) is
    being parsed.
    """
    assert logfile.data.metadata["package_version"] == "4.3.0"
    assert hasattr(logfile.data, "grads")
    assert logfile.data.grads.shape == (1, logfile.data.natom, 3)
    # atom 9, y-coordinate.
    idx = (0, 8, 1)
    assert logfile.data.grads[idx] == 0.00436172


def testQChem_QChem4_3_R_propylene_oxide_freq_ccsd_out(logfile):
    """Check to see that the CCSD (numerical) Hessian is being parsed."""
    assert logfile.data.metadata["package_version"] == "4.3.0"
    # The gradient of the initial geometry in a Hessian calculated
    # from finite difference of gradients should be the same as in a
    # force calculation.
    assert hasattr(logfile.data, "grads")
    ngrads = 1 + 6 * logfile.data.natom
    assert logfile.data.grads.shape == (ngrads, logfile.data.natom, 3)
    # atom 9, y-coordinate.
    idx = (0, 8, 1)
    assert logfile.data.grads[idx] == 0.00584973

    assert hasattr(logfile.data, "hessian")
    assert logfile.data.hessian.shape == (3 * logfile.data.natom, 3 * logfile.data.natom)
    # atom 4, x-coordinate.
    idx = (9, 9)
    assert logfile.data.hessian[idx] == 0.3561243


def testQChem_QChem4_3_R_propylene_oxide_freq_hf_numerical_gradients_out(logfile):
    """Check to see that the HF Hessian (from gradients) is being parsed."""
    assert logfile.data.metadata["package_version"] == "4.3.0"
    # This isn't implemented yet.
    assert not hasattr(logfile.data, "freq")


def testQChem_QChem4_3_R_propylene_oxide_freq_mp2_out(logfile):
    """Check to see that the MP2 (numerical) Hessian is being parsed."""
    assert logfile.data.metadata["package_version"] == "4.3.0"
    # The gradient of the initial geometry in a Hessian calculated
    # from finite difference of gradients should be the same as in a
    # force calculation.
    assert hasattr(logfile.data, "grads")
    ngrads = 1 + 6 * logfile.data.natom
    assert logfile.data.grads.shape == (ngrads, logfile.data.natom, 3)
    # atom 9, y-coordinate.
    idx = (0, 8, 1)
    assert logfile.data.grads[idx] == 0.00436177

    assert hasattr(logfile.data, "hessian")
    assert logfile.data.hessian.shape == (3 * logfile.data.natom, 3 * logfile.data.natom)
    # atom 4, x-coordinate.
    idx = (9, 9)
    assert logfile.data.hessian[idx] == 0.3520255


def testQChem_QChem4_3_R_propylene_oxide_freq_rimp2_out(logfile):
    """Check to see that the RI-MP2 (numerical) Hessian is being parsed."""
    assert logfile.data.metadata["package_version"] == "4.3.0"
    # The gradient of the initial geometry in a Hessian calculated
    # from finite difference of gradients should be the same as in a
    # force calculation.
    assert hasattr(logfile.data, "grads")
    ngrads = 1 + 6 * logfile.data.natom
    assert logfile.data.grads.shape == (ngrads, logfile.data.natom, 3)
    # atom 9, y-coordinate.
    idx = (0, 8, 1)
    # Well, not quite in this case...
    assert logfile.data.grads[idx] == 0.00436167

    assert hasattr(logfile.data, "hessian")
    assert logfile.data.hessian.shape == (3 * logfile.data.natom, 3 * logfile.data.natom)
    # atom 4, x-coordinate.
    idx = (9, 9)
    assert logfile.data.hessian[idx] == 0.3520538


def testQChem_QChem4_4_full_2_out(logfile):
    """The polarizability section may not be parsed due to something
    appearing just beforehand from a frequency-type calculation.
    """
    assert logfile.data.metadata["legacy_package_version"] == "4.4.2"
    assert logfile.data.metadata["package_version"] == "4.4.2"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)
    assert hasattr(logfile.data, "polarizabilities")


def testQChem_QChem4_4_srtlg_out(logfile):
    """Some lines in the MO coefficients require fixed-width parsing. See
    #349 and #381.
    """
    assert logfile.data.metadata["legacy_package_version"] == "4.4.0"
    assert logfile.data.metadata["package_version"] == "4.4.0"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)
    # There is a linear dependence problem.
    nbasis, nmo = 1129, 1115
    assert len(logfile.data.mocoeffs) == 2
    assert logfile.data.mocoeffs[0].shape == (nmo, nbasis)
    assert logfile.data.mocoeffs[1].shape == (nmo, nbasis)
    index_ao = 151 - 1
    indices_mo = [index_mo - 1 for index_mo in (493, 494, 495, 496, 497, 498)]
    # line 306371:
    #  151  C 7   s      -54.24935 -36.37903-102.67529  32.37428-150.40380-103.24478
    ref = numpy.asarray([-54.24935, -36.37903, -102.67529, 32.37428, -150.40380, -103.24478])
    res = logfile.data.mocoeffs[1][indices_mo, index_ao]
    numpy.testing.assert_allclose(ref, res, atol=1.0e-5, rtol=0.0)


def testQChem_QChem4_4_Trp_polar_ideriv0_out(logfile):
    """Ensure that the polarizability section is being parsed, but don't
    compare to reference results as 2nd-order finite difference can have
    large errors.
    """
    assert logfile.data.metadata["package_version"] == "4.4.2"
    assert hasattr(logfile.data, "polarizabilities")


def testQChem_QChem4_4_top_out(logfile):
    """This job has fewer MOs (7) than would normally be printed (15)."""
    assert logfile.data.metadata["package_version"] == "4.4.2"
    nbasis = 7
    nmo = 7
    assert logfile.data.nbasis == nbasis
    assert logfile.data.nmo == nmo
    assert len(logfile.data.mocoeffs) == 1
    assert logfile.data.mocoeffs[0].shape == (nmo, nbasis)
    assert logfile.data.mocoeffs[0].T[6, 5] == 0.8115082


def testQChem_QChem5_0_438_out(logfile):
    """This job has an ECP on Pt, replacing 60 of 78 electrons, and was
    showing the charge as 60.
    """
    assert logfile.data.metadata["legacy_package_version"] == "5.0.0"
    assert logfile.data.metadata["package_version"] == "5.0.0"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)
    assert logfile.data.charge == 0
    assert logfile.data.coreelectrons[0] == 60


def testQChem_QChem5_0_argon_out(logfile):
    """This job has unit specifications at the end of 'Total energy for
    state' lines.
    """
    assert logfile.data.metadata["legacy_package_version"] == "5.0.1"
    assert logfile.data.metadata["package_version"] == "5.0.1"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)
    nroots = 12
    assert len(logfile.data.etenergies) == nroots
    state_0_energy = -526.6323968555
    state_1_energy = -526.14663738
    assert logfile.data.scfenergies[0] == state_0_energy
    assert abs(logfile.data.etenergies[0] - (state_1_energy - state_0_energy)) < 1.0e-1


def testQChem_QChem5_0_Si_out(logfile):
    """
    This job includes MOs as a test for this version. This fist MO coefficient is checked to ensure they were parsed.
    """
    assert logfile.data.metadata["legacy_package_version"] == "5.0.2"
    assert logfile.data.metadata["package_version"] == "5.0.2"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)
    assert logfile.data.mocoeffs[0][0, 0] == 1.00042


def testQChem_QChem5_1_old_final_print_1_out(logfile):
    """This job has was run from a development version."""
    assert logfile.data.metadata["legacy_package_version"] == "5.1.0"
    assert logfile.data.metadata["package_version"] == "5.1.0dev+branches_libresponse-27553"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


def testQChem_QChem5_2_HSO2F_TS_out(logfile):
    """This file consists of three parts: an initial frequency calculation, a
    transition state search, then a final frequency calculation for TS
    confirmation.  The last set of vibrational analysis information should be
    saved, not the first set.
    """
    assert logfile.data.vibfreqs[0] == -1388.93


def testQChem_QChem5_3_ccman2_soc_cisd_out(logfile):
    """This file has its atomcoords in bohr, which need to be converted."""
    convfac = 0.5291772109
    assert logfile.data.atomcoords[0, 0, 2] == -0.24685 * convfac
    assert logfile.data.atomcoords[0, 1, 2] == 1.72795 * convfac


def testQChem_QChem5_3_ts_30_irc_out(logfile):
    """This is a compound frequency -> TS -> freq -> IRC -> opt job,
    originally meant for #1040.

    It incorrectly has (developer) version information appended to metadata
    for each new calculation section.
    """
    assert logfile.data.metadata["package_version"] == "5.3.2dev+trunk-35976"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)


# Turbomole


def testTurbomole_Turbomole7_2_dvb_gopt_b3_lyp_Gaussian__(logfile):
    assert logfile.data.metadata["legacy_package_version"] == "7.2"
    assert logfile.data.metadata["package_version"] == "7.2.r21471"
    assert isinstance(parse_version(logfile.data.metadata["package_version"]), Version)
    assert logfile.data.natom == 20


def testTurbomole_Turbomole7_5_mp2_opt__(logfile):
    assert len(logfile.data.scfenergies) == len(logfile.data.mpenergies)


# These regression tests are for logfiles that are not to be parsed
# for some reason, and the function should start with 'testnoparse'.


def testnoparseADF_ADF2004_01_mo_sp_adfout(filename):
    """This is an ADF file that has a different number of AO functions
    and SFO functions. Currently nbasis parses the SFO count. This will
    be discussed and resolved in the future (see issue #170), and can
    this to get rid of the error in the meantime.
    """
    pass


def testnoparseGaussian_Gaussian09_coeffs_log(filename):
    """This is a test for a Gaussian file with more than 999 basis functions.

    The log file is too big, so we are just including a section. Before
    parsing, we set some attributes of the parser so that it all goes smoothly.
    """

    parser = Gaussian(__filedir__ / filename, loglevel=logging.ERROR)
    parser.nmo = 5
    parser.nbasis = 1128

    data = parser.parse()
    assert data.mocoeffs[0].shape == (5, 1128)
    assert data.aonames[-1] == "Ga71_19D-2"
    assert data.aonames[0] == "Mn1_1S"


# When a unit test is removed or replaced by a newer version, we normally want
# the old logfile to become a regression, namely to run the unit test as part of
# the regression suite. To this end, add the logfile path to the dictionary
# below along with the appropriate unit test class to use, and the appropriate
# regression test function will be created automatically. If modifications
# are necessary due to developments in the unit test class, tweak it here
# and provide the modified version of the test class.


class ADFGeoOptTest_noscfvalues(ADFGeoOptTest):
    @pytest.mark.skip("Cannot parse scfvalues from this file.")
    def testgeovalues_scfvalues(self, data: "ccData") -> None:
        """SCF cycles were not printed here."""

    @pytest.mark.skip("Cannot parse scfvalues from this file.")
    def testscftargetdim(self, data: "ccData") -> None:
        """SCF cycles were not printed here."""

    @pytest.mark.skip("Cannot parse scfvalues from this file.")
    def testscfvaluetype(self, data: "ccData") -> None:
        """SCF cycles were not printed here."""


class ADFSPTest_noscfvalues(ADFSPTest):
    @pytest.mark.skip("Cannot parse scfvalues from this file.")
    def testscftargetdim(self, data: "ccData") -> None:
        """SCF cycles were not printed here."""

    @pytest.mark.skip("Cannot parse scfvalues from this file.")
    def testscfvaluetype(self, data: "ccData") -> None:
        """SCF cycles were not printed here."""

    @pytest.mark.skip("Cannot parse aooverlaps from this file.")
    def testaooverlaps(self, data: "ccData") -> None:
        """AO overlaps were not printed here."""


class ADFSPTest_nosyms(ADFSPTest, GenericSPTest):
    foverlap00 = 1.00000
    foverlap11 = 0.99999
    foverlap22 = 0.99999

    @pytest.mark.skip("Symmetry labels were not printed here")
    def testsymlabels(self, data: "ccData") -> None:
        """Symmetry labels were not printed here."""


class ADFSPTest_nosyms_noscfvalues(ADFSPTest_nosyms):
    @pytest.mark.skip("Cannot parse scfvalues from this file.")
    def testscftargetdim(self, data: "ccData") -> None:
        """SCF cycles were not printed here."""

    @pytest.mark.skip("Cannot parse scfvalues from this file.")
    def testscfvaluetype(self, data: "ccData") -> None:
        """SCF cycles were not printed here."""

    @pytest.mark.skip("Cannot parse aooverlaps from this file.")
    def testaooverlaps(self, data: "ccData") -> None:
        """AO overlaps were not printed here."""

    def testmetadata_symmetry_detected(self, data: "ccData") -> None:
        """Symmetry is completely turned off and not even detected."""
        assert data.metadata["symmetry_detected"] == "c1"

    def testmetadata_symmetry_used(self, data: "ccData") -> None:
        """Symmetry is completely turned off and not even detected."""
        assert data.metadata["symmetry_used"] == "c1"


class ADFSPTest_nosyms_valence(ADFSPTest_nosyms):
    def testlengthmoenergies(self, data: "ccData") -> None:
        """Only valence orbital energies were printed here."""
        assert len(data.moenergies[0]) == 45
        assert numpy.isnan(data.moenergies[0][0])


class ADFSPTest_nosyms_valence_noscfvalues(ADFSPTest_nosyms_valence):
    @pytest.mark.skip("Cannot parse scfvalues from this file.")
    def testscftargetdim(self, data: "ccData") -> None:
        """SCF cycles were not printed here."""

    @pytest.mark.skip("Cannot parse scfvalues from this file.")
    def testscfvaluetype(self, data: "ccData") -> None:
        """SCF cycles were not printed here."""

    @pytest.mark.skip("Cannot parse moenergies from this file.")
    def testfirstmoenergy(self, data: "ccData") -> None:
        """MO energies were not printed here."""

    @pytest.mark.skip("Cannot parse aooverlaps from this file.")
    def testaooverlaps(self, data: "ccData") -> None:
        """AO overlaps were not printed here."""

    def testmetadata_symmetry_detected(self, data: "ccData") -> None:
        """Symmetry is completely turned off and not even detected."""
        assert data.metadata["symmetry_detected"] == "c1"

    def testmetadata_symmetry_used(self, data: "ccData") -> None:
        """Symmetry is completely turned off and not even detected."""
        assert data.metadata["symmetry_used"] == "c1"


# DALTON #


class DALTONBigBasisTest_aug_cc_pCVQZ(GenericBigBasisTest):
    contractions = {6: 29}
    spherical = True


class DALTONSPTest_nosymmetry(GenericSPTest):
    def testsymlabels(self, data: "ccData") -> None:
        """Are all the symmetry labels either Ag/u or Bg/u?"""
        # A calculation without symmetry, meaning it belongs to the C1 point
        # group, only has the `A` irreducible representation.
        sumwronglabels = sum(x not in {"A"} for x in data.mosyms[0])
        assert sumwronglabels == 0

    def testmetadata_symmetry_detected(self, data: "ccData") -> None:
        """Does metadata have expected keys and values?"""
        assert data.metadata["symmetry_detected"] == "c1"

    def testmetadata_symmetry_used(self, data: "ccData") -> None:
        """Does metadata have expected keys and values?"""
        assert data.metadata["symmetry_used"] == "c1"


class DALTONHFSPTest_nosymmetry(DALTONSPTest_nosymmetry, GenericHFSPTest):
    pass


class DALTONTDTest_noetsecs(DALTONTDTest):
    @pytest.mark.skip("etsecs cannot be parsed from this file")
    def testsecs(self, data: "ccData") -> None:
        pass

    @pytest.mark.skip("etsecs cannot be parsed from this file")
    def testsecs_transition(self, data: "ccData") -> None:
        pass


# GAMESS #


class GAMESSUSSPunTest_charge0(GenericSPunTest):
    def testcharge_and_mult(self, data: "ccData") -> None:
        """The charge in the input was wrong."""
        assert data.charge == 0

    def testatomcharges_mulliken(self, data: "ccData") -> None:
        """The charge in the input was wrong."""
        charges = data.atomcharges["mulliken"]
        assert abs(sum(charges)) < 0.001

    @pytest.mark.skip("HOMOs were incorrect due to charge being wrong")
    def testhomos(self, data: "ccData") -> None:
        """HOMOs were incorrect due to charge being wrong."""


class GamessIRTest_old(GamessIRTest):
    entropy_places = 5
    freeenergy_places = 2


class GAMESSUSIRTest_ts(GenericIRimgTest):
    @pytest.mark.skip("This is a transition state with different intensities")
    def testirintens(self, data: "ccData") -> None:
        """This is a transition state with different intensities."""


class GAMESSUSCISTest_dets(GenericCISTest):
    nstates = 10

    @pytest.mark.skip("This gives unexpected coeficcients, also for current unit tests.")
    def testetsecsvalues(self, data: "ccData") -> None:
        """This gives unexpected coeficcients, also for current unit tests."""


class GAMESSSPTest_noaooverlaps(GenericSPTest):
    @pytest.mark.skip("Cannot parse aooverlaps from this file.")
    def testaooverlaps(self, data: "ccData") -> None:
        """aooverlaps were not printed here."""


# Gaussian #


class GaussianSPunTest_nomosyms(GaussianSPunTest):
    @pytest.mark.skip("Cannot parse mosyms from this file.")
    def testmosyms(self, data: "ccData") -> None:
        """mosyms were not printed here."""


class GaussianSPunTest_nonaturalorbitals(GaussianCISTest):
    @pytest.mark.skip("Cannot parse natrual orbitals from this file.")
    def testnocoeffs(self, data: "ccData") -> None:
        """natural orbitals were not printed here."""

    @pytest.mark.skip("Cannot parse natrual orbital occupation numbers from this file.")
    def testnooccnos(self, data: "ccData") -> None:
        """natural orbital occupation numbers were not printed here."""


class GaussianPolarTest(ReferencePolarTest):
    """Customized static polarizability unittest, meant for calculations
    with symmetry enabled.
    """

    # Reference values are from Q-Chem 4.2/trithiolane_freq.out, since
    # with symmetry enabled Q-Chem reorients molecules similarly to
    # Gaussian.
    isotropic = 66.0955766
    principal_components = [46.71020322, 75.50778705, 76.06873953]
    # Make the thresholds looser because these test jobs use symmetry,
    # and the polarizability is orientation dependent.
    isotropic_delta = 2.0
    principal_components_delta = 0.7


# Jaguar #


class JaguarIRTest_v42(JaguarIRTest):
    @pytest.mark.skip("Data file does not contain force constants")
    def testvibfconsts(self, data: "ccData") -> None:
        pass


class JaguarSPTest_noatomcharges(JaguarSPTest):
    """Atomic partial charges were not printed in old Jaguar unit tests."""

    @pytest.mark.skip("Cannot parse atomcharges from this file.")
    def testatomcharges(self, data: "ccData") -> None:
        """Are atomic charges consistent with natom?"""

    @pytest.mark.skip("Cannot parse atomcharges from this file.")
    def testatomcharges_mulliken(self, data: "ccData") -> None:
        """Do Mulliken atomic charges sum to zero?"""

    @pytest.mark.skip("Cannot parse atomcharges from this file.")
    def testatomcharges_lowdin(self, data: "ccData") -> None:
        """Do Lowdin atomic charges sum to zero?"""


class JaguarSPTest_6_31gss(JaguarSPTest_noatomcharges):
    """AO counts and some values are different in 6-31G** compared to STO-3G."""

    nbasisdict = {1: 5, 6: 15}
    scfenergy = -387.06414443
    moenergy = -10.20198
    overlap01 = 0.22

    def testmetadata_basis_set(self, data: "ccData") -> None:
        """This calculation did not use STO-3G for the basis set."""
        assert data.metadata["basis_set"].lower() == "6-31g**"


class JaguarSPTest_6_31gss_nomosyms(JaguarSPTest_6_31gss):
    @pytest.mark.skip("Cannot parse mosyms from this file.")
    def testsymlabels(self, data: "ccData") -> None:
        """mosyms were not printed here."""

    def testmetadata_symmetry_detected(self, data: "ccData") -> None:
        """This calculation has symmetry detected but disabled."""
        assert data.metadata["symmetry_detected"] == "c2h"

    def testmetadata_symmetry_used(self, data: "ccData") -> None:
        """This calculation has symmetry detected but disabled."""
        assert data.metadata["symmetry_used"] == "c1"


class JaguarSPunTest_nomosyms(JaguarSPunTest):
    @pytest.mark.skip("Cannot parse mosyms from this file.")
    def testmosyms(self, data: "ccData") -> None:
        """mosyms were not printed here."""


class JaguarSPunTest_nmo_all(JaguarSPunTest):
    def testmoenergies(self, data: "ccData") -> None:
        """Some tests printed all MO energies apparently."""
        assert len(data.moenergies[0]) == data.nmo


class JaguarSPunTest_nmo_all_nomosyms(JaguarSPunTest_nmo_all):
    @pytest.mark.skip("Cannot parse mosyms from this file.")
    def testmosyms(self, data: "ccData") -> None:
        """mosyms were not printed here."""


class JaguarGeoOptTest_nmo45(GenericGeoOptTest):
    def testlengthmoenergies(self, data: "ccData") -> None:
        """Without special options, Jaguar only print Homo+10 orbital energies."""
        assert len(data.moenergies[0]) == 45


class JaguarSPTest_nmo45(JaguarSPTest_noatomcharges):
    def testlengthmoenergies(self, data: "ccData") -> None:
        """Without special options, Jaguar only print Homo+10 orbital energies."""
        assert len(data.moenergies[0]) == 45

    @pytest.mark.skip("Cannot parse mos from this file.")
    def testfornoormo(self, data: "ccData") -> None:
        """mos were not printed here."""

    @pytest.mark.skip("Cannot parse scftargets from this file.")
    def testscftargets(self, data: "ccData") -> None:
        """scftargets were not parsed correctly here."""

    @pytest.mark.skip("Cannot parse atomcharges from this file.")
    def testatomcharges(self, data: "ccData") -> None:
        """atomcharges were not parsed correctly here."""

    @pytest.mark.skip("Cannot parse atombasis from this file.")
    def testatombasis(self, data: "ccData") -> None:
        """atombasis was not parsed correctly here."""


class JaguarGeoOptTest_6_31gss(GenericGeoOptTest):
    nbasisdict = {1: 5, 6: 15}
    scfenergy = -387.064207


class MolcasBigBasisTest_nogbasis(MolcasBigBasisTest):
    @pytest.mark.skip("gbasis was not printed in this output file")
    def testgbasis(self, data: "ccData") -> None:
        """gbasis was not parsed for this file"""

    @pytest.mark.skip("gbasis was not printed in this output file")
    def testnames(self, data: "ccData") -> None:
        """gbasis was not parsed for this file"""

    @pytest.mark.skip("gbasis was not printed in this output file")
    def testprimitives(self, data: "ccData") -> None:
        """gbasis was not parsed for this file"""

    @pytest.mark.skip("gbasis was not printed in this output file")
    def testsizeofbasis(self, data: "ccData") -> None:
        """gbasis was not parsed for this file"""


# Molpro #


class MolproBigBasisTest_cart(MolproBigBasisTest):
    spherical = False


# ORCA #


class OrcaRelaxedScanTest(GenericRelaxedScanTest):
    """Customized relaxed potential energy surface scan unittest"""

    @pytest.fixture
    def extra(self) -> int:
        """extra indices"""
        return 1


class OrcaROCIS40Test(OrcaROCISTest):
    """Customized test for ROCIS"""

    # In ORCA 4.0, an additional spectrum ("COMBINED ELECTRIC DIPOLE +
    # MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM (Origin Independent,
    # Length Representation)") was present that is not in ORCA 4.1.
    # Changed via 1085. VELOCITY DIPOLE MOMENTS are not parsed.
    n_spectra = 8


class OrcaSPTest_nohirshfeld(OrcaSPTest):
    """Versions pre-5.0 did not specify calculating Hirshfeld atomic charges."""

    @pytest.mark.skip("atomcharges['hirshfeld'] were not calculated")
    def testatomcharges_hirshfeld(self, data: "ccData") -> None:
        """Hirshfeld atomic charges were not calculated"""


class OrcaSPTest_nobasis(OrcaSPTest_nohirshfeld):
    """Versions pre-4.0 do not concretely print the basis set(s) used aside
    from repeating the input file.
    """

    def testmetadata_basis_set(self, data: "ccData") -> None:
        assert "basis_set" not in data.metadata


class OrcaSPunTest_charge0(GenericSPunTest):
    def testcharge_and_mult(self, data: "ccData") -> None:
        """The charge in the input was wrong."""
        assert data.charge == 0

    def testatomcharges_mulliken(self, data: "ccData") -> None:
        """The charge in the input was wrong."""
        charges = data.atomcharges["mulliken"]
        assert abs(sum(charges)) < 0.001

    @pytest.mark.skip("HOMOs were incorrect due to charge being wrong.")
    def testhomos(self, data: "ccData") -> None:
        """HOMOs were incorrect due to charge being wrong."""

    def testorbitals(self, data: "ccData") -> None:
        """Closed-shell calculation run as open-shell."""
        assert data.closed_shell


class OrcaTDDFTTest_pre5(OrcaTDDFTTest):
    symmetries = [
        "Triplet",
        "Triplet",
        "Triplet",
        "Triplet",
        "Triplet",
        "Singlet",
        "Singlet",
        "Singlet",
        "Singlet",
        "Singlet",
    ]


class OrcaTDDFTTest_pre1085(OrcaTDDFTTest_pre5):
    def testoscs(self, data: "ccData") -> None:
        """These values changed in the electric dipole osc strengths prior to Orca 4.0. See PR1085"""
        assert len(data.etoscs) == self.number
        assert abs(max(data.etoscs) - 0.94) < 0.2


class OrcaIRTest_pre4(OrcaIRTest):
    """Customized vibrational frequency unittest"""

    # ORCA has a bug in the intensities for version < 4.0
    max_IR_intensity = 215
    zpve = 0.1921

    entropy = 0.00012080325339594164
    enthalpy = -381.85224835
    freeenergy = -381.88826585

    enthalpy_places = 3
    entropy_places = 6
    freeenergy_places = 3


class OrcaIRTest_old(OrcaIRTest_pre4):
    """The frequency part of this calculation didn't finish, but went ahead and
    printed incomplete and incorrect results anyway.
    """

    zpve = 0.0200

    enthalpy_places = -1
    entropy_places = 2
    freeenergy_places = -1

    @pytest.mark.skip("These values were wrong due to wrong input coordinates.")
    def testfreqval(self, data: "ccData") -> None:
        """These values were wrong due to wrong input coordinates."""

    @pytest.mark.skip("These values were wrong due to wrong input coordinates.")
    def testirintens(self, data: "ccData") -> None:
        """These values were wrong due to wrong input coordinates."""


# PSI3 #


class Psi3SPTest(GenericHFSPTest):
    """Customized restricted single point HF/KS unittest"""

    # The final energy is also a bit higher here, I think due to the fact
    # that a SALC calculation is done instead of a full LCAO.
    scfenergy = -378.895220030994
    moenergy = -11.085851

    @pytest.mark.skip("atomcharges not implemented for Psi3")
    def testatomcharges_mulliken(self, data: "ccData") -> None:
        pass

    @pytest.mark.skip("atomcharges not implemented for Psi3")
    def testatomcharges_lowdin(self, data: "ccData") -> None:
        pass

    @pytest.mark.skip("atomcharges not implemented for Psi3")
    def testatomcharges_hirshfeld(self, data: "ccData") -> None:
        pass

    @pytest.mark.skip("atommasses not implemented yet")
    def testatommasses(self, data: "ccData") -> None:
        pass

    @pytest.mark.skip("atomcharges not implemented for Psi3")
    def testatomcharges(self, data: "ccData") -> None:
        pass

    @pytest.mark.skip("MO coefficients are printed separately for each SALC")
    def testfornoormo(self, data: "ccData") -> None:
        pass

    @pytest.mark.skip("MO coefficients are printed separately for each SALC")
    def testdimmocoeffs(self, data: "ccData") -> None:
        pass

    @pytest.mark.skip("Psi3 does not currently have the option to print the overlap matrix")
    def testaooverlaps(self, data: "ccData") -> None:
        pass

    @pytest.mark.skip("reading rotational constants is not implemented for Psi3")
    def testrotconsts(self, data: "ccData") -> None:
        pass

    @pytest.mark.skip("reading basis set is not implemented for Psi3")
    def testmetadata_basis_set(self, data: "ccData") -> None:
        pass

    @pytest.mark.skip("reading input file is not implemented for Psi3")
    def testmetadata_input_file(self, data: "ccData") -> None:
        pass

    @pytest.mark.skip("reading unformatted package version is not implemented for Psi3")
    def testmetadata_legacy_package_version(self, data: "ccData") -> None:
        pass

    @pytest.mark.skip("reading cpu/wall time is not implemented for Psi3")
    def testmetadata_times(self, data: "ccData") -> None:
        pass


# PSI4 #


class PsiSPTest_noatommasses(PsiSPTest):
    @pytest.mark.skip("atommasses were not printed in this file.")
    def testatommasses(self, data: "ccData") -> None:
        """These values are not present in this output file."""


class PsiHFSPTest_noatommasses(PsiHFSPTest):
    @pytest.mark.skip("atommasses were not printed in this file.")
    def testatommasses(self, data: "ccData") -> None:
        """These values are not present in this output file."""


class Psi4HFIRTest_v1(Psi4HFIRTest):
    max_force_constant = 7.981

    enthalpy_places = 2
    freeenergy_places = 2

    @pytest.mark.skip("not implemented in versions older than 1.7")
    def testirintens(self, data: "ccData") -> None:
        pass

    @pytest.mark.skip("not implemented in versions older than 1.2")
    def testvibrmasses(self, data: "ccData") -> None:
        pass
