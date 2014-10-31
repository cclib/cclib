# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""A combined test framework for regression, ccopen and parsing.

The design here is supposed to make it easy to add new tests or datafiles.
To run the doctest, just use "python regression.py test".
"""

from __future__ import print_function

import glob
import importlib
import inspect
import logging
import os
import sys
import traceback
import unittest

import numpy

from cclib.parser import ccopen
from cclib.parser import ADF, GAMESS, GAMESSUK, Gaussian, Jaguar, Molpro, NWChem, ORCA, Psi, QChem

import testall


# The following regression test functions were manually written, because they
# contain custom checks that were determined on a per-file basis. Care needs to be taken
# that the function name corresponds to the path of the logfile, with some characters
# changed according to normalisefilename().

# ADF #

def testADF_ADF2004_01_Fe_ox3_final_out(logfile):
    """Make sure HOMOS are correct."""
    assert logfile.data.homos[0]==59 and logfile.data.homos[1]==54

def testADF_ADF2013_01_dvb_gopt_b_unconverged_adfout(logfile):
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert hasattr(logfile.data, 'optdone') and not logfile.data.optdone

# Firefly #

def testGAMESS_Firefly8_0_dvb_gopt_a_unconverged_out(logfile):
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert hasattr(logfile.data, 'optdone') and not logfile.data.optdone

# GAMESS #

def testGAMESS_GAMESS_US2008_N2_UMP2_out(logfile):
    """Check that the new format for GAMESS MP2 is parsed."""
    assert hasattr(logfile.data, "mpenergies")
    assert len(logfile.data.mpenergies) == 1
    assert abs(logfile.data.mpenergies[0] + 2975.97) < 0.01

def testGAMESS_GAMESS_US2008_N2_ROMP2_out(logfile):
    """Check that the new format for GAMESS MP2 is parsed."""
    assert hasattr(logfile.data, "mpenergies")
    assert len(logfile.data.mpenergies) == 1
    assert abs(logfile.data.mpenergies[0] + 2975.97) < 0.01    

def testGAMESS_GAMESS_US2009_open_shell_ccsd_test_log(logfile):
    """Parse ccenergies from open shell CCSD calculations."""
    assert hasattr(logfile.data, "ccenergies")
    assert len(logfile.data.ccenergies) == 1
    assert abs(logfile.data.ccenergies[0] + 3501.50) < 0.01

def testGAMESS_GAMESS_US2009_paulo_h2o_mp2_out(logfile):
    """Check that the new format for GAMESS MP2 is parsed."""
    assert hasattr(logfile.data, "mpenergies")
    assert len(logfile.data.mpenergies) == 1
    assert abs(logfile.data.mpenergies[0] + 2072.13) < 0.01

def testGAMESS_GAMESS_US2012_dvb_gopt_a_unconverged_out(logfile):
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert hasattr(logfile.data, 'optdone') and not logfile.data.optdone

def testGAMESS_GAMESS_US2013_N_UHF_out(logfile):
    """An UHF job that has an LZ value analysis between the alpha and beta orbitals."""
    assert len(logfile.data.moenergies) == 2

def testGAMESS_WinGAMESS_dvb_td_trplet_2007_03_24_r1_out(logfile):
    """Do some basic checks for this old unit test that was failing.

    The unit tests are not run automatically on this old unit logfile,
    because we know the output has etsecs whose sum is way off.
    So, perform a subset of the basic assertions for GenericTDTesttrp.
    """
    number = 5
    assert len(logfile.data.etenergies) == number
    idx_lambdamax = [i for i, x in enumerate(logfile.data.etoscs) if x == max(logfile.data.etoscs)][0]
    assert abs(logfile.data.etenergies[idx_lambdamax] - 24500) < 100
    assert len(logfile.data.etoscs) == number
    assert abs(max(logfile.data.etoscs) - 0.0) < 0.01
    assert len(logfile.data.etsecs) == number

# GAMESS-UK #

def testGAMESS_UK_GAMESS_UK8_0_dvb_gopt_hf_unconverged_out(logfile):
    assert hasattr(logfile.data, 'optdone') and not logfile.data.optdone

# Gaussian #

def testGaussian_Gaussian98_C_bigmult_log(logfile):
    """
    This file failed first becuase it had a double digit multiplicity.
    Then it failed because it had no alpha virtual orbitals.
    """
    assert logfile.data.charge == -3
    assert logfile.data.mult == 10
    assert logfile.data.homos[0] == 8
    assert logfile.data.homos[1] == -1 # No occupied beta orbitals

def testGaussian_Gaussian98_test_Cu2_log(logfile):
    """An example of the number of basis set function changing."""
    assert logfile.data.nbasis == 38

def testGaussian_Gaussian98_test_H2_log(logfile):
    """
    The atomic charges from a natural population analysis were
    not parsed correctly, and they should be zero for dihydrogen.
    """
    assert logfile.data.atomcharges['natural'][0] == 0.0
    assert logfile.data.atomcharges['natural'][1] == 0.0

def testGaussian_Gaussian98_water_zmatrix_nosym_log(logfile):
    """This file is missing natom.

    This file had no atomcoords as it did not contain either an
    "Input orientation" or "Standard orientation section".
    As a result it failed to parse. Fixed in r400.
    """
    assert len(logfile.data.atomcoords)==1
    assert logfile.data.natom == 3

def testGaussian_Gaussian03_AM1_SP_out(logfile):
    """Previously, caused scfvalue parsing to fail."""
    assert len(logfile.data.scfvalues[0])==13

def testGaussian_Gaussian03_anthracene_log(logfile):
    """This file exposed a bug in extracting the vibsyms."""
    assert len(logfile.data.vibsyms) == len(logfile.data.vibfreqs)

def testGaussian_Gaussian03_borane_opt_log(logfile):
    """An example of changing molecular orbital count."""
    assert logfile.data.nmo == 609

def testGaussian_Gaussian03_chn1_log(logfile):
    """
    This file failed to parse, due to the use of 'pop=regular'.
    We have decided that mocoeffs should not be defined for such calculations.
    """
    assert not hasattr(logfile.data, "mocoeffs")

def testGaussian_Gaussian03_cyclopropenyl_rhf_g03_cut_log(logfile):
    """
    Not using symmetry at all (option nosymm) means standard orientation
    is not printed. In this case inputcoords are copied by the parser,
    which up till now stored the last coordinates.
    """
    assert len(logfile.data.atomcoords)==len(logfile.data.geovalues)

def testGaussian_Gaussian03_DCV4T_C60_log(logfile):
    """This is a test for a very large Gaussian file with > 99 atoms.

    The log file is too big, so we are just including the start.
    Previously, parsing failed in the pseudopotential section.
    """
    assert len(logfile.data.coreelectrons) == 102
    assert logfile.data.coreelectrons[101] == 2

def testGaussian_Gaussian03_dvb_gopt_symmfollow_log(logfile):
    """Non-standard treatment of symmetry.

    In this case the Standard orientation is also printed non-standard,
    which caused only the first coordinates to be read previously.
    """
    assert len(logfile.data.atomcoords) == len(logfile.data.geovalues)

def testGaussian_Gaussian03_mendes_out(logfile):
    """Previously, failed to extract coreelectrons."""
    centers = [9, 10, 11, 27]
    for i, x in enumerate(logfile.data.coreelectrons):
        if i in centers:
            assert x == 10
        else:
            assert x == 0

def testGaussian_Gaussian03_Mo4OSibdt2_opt_log(logfile):
    """
    This file had no atomcoords as it did not contain any
    "Input orientation" sections, only "Standard orientation".
    """
    assert hasattr(logfile.data, "atomcoords")

def testGaussian_Gaussian03_orbgs_log(logfile):
    """Check that the pseudopotential is being parsed correctly."""
    assert hasattr(logfile.data, "coreelectrons"), "Missing coreelectrons"
    assert logfile.data.coreelectrons[0] == 28
    assert logfile.data.coreelectrons[15] == 10
    assert logfile.data.coreelectrons[20] == 10
    assert logfile.data.coreelectrons[23] == 10

def testGaussian_Gaussian09_25DMF_HRANH_log(logfile):
    """Check that the anharmonicities are being parsed correctly."""
    assert hasattr(logfile.data, "vibanharms"), "Missing vibanharms"
    anharms = logfile.data.vibanharms
    N = len(logfile.data.vibfreqs)
    assert 39 == N == anharms.shape[0] == anharms.shape[1]
    assert abs(anharms[0][0] + 43.341) < 0.01
    assert abs(anharms[N-1][N-1] + 36.481) < 0.01 
    
def testGaussian_Gaussian09_534_out(logfile):
    """Previously, caused etenergies parsing to fail."""
    assert logfile.data.etsyms[0] == "Singlet-?Sym"
    assert abs(logfile.data.etenergies[0] - 20920.55328) < 1.0

def testGaussian_Gaussian09_dvb_gopt_unconverged_log(logfile):
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert hasattr(logfile.data, 'optdone') and not logfile.data.optdone

def testGaussian_Gaussian09_dvb_lowdin_log(logfile):
    """Check if both Mulliken and Lowdin charges are parsed."""
    assert "mulliken" in logfile.data.atomcharges
    assert "lowdin" in logfile.data.atomcharges

def testGaussian_Gaussian09_Dahlgren_TS_log(logfile):
    """Failed to parse ccenergies for a variety of reasons"""
    assert hasattr(logfile.data, "ccenergies")
    assert abs(logfile.data.ccenergies[0] - (-11819.96506609)) < 0.001

def testGaussian_Gaussian09_irc_point_log(logfile):
    """Failed to parse vibfreqs except for 10, 11"""
    assert hasattr(logfile.data, "vibfreqs")
    assert len(logfile.data.vibfreqs) == 11

def testGaussian_Gaussian09_OPT_td_g09_out(logfile):
    """Couldn't find etrotats as G09 has different output than G03."""
    assert len(logfile.data.etrotats) == 10
    assert logfile.data.etrotats[0] == -0.4568

def testGaussian_Gaussian09_OPT_td_out(logfile):
    """Working fine - adding to ensure that CD is parsed correctly."""
    assert len(logfile.data.etrotats) == 10
    assert logfile.data.etrotats[0] == -0.4568

def testGaussian_Gaussian09_Ru2bpyen2_H2_freq3_log(logfile):
    """Here atomnos wans't added to the gaussian parser before."""
    assert len(logfile.data.atomnos) == 69

# Jaguar #

# It would be good to have an unconverged geometry optimization so that
# we can test that optdone is set properly.
#def testJaguarX.X_dvb_gopt_unconverged:
#    assert hasattr(logfile.data, 'optdone') and not logfile.data.optdone

# Molpro #

def testMolpro_Molpro2012_dvb_gopt_unconverged_out(logfile):
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert hasattr(logfile.data, 'optdone') and not logfile.data.optdone

# NWChem #

def testNWChem_NWChem6_0_dvb_gopt_hf_unconverged_out(logfile):
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert hasattr(logfile.data, 'optdone') and not logfile.data.optdone

def testNWChem_NWChem6_0_dvb_sp_hf_moments_only_quadrupole_out(logfile):
    """Quadrupole moments are printed/parsed, but not lower moments (no shape)."""
    assert hasattr(logfile.data, 'moments') and len(logfile.data.moments) == 3
    assert len(logfile.data.moments[0]) == 3
    assert not logfile.data.moments[1].shape
    assert len(logfile.data.moments[2]) == 6

def testNWChem_NWChem6_0_dvb_sp_hf_moments_only_octupole_out(logfile):
    """Quadrupole moments are printed/parsed, but not lower moments (no shape)."""
    assert hasattr(logfile.data, 'moments') and len(logfile.data.moments) == 4
    assert len(logfile.data.moments[0]) == 3
    assert not logfile.data.moments[1].shape
    assert not logfile.data.moments[2].shape
    assert len(logfile.data.moments[3]) == 10

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

def testORCA_ORCA2_9_job_out(logfile):
    """First output file and request to parse atomic spin densities.

    Make sure that the sum of such densities is one in this case (or reasonaby close),
    but remember that this attribute is a dictionary, so we must iterate.
    """
    assert all([abs(sum(v)-1.0) < 0.0001 for k,v in logfile.data.atomspins.items()])

def testORCA_ORCA3_0_dvb_gopt_unconverged_out(logfile):
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert hasattr(logfile.data, 'optdone') and not logfile.data.optdone

# PSI #

def testPsi_Psi4_dvb_gopt_hf_unconverged_out(logfile):
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert hasattr(logfile.data, 'optdone') and not logfile.data.optdone

# Q-Chem #

def testQChem_QChem4_2_dvb_gopt_unconverged_out(logfile):
    """An unconverged geometry optimization to test for empty optdone (see #103 for details)."""
    assert hasattr(logfile.data, 'optdone') and not logfile.data.optdone

def testQChem_QChem4_2_dvb_sp_multipole_10_out(logfile):
    """Multipole moments up to the 10-th order.

    Since this example has various formats for the moment ranks, we can test
    the parser by making sure the first moment (pure X) is as expected.
    """
    assert hasattr(logfile.data, 'moments') and len(logfile.data.moments) == 11
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

def testQChem_QChem4_2_CH3___Na__out(logfile):
    """An unrestricted fragment job with BSSE correction.

    For now, we only keep the final block of MO energies. This corresponds to
    the section titled 'Done with counterpoise correction on fragments'.
    """
    assert len(logfile.data.atomnos) == 5
    assert len(logfile.data.scfenergies) == 7
    assert len(logfile.data.moenergies) == 2
    # There are 40 MOs in the supersystem.
    assert len(logfile.data.moenergies[0]) == 40
    assert len(logfile.data.moenergies[1]) == 40
    assert logfile.data.nmo == 40
    assert type(logfile.data.moenergies) == type([])
    assert type(logfile.data.moenergies[0]) == type(numpy.array([]))
    assert type(logfile.data.moenergies[1]) == type(numpy.array([]))

def testQChem_QChem4_2_CH4___Na__out(logfile):
    """A restricted fragment job with no BSSE correction.

    For now, we only keep the final block of MO energies. This corresponds to
    the section titled 'Done with SCF on isolated fragments'.
    """
    assert len(logfile.data.atomnos) == 6
    assert len(logfile.data.scfenergies) == 3
    assert len(logfile.data.moenergies) == 1
    # There are 42 MOs in the supersystem.
    assert len(logfile.data.moenergies[0]) == 42
    assert logfile.data.nmo == 42
    assert type(logfile.data.moenergies) == type([])
    assert type(logfile.data.moenergies[0]) == type(numpy.array([]))

# These regression tests are for logfiles that are not to be parsed
# for some reason, and the function should start with 'testnoparse'.

def testnoparseGaussian_Gaussian09_coeffs_log(filename):
    """This is a test for a Gaussian file with more than 999 basis functions.

    The log file is too big, so we are just including a section. Before
    parsing, we set some attributes of the parser so that it all goes smoothly.
    """

    parser = Gaussian(filename)
    parser.logger.setLevel(logging.ERROR)
    parser.nmo = 5
    parser.nbasis  = 1128
    
    data = parser.parse()
    assert data.mocoeffs[0].shape == (5, 1128)
    assert data.aonames[-1] == "Ga71_19D-2"
    assert data.aonames[0] == "Mn1_1S"


def flatten(seq):
    """Converts a list of lists [of lists] to a single flattened list.

    Taken from the web.
    """
    res = []
    for item in seq:
        if (isinstance(item, (tuple, list))):
            res.extend(flatten(item))
        else:
            res.append(item)
    return res

def normalisefilename(filename):
    """Replace all non-alphanumeric symbols by underscores.

    >>> import regression
    >>> for x in [ "Gaussian/Gaussian03/Mo4OSibdt2-opt.log" ]:
    ...     print(regression.normalisefilename(x))
    ...
    Gaussian_Gaussian03_Mo4OSibdt2_opt_log
    """
    ans = []
    for y in filename:
        x = y.lower()
        if (x >= 'a' and x <= 'z') or (x >= '0' and x <='9'):
            ans.append(y)
        else:
            ans.append("_")
    return "".join(ans)


# When a unit test is removed or replaced by a newer version, we normally want
# the old logfile to become a regression, namely to run the unit test as part of
# the regression suite. To this end, add the logfile path to the dictionary
# below along with the appropriate unit test class to use, and the appropriate
# regression test function will be created automatically. If modifications
# are necessary due to developments in the unit test class, tweak it here
# and provide the modified version of the test class.

# We're going to need to import all of the unit test modules.
test_modules = { m : importlib.import_module('test'+m) for m in testall.test_modules }

# Although there is probably a cleaner way to do this, making the unit class test names
# global makes reading the dictionary of old unit tests much easier, especially it
# will contain some classes defined here.
for m, module in test_modules.items():
    for name in dir(module):
        if name[-4:] == "Test":
            globals()[name] = getattr(module, name)

class ADFSPTest_nosyms(test_modules['SP'].ADFSPTest):
    def testsymlabels(self):
        """Symmetry labels were not printed here. PASS"""

class ADFSPTest_nosyms_valence(ADFSPTest_nosyms):
    def testlengthmoenergies(self):
        """Only valence orbital energies were printed here."""
        self.assertEquals(len(self.data.moenergies[0]), 45)
        self.assertEquals(self.data.moenergies[0][0], 99999.0)

class GAMESSUSSPunTest_charge0(GenericSPunTest):
    def testcharge_and_mult(self):
        """The charge in the input was wrong."""
        self.assertEquals(self.data.charge, 0)
    def testhomos(self):
        """HOMOs were incorrect due to charge being wrong. PASS"""

class GAMESSUSIRTest_ts(GenericIRTest):
    def testirintens(self):
        """This is a transition state with different intensities. PASS"""

class GAMESSUSCISTest_dets(GenericCISTest):
    nstates = 10
    def testetsecsvalues(self):
        """This gives unexpected coeficcients, also for current unit tests. PASS"""

class JaguarSPTest_6_31gss(GenericSPTest):
    b3lyp_energy = -10530
    nbasisdict = {1:5, 6:15}
    def testnbasis(self):
        """The AO count is larger in 6-31G** than STO-3G."""
        self.assertEquals(self.data.nbasis, 200)
    def testlengthmoenergies(self):
        """Some tests printed all MO energies apparently."""
        self.assertEquals(len(self.data.moenergies[0]), self.data.nmo)

class JaguarSPunTest_nmo_all(JaguarSPunTest):
    def testmoenergies(self):
        """Some tests printed all MO energies apparently."""
        self.assertEquals(len(self.data.moenergies[0]), self.data.nmo)

class JaguarGeoOptTest_nmo45(JaguarGeoOptTest):
    def testlengthmoenergies(self):
        """Without special options, Jaguar only print Homo+10 orbital energies."""
        self.assertEquals(len(self.data.moenergies[0]), 45)

class JaguarGeoOptTest_6_31gss(JaguarGeoOptTest):
    b3lyp_energy = -10530
    def testnbasis(self):
        """The AO count is larger in 6-31G** than STO-3G."""
        self.assertEquals(self.data.nbasis, 200)

class MolproBigBasisTest_cart(MolproBigBasisTest):
    spherical = False

class OrcaSPTest_3_21g(OrcaSPTest):
    b3lyp_energy = -10460
    def testatombasis(self):
        """The basis set here was 3-21G instead of STO-3G. PASS"""
    def testnbasis(self):
        """The basis set here was 3-21G instead of STO-3G."""
        self.assertEquals(self.data.nbasis, 110)

class OrcaGeoOptTest_3_21g(OrcaGeoOptTest):
    b3lyp_energy = -10460
    def testatombasis(self):
        """The basis set here was 3-21G instead of STO-3G. PASS"""
    def testnbasis(self):
        """The basis set here was 3-21G instead of STO-3G."""
        self.assertEquals(self.data.nbasis, 110)

class OrcaSPunTest_charge0(OrcaSPunTest):
    def testcharge_and_mult(self):
        """The charge in the input was wrong."""
        self.assertEquals(self.data.charge, 0)
    def testhomos(self):
        """HOMOs were incorrect due to charge being wrong. PASS"""

class OrcaTDDFTTest_error(OrcaTDDFTTest):
    def testoscs(self):
        """These values used to be less accurate, probably due to wrong coordinates."""
        self.assertEqual(len(self.data.etoscs), self.number)
        self.assertInside(max(self.data.etoscs), 1.0, 0.2)

class OrcaIRTest_old(OrcaIRTest):
    def testfreqval(self):
        """These values were wrong due to wrong input coordinates. PASS"""
    def testirintens(self):
        """These values were wrong due to wrong input coordinates. PASS"""

old_unittests = {

    "ADF/ADF2004.01/MoOCl4-sp.adfout":      ADFCoreTest,
    "ADF/ADF2004.01/dvb_gopt.adfout":       ADFGeoOptTest,
    "ADF/ADF2004.01/dvb_gopt_b.adfout":     ADFGeoOptTest,
    "ADF/ADF2004.01/dvb_sp.adfout":         ADFSPTest,
    "ADF/ADF2004.01/dvb_sp_b.adfout":       ADFSPTest,
    "ADF/ADF2004.01/dvb_sp_c.adfout":       ADFSPTest_nosyms_valence,
    "ADF/ADF2004.01/dvb_sp_d.adfout":       ADFSPTest_nosyms,
    "ADF/ADF2004.01/dvb_un_sp.adfout":      ADFSPunTest,
    "ADF/ADF2004.01/dvb_un_sp_c.adfout":    ADFSPunTest,
    "ADF/ADF2004.01/dvb_ir.adfout":         GenericIRTest,

    "ADF/ADF2006.01/dvb_gopt.adfout":       ADFGeoOptTest,

    "ADF/ADF2013.01/dvb_gopt_b_fullscf.adfout":       ADFGeoOptTest,
    "ADF/ADF2014.01/dvb_gopt_b_fullscf.out":       ADFGeoOptTest,

    "GAMESS/GAMESS-US2005/water_ccd_2005.06.27.r3.out":         GenericCCTest,
    "GAMESS/GAMESS-US2005/water_ccsd_2005.06.27.r3.out":        GenericCCTest,
    "GAMESS/GAMESS-US2005/water_ccsd(t)_2005.06.27.r3.out":     GenericCCTest,
    "GAMESS/GAMESS-US2005/water_cis_dets_2005.06.27.r3.out":    GAMESSUSCISTest_dets,
    "GAMESS/GAMESS-US2005/water_cis_saps_2005.06.27.r3.out":    GenericCISTest,
    "GAMESS/GAMESS-US2005/MoOCl4-sp_2005.06.27.r3.out":         GenericCoreTest,
    "GAMESS/GAMESS-US2005/water_mp2_2005.06.27.r3.out":         GenericMP2Test,

    "GAMESS/GAMESS-US2006/C_bigbasis_2006.02.22.r3.out":    GenericBigBasisTest,
    "GAMESS/GAMESS-US2006/dvb_gopt_a_2006.02.22.r2.out":    GenericGeoOptTest,
    "GAMESS/GAMESS-US2006/dvb_sp_2006.02.22.r2.out":        GenericSPTest,
    "GAMESS/GAMESS-US2006/dvb_un_sp_2006.02.22.r2.out":     GenericSPunTest,
    "GAMESS/GAMESS-US2006/dvb_ir.2006.02.22.r2.out":        GenericIRTest,
    "GAMESS/GAMESS-US2006/nh3_ts_ir.2006.2.22.r2.out":      GAMESSUSIRTest_ts,

    "GAMESS/GAMESS-US2010/dvb_gopt.log":    GenericGeoOptTest,
    "GAMESS/GAMESS-US2010/dvb_sp.log":      GenericSPTest,
    "GAMESS/GAMESS-US2010/dvb_sp_un.log":   GAMESSUSSPunTest_charge0,
    "GAMESS/GAMESS-US2010/dvb_td.log":      GAMESSUSTDDFTTest,
    "GAMESS/GAMESS-US2010/dvb_ir.log":      GenericIRTest,

    "GAMESS/PCGAMESS/C_bigbasis.out":       GenericBigBasisTest,
    "GAMESS/PCGAMESS/dvb_gopt_b.out":       GenericGeoOptTest,
    "GAMESS/PCGAMESS/dvb_ir.out":           FireflyIRTest,
    "GAMESS/PCGAMESS/dvb_raman.out":        GenericRamanTest,
    "GAMESS/PCGAMESS/dvb_sp.out":           GenericSPTest,
    "GAMESS/PCGAMESS/dvb_td.out":           GenericTDTest,
    "GAMESS/PCGAMESS/dvb_td_trplet.out":    GenericTDDFTtrpTest,
    "GAMESS/PCGAMESS/dvb_un_sp.out":        GenericSPunTest,
    "GAMESS/PCGAMESS/water_mp2.out":        GenericMP2Test,
    "GAMESS/PCGAMESS/water_mp3.out":        GenericMP3Test,
    "GAMESS/PCGAMESS/water_mp4.out":        GenericMP4SDQTest,
    "GAMESS/PCGAMESS/water_mp4_sdtq.out":   GenericMP4SDTQTest,

    "GAMESS/WinGAMESS/dvb_td_2007.03.24.r1.out":    GAMESSUSTDDFTTest,

    "Gaussian/Gaussian09/dvb_gopt_revA.02.out":         GaussianGeoOptTest,
    "Gaussian/Gaussian09/dvb_ir_revA.02.out":           GaussianIRTest,
    "Gaussian/Gaussian09/dvb_raman_revA.02.out":        GaussianRamanTest,
    "Gaussian/Gaussian09/dvb_scan_revA.02.log":         GaussianScanTest,
    "Gaussian/Gaussian09/dvb_sp_basis_b_gfprint.log":   GenericBasisTest,
    "Gaussian/Gaussian09/dvb_sp_basis_gfinput.log":     GenericBasisTest,
    "Gaussian/Gaussian09/dvb_sp_revA.02.out":           GaussianSPTest,
    "Gaussian/Gaussian09/dvb_td_revA.02.out":           GaussianTDDFTTest,
    "Gaussian/Gaussian09/dvb_un_sp_revA.02.log":        GaussianSPunTest,
    "Gaussian/Gaussian09/dvb_un_sp_b_revA.02.log":      GaussianSPunTest,

    "Jaguar/Jaguar4.2/dvb_gopt.out":    JaguarGeoOptTest_nmo45,
    "Jaguar/Jaguar4.2/dvb_gopt_b.out":  JaguarGeoOptTest,
    "Jaguar/Jaguar4.2/dvb_sp.out":      JaguarGeoOptTest_nmo45,
    "Jaguar/Jaguar4.2/dvb_sp_b.out":    JaguarGeoOptTest_nmo45,
    "Jaguar/Jaguar4.2/dvb_un_sp.out":   JaguarSPunTest_nmo_all,
    "Jaguar/Jaguar4.2/dvb_ir.out":      JaguarIRTest,

    "Jaguar/Jaguar6.0/dvb_gopt.out":    JaguarGeoOptTest_6_31gss,
    "Jaguar/Jaguar6.0/dvb_sp.out":      JaguarSPTest_6_31gss,
    "Jaguar/Jaguar6.0/dvb_un_sp.out" :  JaguarSPunTest_nmo_all,

    "Jaguar/Jaguar6.5/dvb_gopt.out":    JaguarGeoOptTest_nmo45,
    "Jaguar/Jaguar6.5/dvb_sp.out":      JaguarGeoOptTest_nmo45,
    "Jaguar/Jaguar6.5/dvb_un_sp.out":   JaguarSPunTest,
    "Jaguar/Jaguar6.5/dvb_ir.out":      JaguarIRTest,

    "Molpro/Molpro2006/C_bigbasis_cart.out":    MolproBigBasisTest_cart,

    "ORCA/ORCA2.6/dvb_gopt.out":    OrcaGeoOptTest_3_21g,
    "ORCA/ORCA2.6/dvb_sp.out":      OrcaSPTest_3_21g,
    "ORCA/ORCA2.6/dvb_td.out":      OrcaTDDFTTest_error,
    "ORCA/ORCA2.6/dvb_ir.out":      OrcaIRTest,

    "ORCA/ORCA2.8/dvb_gopt.out":    OrcaGeoOptTest,
    "ORCA/ORCA2.8/dvb_sp.out":      OrcaSPTest,
    "ORCA/ORCA2.8/dvb_sp_un.out":   OrcaSPunTest_charge0,
    "ORCA/ORCA2.8/dvb_td.out":      OrcaTDDFTTest,
    "ORCA/ORCA2.8/dvb_ir.out":      OrcaIRTest_old,
}

def make_regression_from_old_unittest(test_class):
    """Return a regression test function from an old unit test logfile."""

    def old_unit_test(logfile):
        test_class.logfile = logfile
        test_class.data = logfile.data
        devnull = open(os.devnull, 'w')
        return unittest.TextTestRunner(stream=devnull).run(unittest.makeSuite(test_class))

    return old_unit_test


def main(which=[], opt_traceback=False, opt_status=False):

    # It would be nice to fix the structure of this nested list,
    # because in its current form it is not amenable to tweaks.
    regdir = os.path.join("..", "data", "regression")
    programs = [os.path.join(regdir,testall.get_program_dir(p)) for p in testall.parsers]
    try:
        filenames = [[os.path.join(p,version,fn) for version in os.listdir(p) for fn in os.listdir(os.path.join(p,version))] for p in programs]
    except OSError as e:
        print(e)
        print("\nERROR: At least one program direcory is missing.")
        print("Run regression_download.sh in the ../data directory to update.")
        sys.exit(1)

    # This file should contain the paths to all regresssion test files we have gathered
    # over the years. It is not really necessary, since we can discover them on the disk,
    # but we keep it as a legacy and a way to track double check the regression tests.
    regfile = open(os.path.join("..", "data", "regressionfiles.txt"), "r")
    regfilenames = [os.sep.join(x.strip().split("/")) for x in regfile.readlines()]
    regfile.close()

    # We will want to print a warning if you haven't downloaded all of the regression
    # test files, or when, vice versa, not all of the regression test files found on disk
    # are included in filenames. However, gather that data here and print the warnings
    # at the end so that we test all available files and the messages are displayed
    # prominently at the end.
    missing_on_disk = []
    missing_in_list = []
    for fn in regfilenames:
        if not os.path.isfile(os.path.join("..", "data", "regression", fn)):
            missing_on_disk.append(fn)
    for fn in glob.glob(os.path.join('..', 'data', 'regression', '*', '*', '*')):
        if os.path.join(*fn.split(os.path.sep)[3:]) not in regfilenames:
            missing_in_list.append(fn)

    # Create the regression test functions from logfiles that were old unittests.
    for path, test_class in old_unittests.items():
        funcname = "test" + normalisefilename(path)
        func = make_regression_from_old_unittest(test_class)
        globals()[funcname] = func

    # Gather orphaned tests - functions starting with 'test' and not corresponding
    # to any regression file name.
    orphaned_tests = []
    for ip, parser in enumerate(testall.parsers):
        prefix = "test%s_%s" % (parser, parser)
        tests = [fn for fn in globals() if fn[:len(prefix)] == prefix]
        normalize = lambda fn: normalisefilename("_".join(fn.split(os.sep)[3:]))
        normalized = [normalize(fname) for fname in filenames[ip]]
        orphaned = [t for t in tests if t[4:] not in normalized]
        orphaned_tests.extend(orphaned)

    failures = errors = total = 0
    for iname, name in enumerate(testall.parsers):

        parser_class = eval(name)

        # Continue to next iteration if we are limiting the regression and the current
        #   name was not explicitely chosen (that is, passed as an argument).
        if len(which) > 0 and not name in which:
            continue;

        print("Are the %s files ccopened and parsed correctly?" % name)
        current_filenames = filenames[iname]
        current_filenames.sort()
        for fname in current_filenames:
            total += 1
            print("  %s..."  % fname, end=" ")

            # Check if there is a test (needs to be an appropriately named function).
            # If not, there can also be a test that does not assume the file is
            # correctly parsed (for fragments, for example), and these test need
            # to be additionaly prepended with 'testnoparse'.
            test_this = test_noparse = False
            fname_norm = normalisefilename("_".join(fname.split(os.sep)[3:]))

            funcname = "test" + fname_norm
            test_this = funcname in globals()

            funcname_noparse = "testnoparse" + fname_norm
            test_noparse = not test_this and funcname_noparse in globals()

            if not test_noparse:
                try:
                    logfile  = ccopen(fname)
                except:
                    errors += 1
                    print("ccopen error")
                else:
                    if type(logfile) == parser_class:
                        try:
                            logfile.logger.setLevel(logging.ERROR)
                            logfile.data = logfile.parse()
                        except KeyboardInterrupt:
                            sys.exit(1)
                        except Exception as e:
                            print("parse error")
                            errors += 1
                            if opt_traceback:
                                print(traceback.format_exc())
                        else:
                            if test_this:
                                try:
                                    res = eval(funcname)(logfile)
                                    if res and len(res.failures) > 0:
                                        failures += len(res.failures)
                                        print("%i test(s) failed" % len(res.failures))
                                        if opt_traceback:
                                            for f in res.failures:
                                                print("Failure for", f[0])
                                                print(f[1])
                                        continue
                                except AssertionError:
                                    print("test failed")
                                    failures += 1
                                    if opt_traceback:
                                        print(traceback.format_exc())
                                else:
                                    print("parsed and tested")
                            else:
                                print("parsed")
                    else:
                        print("ccopen failed")
                        failures += 1
            else:
                try:
                    eval(funcname_noparse)(fname)
                except AssertionError:
                    print("test failed")
                    failures += 1
                except:
                    print("parse error")
                    errors += 1
                else:
                    print("test passed")
                
        print()
            
    print("Total: %d   Failed: %d  Errors: %d" % (total, failures, errors))
    if not opt_traceback and failures + errors > 0:
        print("\nFor more information on failures/errors, add --traceback as an argument.")

    # Show these warnings at the end, so that they're easy to notice. Notice that the lists
    # were populated at the beginning of this function.
    if len(missing_on_disk) > 0:
        print("\nWARNING: You are missing %d regression file(s)." % len(missing_on_disk))
        print("Run regression_download.sh in the ../data directory to update.")
        print("Missing files:")
        print("\n".join(missing_on_disk))
    if len(missing_in_list) > 0:
        print("\nWARNING: The list in 'regressionfiles.txt' is missing %d file(s)." % len(missing_in_list))
        print("Add these files paths to the list and commit the change.")
        print("Missing files:")
        print("\n".join(missing_in_list))
    if len(orphaned_tests) > 0:
        print("\nWARNING: There are %d orphaned regression test functions." % len(orphaned_tests))
        print("Please make sure these function names correspond to regression files:")
        print("\n".join(orphaned_tests))

    if opt_status and failures+errors > 0:
        sys.exit(1)

if __name__=="__main__":

    # If 'test' is passed as the first argument, do a doctest on this module.
    # Otherwise, any arguments are used to limit the test to the packages/parsers
    # passed as arguments. No arguments implies all parsers.
    # In general, it would be best to replace this trickery by argparse magic.
    if len(sys.argv) == 2 and sys.argv[1] == "test":
        import doctest
        doctest.testmod()
    else:
        opt_traceback = "--traceback" in sys.argv
        opt_status = "--status" in sys.argv
        which = [arg for arg in sys.argv[1:] if not arg in ["--status", "--traceback"]]
        main(which, opt_traceback, opt_status)
