# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""A combined test framework for regression, ccopen and parsing which is
designed to make it easy to add new tests or datafiles.

To run the doctest, just use "python regression.py test".
"""

__revision__ = "$Revision$"

import os
import sys
import inspect
import logging
import unittest

from glob import glob
from io import StringIO

from cclib.parser import ccopen
from cclib.parser import ADF, GAMESS, GAMESSUK, Gaussian, Jaguar, Molpro, ORCA

import testall
parsers = testall.parsers
test_modules = testall.test_modules

# Edit the following variable definitions to add new parsers or new datafile patterns.

data = os.path.join("..", "data")
dummyfiles = [eval(n)("") for n in parsers]

filenames = [glob(os.path.join(data, "ADF", "basicADF2006.01", "*.adfout")) +
             glob(os.path.join(data, "ADF", "ADF2004.01", "*.gz")) +
             glob(os.path.join(data, "ADF", "ADF2004.01", "*.bz2")) +
             glob(os.path.join(data, "ADF", "ADF2005.01", "*.zip")) +
             glob(os.path.join(data, "ADF", "ADF2006.01", "*.out")) +
             glob(os.path.join(data, "ADF", "ADF2006.01", "*.bz2")) +
             glob(os.path.join(data, "ADF", "ADF2009.01", "*.out")),
                          
             glob(os.path.join(data, "GAMESS", "basicGAMESS-US", "*.out")) +
             glob(os.path.join(data, "GAMESS", "basicPCGAMESS", "*.out")) +
             glob(os.path.join(data, "GAMESS", "GAMESS-US", "*.out")) +
             glob(os.path.join(data, "GAMESS", "GAMESS-US", "*.bz2")) +
             glob(os.path.join(data, "GAMESS", "GAMESS-US", "*.gz")) +
             glob(os.path.join(data, "GAMESS", "GAMESS-US", "*.zip")) +
             glob(os.path.join(data, "GAMESS", "PCGAMESS", "*.*.bz2")) +
             glob(os.path.join(data, "GAMESS", "PCGAMESS", "*.*.gz")) +             
             glob(os.path.join(data, "GAMESS", "WinGAMESS", "*.gz")),
             
             glob(os.path.join(data, "GAMESS-UK", "basicGAMESS-UK", "*.out")) +
             glob(os.path.join(data, "GAMESS-UK", "GAMESS-UK6.0", "*.out.gz")) +
             glob(os.path.join(data, "GAMESS-UK", "GAMESS-UK7.0", "*.out.gz")),
             
             glob(os.path.join(data, "Gaussian", "basicGaussian03", "*.out")) +  
             glob(os.path.join(data, "Gaussian", "basicGaussian03", "*.log")) +
             glob(os.path.join(data, "Gaussian", "basicGaussian09", "*.log")) +
             glob(os.path.join(data, "Gaussian", "Gaussian09", "*.gz")) +
             glob(os.path.join(data, "Gaussian", "Gaussian09", "*.zip")) +
             glob(os.path.join(data, "Gaussian", "Gaussian09", "*.bz2")) +
             glob(os.path.join(data, "Gaussian", "Gaussian03", "*.out")) +
             glob(os.path.join(data, "Gaussian", "Gaussian03", "*.bz2")) +
             glob(os.path.join(data, "Gaussian", "Gaussian03", "*.zip")) +
             glob(os.path.join(data, "Gaussian", "Gaussian03", "*.gz")) +
             glob(os.path.join(data, "Gaussian", "Gaussian98", "*.bz2")) +
             glob(os.path.join(data, "Gaussian", "Gaussian98", "*.gz")),
             
             glob(os.path.join(data, "Jaguar", "Jaguar4.2", "*.bz2")) +
             glob(os.path.join(data, "Jaguar", "Jaguar6.0", "*.bz2")) +
             glob(os.path.join(data, "Jaguar", "Jaguar6.5", "*.bz2")) +
             glob(os.path.join(data, "Jaguar", "Jaguar7.0", "*.bz2")) +
             glob(os.path.join(data, "Jaguar", "basicJaguar7.0", "*.out")),

             glob(os.path.join(data, "Molpro", "basicMolpro2006", "*.out")) +
             glob(os.path.join(data, "Molpro", "Molpro2006", "*.bz2")),

             glob(os.path.join(data, "ORCA", "basicORCA2.6", "*.out")) +
             glob(os.path.join(data, "ORCA", "basicORCA2.8", "*.out")) +
             glob(os.path.join(data, "ORCA", "basicORCA2.9", "*.out")) +
             glob(os.path.join(data, "ORCA", "ORCA2.8", "*.out")) +
             glob(os.path.join(data, "ORCA", "ORCA2.8", "*.out.gz")) +
             glob(os.path.join(data, "ORCA", "ORCA2.9", "*.out.gz")) +
             glob(os.path.join(data, "ORCA", "ORCA2.9", "*.out.bz2"))
             ]


# The regression test functions defined below should be named according to the path
# of the logfile, with some characters changed according to normalisefilename().

def testADF_basicADF2004_01_dvb_sp_c_adfout(logfile):
    """Had homo[0] as 35, when it should be 34."""
    assert logfile.data.homos[0] == 34


def testADF_ADF2004_01_Fe_ox3_final_out_gz(logfile):
    """Make sure HOMOS are correct."""
    assert logfile.data.homos[0]==59 and logfile.data.homos[1]==54


def testGAMESS_basicGAMESS_US_water_ccd_out(logfile):
    """Keep parsing ccenergies correctly."""
    assert len(logfile.data.ccenergies) == 1
    assert abs(logfile.data.ccenergies[0] + 2074.22) < 0.01

def testGAMESS_basicGAMESS_US_water_ccsd_out(logfile):
    """Keep parsing ccenergies correctly"""
    assert len(logfile.data.ccenergies) == 1
    assert abs(logfile.data.ccenergies[0] + 2074.24) < 0.01

def testGAMESS_basicGAMESS_US_water_ccsd_t__out(logfile):
    """Keep parsing ccenergies correctly."""
    assert len(logfile.data.ccenergies) == 1
    assert abs(logfile.data.ccenergies[0] + 2074.32) < 0.01


def testGAMESS_basicPCGAMESS_dvb_td_out(logfile):
    """Previously, etoscs was not extracted for this TD DFT calculation."""
    assert len(logfile.data.etoscs) == 5


def testGAMESS_GAMESS_US_N2_UMP2_zip(logfile):
    """Check that the new format for GAMESS MP2 is parsed."""
    assert hasattr(logfile.data, "mpenergies")
    assert len(logfile.data.mpenergies) == 1
    assert abs(logfile.data.mpenergies[0] + 2975.97) < 0.01

def testGAMESS_GAMESS_US_N2_ROMP2_zip(logfile):
    """Check that the new format for GAMESS MP2 is parsed."""
    assert hasattr(logfile.data, "mpenergies")
    assert len(logfile.data.mpenergies) == 1
    assert abs(logfile.data.mpenergies[0] + 2975.97) < 0.01    

def testGAMESS_GAMESS_US_open_shell_ccsd_test_log_gz(logfile):
    """Parse ccenergies from open shell CCSD calculations."""
    assert hasattr(logfile.data, "ccenergies")
    assert len(logfile.data.ccenergies) == 1
    assert abs(logfile.data.ccenergies[0] + 3501.50) < 0.01

def testGAMESS_GAMESS_US_paulo_h2o_mp2_zip(logfile):
    """Check that the new format for GAMESS MP2 is parsed."""
    assert hasattr(logfile.data, "mpenergies")
    assert len(logfile.data.mpenergies) == 1
    assert abs(logfile.data.mpenergies[0] + 2072.13) < 0.01

def testGAMESS_WinGAMESS_dvb_td_trplet_2007_03_24_r1_out_gz(logfile):
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

def testGaussian_Gaussian09_Dahlgren_TS_zip(logfile):
    """Failed to parse ccenergies for a variety of reasons"""
    assert hasattr(logfile.data, "ccenergies")
    assert abs(logfile.data.ccenergies[0] - (-11819.96506609)) < 0.001

def testGaussian_basicGaussian03_water_ccd_log(logfile):
    """Ensure that ccenergies continues to be parsed correctly"""
    assert hasattr(logfile.data, "ccenergies")
    assert abs(logfile.data.ccenergies[0] - (-2041.32671705)) < 0.001
def testGaussian_basicGaussian03_water_ccsd_log(logfile):
    """Ensure that ccenergies continues to be parsed correctly"""
    assert hasattr(logfile.data, "ccenergies")
    assert abs(logfile.data.ccenergies[0] - (-2041.33503383)) < 0.001
def testGaussian_basicGaussian03_water_ccsd_t__log(logfile):
    """Ensure that ccenergies continues to be parsed correctly"""
    assert hasattr(logfile.data, "ccenergies")
    assert abs(logfile.data.ccenergies[0] - (-2041.3371232)) < 0.001

    ##### ADD TESTS FOR THE OTHERS ######

def testGaussian_Gaussian09_dvb_lowdin_log_gz(logfile):
    """Check if both Mulliken and Lowdin charges are parsed."""
    assert "mulliken" in logfile.data.atomcharges
    assert "lowdin" in logfile.data.atomcharges

def testGaussian_basicGaussian03_dvb_gopt_out(logfile):
    """Example regression test for Gaussian/basicGaussian03/dvb_gopt.out

    Note: the name of the test must match the full path to the datafile
    exactly, except that all periods, hyphens, path separators and
    parentheses, which should be2 replaced by underscores.
    """
    assert len(logfile.data.homos)==1

def testGaussian_basicGaussian09_dvb_raman_log(logfile):
    """Was not extracting vibdisps"""
    assert hasattr(logfile.data, "vibdisps")
    assert len(logfile.data.vibdisps) == 54

def testGaussian_basicGaussian03_dvb_raman_out(logfile):
    """Was extracting the "Depolar P" instead of the "Raman activity". Oops!"""
    assert logfile.data.vibramans[1] - 2.6872 < 0.0001
    assert hasattr(logfile.data, "vibdisps")
    assert len(logfile.data.vibdisps) == 54

def testGaussian_basicGaussian03_dvb_un_sp_out(logfile):
    """
    This file had no atomcoords at all at all, due to only having an
    Input Orientation section and no Standard Orientation.
    """
    assert len(logfile.data.atomnos) == 20
    assert logfile.data.atomcoords.shape == (1, 20, 3)


def testGaussian_basicGaussian09_dvb_gopt_log(logfile):
    """Check that the atomnos is being parsed correctly."""
    assert hasattr(logfile.data, "atomnos"), "Missing atomnos"
    assert len(logfile.data.atomnos) == logfile.data.natom == 20


def testGaussian_Gaussian98_C_bigmult_log_gz(logfile):
    """
    This file failed first becuase it had a double digit multiplicity.
    Then it failed because it had no alpha virtual orbitals.
    """
    assert logfile.data.charge == -3
    assert logfile.data.mult == 10
    assert logfile.data.homos[0] == 8
    assert logfile.data.homos[1] == -1 # No occupied beta orbitals

def testGaussian_Gaussian98_water_zmatrix_nosym_log_gz(logfile):
    """This file is missing natom.

    This file had no atomcoords as it did not contain either an
    "Input orientation" or "Standard orientation section".
    As a result it failed to parse. Fixed in r400.
    """
    assert len(logfile.data.atomcoords)==1
    assert logfile.data.natom == 3


def testGaussian_Gaussian03_AM1_SP_out_gz(logfile):
    """Previously, caused scfvalue parsing to fail."""
    assert len(logfile.data.scfvalues[0])==12

def testGaussian_Gaussian03_anthracene_log_gz(logfile):
    """This file exposed a bug in extracting the vibsyms."""
    assert len(logfile.data.vibsyms) == len(logfile.data.vibfreqs)

def testGaussian_Gaussian03_chn1_log_gz(logfile):
    """
    This file failed to parse, due to the use of 'pop=regular'.
    We have decided that mocoeffs should not be defined for such calculations.
    """
    assert not hasattr(logfile.data, "mocoeffs")

def testGaussian_Gaussian03_cyclopropenyl_rhf_g03_cut_log_bz2(logfile):
    """
    Not using symmetry at all (option nosymm) means standard orientation
    is not printed. In this case inputcoords are copied by the parser,
    which up till now stored the last coordinates.
    """
    assert len(logfile.data.atomcoords)==len(logfile.data.geovalues)

def testGaussian_Gaussian03_DCV4T_C60_start_zip(logfile):
    """This is a test for a very large Gaussian file with > 99 atoms.

    The log file is too big, so we are just including the start.
    Previously, parsing failed in the pseudopotential section.
    """
    assert len(logfile.data.coreelectrons) == 102
    assert logfile.data.coreelectrons[101] == 2

def testGaussian_Gaussian03_dvb_gopt_symmfollow_log_bz2(logfile):
    """Non-standard treatment of symmetry.

    In this case the Standard orientation is also printed non-standard,
    which caused only the first coordinates to be read previously.
    """
    assert len(logfile.data.atomcoords) == len(logfile.data.geovalues)

def testGaussian_Gaussian03_mendes_zip(logfile):
    """Previously, failed to extract coreelectrons."""
    centers = [9, 10, 11, 27]
    for i, x in enumerate(logfile.data.coreelectrons):
        if i in centers:
            assert x == 10
        else:
            assert x == 0

def testGaussian_Gaussian03_Mo4OSibdt2_opt_log_bz2(logfile):
    """
    This file had no atomcoords as it did not contain any
    "Input orientation" sections, only "Standard orientation".
    """
    assert hasattr(logfile.data, "atomcoords")

def testGaussian_Gaussian03_orbgs_log_bz2(logfile):
    """Check that the pseudopotential is being parsed correctly."""
    assert hasattr(logfile.data, "coreelectrons"), "Missing coreelectrons"
    assert logfile.data.coreelectrons[0] == 28
    assert logfile.data.coreelectrons[15] == 10
    assert logfile.data.coreelectrons[20] == 10
    assert logfile.data.coreelectrons[23] == 10


def testGaussian_Gaussian09_25DMF_HRANH_zip(logfile):
    """Check that the anharmonicities are being parsed correctly."""
    assert hasattr(logfile.data, "vibanharms"), "Missing vibanharms"
    anharms = logfile.data.vibanharms
    N = len(logfile.data.vibfreqs)
    assert 39 == N == anharms.shape[0] == anharms.shape[1]
    assert abs(anharms[0][0] + 43.341) < 0.01
    assert abs(anharms[N-1][N-1] + 36.481) < 0.01 
    
def testGaussian_Gaussian09_534_out_zip(logfile):
    """Previously, caused etenergies parsing to fail."""
    assert logfile.data.etsyms[0] == "Singlet-?Sym"
    assert logfile.data.etenergies[0] == 20920.55328

def testGaussian_Gaussian09_OPT_td_g09_zip(logfile):
    """Couldn't find etrotats as G09 has different output than G03."""
    assert len(logfile.data.etrotats) == 10
    assert logfile.data.etrotats[0] == -0.4568

def testGaussian_Gaussian09_OPT_td_zip(logfile):
    """Working fine - adding to ensure that CD is parsed correctly."""
    assert len(logfile.data.etrotats) == 10
    assert logfile.data.etrotats[0] == -0.4568

def testGaussian_Gaussian09_Ru2bpyen2_H2_freq3_log_bz2(logfile):
    """Here atomnos wans't added to the gaussian parser before."""
    assert len(logfile.data.atomnos) == 69


def testORCA_ORCA2_8_co_cosmo_out_gz(logfile):
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


def testORCA_ORCA2_9_job_out_gz(logfile):
    """First output file and request to parse atomic spin densities.

    Make sure that the sum of such densities is one in this case (or reasonaby close),
    but remember that this attribute is a dictionary, so we must iterate.
    """
    assert all([abs(sum(v)-1.0) < 0.0001 for k,v in logfile.data.atomspins.items()])


# These regression tests are for logfiles that are not to be parsed
# for some reason, and the function should start with 'testnoparse'.

def testnoparseGaussian_Gaussian09_coeffs_zip(filename):
    """This is a test for a Gaussian file with more than 999 basis functions.

    The log file is too big, so we are just including a section. Before
    parsing, we set some attributes of the parser so that it all goes smoothly.
    """

    d = Gaussian(filename)
    d.logger.setLevel(logging.ERROR)
    d.nmo = 5
    d.nbasis  = 1128
    
    logfile = d.parse()
    assert logfile.data.mocoeffs[0].shape == (5, 1128)
    assert logfile.data.aonames[-1] == "Ga71_19D-2"
    assert logfile.data.aonames[0] == "Mn1_1S"


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
    >>> for x in [ "Gaussian_Gaussian03_Mo4OSibdt2-opt.log" ]:
    ...     print regression.normalisefilename(x)
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


def make_regression_from_old_unittest(filename, module_name, test_name):
    """Return a regression test function from an old unit test logfile.

    When using this, take precaution to ensure that test_name is
    a valid unit test class within module_name.
    """

    def old_unit_test(logfile):
        test_class = getattr(__import__(module_name), test_name)
        test_class.logfile = logfile
        test_class.data = logfile.data
        devnull = open(os.devnull, 'w')
        return unittest.TextTestRunner(stream=devnull).run(unittest.makeSuite(test_class))

    return old_unit_test


def main(which=[], traceback=False):

    # Print a warning if you haven't downloaded all of the regression test files,
    # or an error if not all of the regression test files are included in filenames.
    regfile = open(os.path.join("..", "data", "regressionfiles.txt"), "r")
    regfilenames = [os.sep.join(x.strip().split("/")) for x in regfile.readlines()]
    regfile.close()

    missing = 0
    for x in regfilenames:
        if not os.path.isfile(os.path.join("..", "data", x)):
            missing += 1
        elif os.path.join("..", "data", x) not in flatten(filenames):
            print("\nERROR: The regression file %s is present, but not included in " \
                  "the 'filenames' variable.\n\nPlease add a new glob statement." % x)
            sys.exit(1)       

    if missing > 0:
        print("\nWARNING: You are missing %d regression file(s).\n" \
              "         Run wget.sh in the ../data directory to update.\n" % missing)
        try:
            input("(Press ENTER to continue or CTRL+C to exit)")
        except KeyboardInterrupt:
            print("\n")
            sys.exit(0)

    # When a unit test is removed or replaced by a newer version, the old logfile
    # typically becomes a regression, and normally we still want to run the unit test
    # within the regression suite. To this end, add the logfile location to a list
    # called 'old_tests' in the appropriate unit test class, in which case the
    # following code will find it and create the becessary regression test function.
    for mod in test_modules:
        mod_name = "test" + mod
        tests = inspect.getmembers(__import__(mod_name), inspect.isclass)
        tests = [tc for tc in tests if tc[0][-4:] == "Test"]
        for test_name, test_class in tests:
            for old in getattr(test_class, "old_tests", []):
                funcname = "test" + normalisefilename(old)
                func = make_regression_from_old_unittest(old, mod_name, test_name)
                globals()[funcname] = func

    failures = errors = total = 0
    for iname, name in enumerate(parsers):

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
            fname_norm = normalisefilename("_".join(fname.split(os.sep)[2:]))

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
                    if type(logfile) == type(dummyfiles[iname]):
                        try:
                            logfile.logger.setLevel(logging.ERROR)
                            logfile.data = logfile.parse()
                        except KeyboardInterrupt:
                            sys.exit(1)
                        except:
                            print("parse error")
                            errors += 1
                        else:
                            if test_this:
                                try:
                                    res = eval(funcname)(logfile)
                                    if res and len(res.failures) > 0:
                                        failures += len(res.failures)
                                        print("%i test(s) failed" % len(res.failures))
                                        if traceback:
                                            for f in res.failures:
                                                print("Failure for", f[0])
                                                print(f[1])
                                        continue
                                except AssertionError:
                                    print("test failed")
                                    failures += 1
                                else:
                                    print("parsed and tested")
                            else:
                                print("parsed")
                    else:
                        print("ccopen failed")
                        failures += 1
            else:
                try:
                    eval(funcname_noparse)(filename)
                except AssertionError:
                    print("test failed")
                    failures += 1
                except:
                    print("parse error")
                    errors += 1
                else:
                    print("test passed")
                
        print
            
    print("Total: %d   Failed: %d  Errors: %d" % (total, failures, errors))
    if not traceback and failures + errors > 0:
        print("\nFor more information on failures/errors, add 'traceback' as argument.")


if __name__=="__main__":

    # If 'test' is passed as the first argument, do a doctest on this module.
    # Otherwise, any arguments are used to limit the test to the packages/parsers
    #   passed as arguments. Not arguments implies all parsers.
    if len(sys.argv) == 2 and sys.argv[1] == "test":
        import doctest
        doctest.testmod()
    else:
        traceback = "traceback" in sys.argv or "tb" in sys.argv
        main(sys.argv[1:], traceback)
