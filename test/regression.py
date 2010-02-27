"""
A combined test framework for regression, ccopen and parsing which is
designed to make it easy to add new tests or datafiles.

To run the doctest, just use "python regression.py test".
"""

__revision__ = "$Revision$"

import os
import sys
import logging

from glob import glob

from cclib.parser import ccopen
from cclib.parser import Gaussian, GAMESS, GAMESSUK, Jaguar, ADF, Molpro

# Regression tests
def testGaussian_basicGaussian03_dvb_gopt_out(logfile):
    """Example regression test for Gaussian/basicGaussian03/dvb_gopt.out

    Note: the name of the test must match the full path to the datafile
    exactly, except that all periods, hyphens, path separators and
    parentheses are replaced by underscores.
    """
    assert len(logfile.homos)==1
    
def testGAMESS_GAMESS_US_paulo_h2o_mp2_zip(logfile):
    """Check that the new format for GAMESS MP2 is parsed"""
    assert hasattr(logfile, "mpenergies")
    assert len(logfile.mpenergies) == 1
    assert abs(logfile.mpenergies[0] + 2072.13) < 0.01
    
def testGAMESS_GAMESS_US_open_shell_ccsd_test_log_gz(logfile):
    """Parse ccenergies from open shell CCSD calculations"""
    assert hasattr(logfile, "ccenergies")
    assert len(logfile.ccenergies) == 1
    assert abs(logfile.ccenergies[0] + 3501.50) < 0.01
def testGAMESS_basicGAMESS_US_water_ccd_out(logfile):
    """Keep parsing ccenergies correctly"""
    assert len(logfile.ccenergies) == 1
    assert abs(logfile.ccenergies[0] + 2074.22) < 0.01
def testGAMESS_basicGAMESS_US_water_ccsd_out(logfile):
    """Keep parsing ccenergies correctly"""
    assert len(logfile.ccenergies) == 1
    assert abs(logfile.ccenergies[0] + 2074.24) < 0.01
def testGAMESS_basicGAMESS_US_water_ccsd_t__out(logfile):
    """Keep parsing ccenergies correctly"""
    assert len(logfile.ccenergies) == 1
    assert abs(logfile.ccenergies[0] + 2074.32) < 0.01

def testGaussian_Gaussian09_OPT_td_g09_zip(logfile):
    """
    Previously, couldn't find etrotats as G09 has different output than G03
    """
    assert len(logfile.etrotats) == 10
    assert logfile.etrotats[0] == -0.4568

def testGaussian_Gaussian09_OPT_td_zip(logfile):
    """
    Working fine - adding to ensure that CD is parsed correctly
    """
    assert len(logfile.etrotats) == 10
    assert logfile.etrotats[0] == -0.4568

def testGaussian_Gaussian09_534_out_zip(logfile):
    """
    Previously, caused etenergies parsing to fail
    """
    assert logfile.etsyms[0] == "Singlet-?Sym"
    assert logfile.etenergies[0] == 20920.55328

def testGaussian_Gaussian03_AM1_SP_out_gz(logfile):
    """
    Previously, caused scfvalue parsing to fail.
    """
    assert len(logfile.scfvalues[0])==12

def testGaussian_Gaussian03_mendes_zip(logfile):
    """
    Previously, failed to extract coreelectrons.
    """
    centers = [9, 10, 11, 27]
    for i, x in enumerate(logfile.coreelectrons):
        if i in centers:
            assert x == 10
        else:
            assert x == 0

def testGAMESS_basicPCGAMESS_dvb_td_out(logfile):
    """
    Previously, etoscs was not extracted for this TD DFT calculation.
    """
    assert len(logfile.etoscs) == 5
    
def testGaussian_Gaussian98_C_bigmult_log_gz(logfile):
    """
    This file failed first becuase it had a double digit multiplicity.
    Then it failed because it had no alpha virtual orbitals.
    """
    assert logfile.charge == -3
    assert logfile.mult == 10
    assert logfile.homos[0] == 8
    assert logfile.homos[1] == -1 # No occupied beta orbitals

def testGaussian_Gaussian03_chn1_log_gz(logfile):
    """
    This file failed to parse, due to the use of 'pop=regular'. We have
    decided that mocoeffs should not be defined for such calculations.
    """
    assert not hasattr(logfile, "mocoeffs")

def testGaussian_basicGaussian03_dvb_un_sp_out(logfile):
    """
    This file had no atomcoords at all at all, due to only having an Input
    Orientation section and no Standard Orientation.
    """
    assert len(logfile.atomnos) == 20
    assert logfile.atomcoords.shape == (1,20,3)

def testGaussian_Gaussian03_Mo4OSibdt2_opt_log_bz2(logfile):
    """
    This file had no atomcoords as it did not contain any
    "Input orientation" sections, only "Standard orientation" sections
    """
    assert hasattr(logfile,"atomcoords")

def testGaussian_basicGaussian03_dvb_raman_out(logfile):
    """
    Was extracting the "Depolar P" instead of the "Raman activity". Oops!
    """
    assert logfile.vibramans[1] - 2.6872 < 0.0001

def testADF_basicADF2004_01_dvb_sp_c_adfout(logfile):
    """
    Had homo[0] as 35, when it should be 34
    """
    assert logfile.homos[0] == 34
    
def testADF_ADF2004_01_Fe_ox3_final_out_gz(logfile):
    """
    Make sure HOMOS are correct
    """
    assert logfile.homos[0]==59 and logfile.homos[1]==54

def testGaussian_Gaussian98_water_zmatrix_nosym_log_gz(logfile):
    """
    This file had no atomcoords as it did not contain either an
    "Input orientation" or "Standard orientation section". As
    a result it failed to parse. Fixed in r400.

    This file is missing natom.
    """
    assert len(logfile.atomcoords)==1
    assert logfile.natom == 3

def testGaussian_Gaussian03_anthracene_log_gz(logfile):
    """
    This file exposed a bug in extracting the vibsyms.
    """
    assert len(logfile.vibsyms)==len(logfile.vibfreqs)

def testGaussian_Gaussian03_dvb_gopt_symmfollow_log_bz2(logfile):
    """
    Non-standard treatment of symmetry and thus non-standard printing of Standard orientation.
    The formatting of Standard orientation is a bit different, causing only the first coordinates to be read.
    """
    assert len(logfile.atomcoords)==len(logfile.geovalues)

def testGaussian_Gaussian03_cyclopropenyl_rhf_g03_cut_log_bz2(logfile):
    """
    Not using symmetry at all (option nosymm) means standard orientation is not printed.
    In this case inputcoords are copied by the parser, which up till now stored the last coordinates.
    """
    assert len(logfile.atomcoords)==len(logfile.geovalues)

# Edit the following variable definitions to add new parsers
# or new datafiles

data = os.path.join("..","data")
names = [ "Gaussian", "GAMESS", "ADF", "GAMESS UK", "Jaguar", "Molpro" ]
dummyfiles = [ Gaussian(""), GAMESS(""), ADF(""), GAMESSUK(""), Jaguar(""),
               Molpro("") ]

filenames = [glob(os.path.join(data, "Gaussian", "basicGaussian03", "*.out")) +  
             glob(os.path.join(data, "Gaussian", "basicGaussian03", "*.log")) +
             glob(os.path.join(data, "Gaussian", "basicGaussian09", "*.log")) +
             glob(os.path.join(data, "Gaussian", "Gaussian09", "*.zip")) +
             glob(os.path.join(data, "Gaussian", "Gaussian03", "*.bz2")) +
             glob(os.path.join(data, "Gaussian", "Gaussian03", "*.zip")) +
             glob(os.path.join(data, "Gaussian", "Gaussian03", "*.gz")) +
             glob(os.path.join(data, "Gaussian", "Gaussian98", "*.bz2")) +
             glob(os.path.join(data, "Gaussian", "Gaussian98", "*.gz")),
             
             glob(os.path.join(data, "GAMESS", "basicGAMESS-US", "*.out")) +
             glob(os.path.join(data, "GAMESS", "basicPCGAMESS", "*.out")) +
             glob(os.path.join(data, "GAMESS", "GAMESS-US", "*.out")) +
             glob(os.path.join(data, "GAMESS", "GAMESS-US", "*.bz2")) +
             glob(os.path.join(data, "GAMESS", "GAMESS-US", "*.gz")) +
             glob(os.path.join(data, "GAMESS", "GAMESS-US", "*.zip")) +
             glob(os.path.join(data, "GAMESS", "PCGAMESS", "*.*.bz2")) +
             glob(os.path.join(data, "GAMESS", "PCGAMESS", "*.*.gz")) +             
             glob(os.path.join(data, "GAMESS", "WinGAMESS", "*.gz")),
             
             glob(os.path.join(data, "ADF", "basicADF2006.01", "*.adfout")) +
             glob(os.path.join(data, "ADF", "ADF2004.01", "*.gz")) +
             glob(os.path.join(data, "ADF", "ADF2004.01", "*.bz2")) +
             glob(os.path.join(data, "ADF", "ADF2005.01", "*.zip")) +
             glob(os.path.join(data, "ADF", "ADF2006.01", "*.out")) +
             glob(os.path.join(data, "ADF", "ADF2006.01", "*.bz2")),
                          
             glob(os.path.join(data, "GAMESS-UK", "basicGAMESS-UK", "*.out")) +
             glob(os.path.join(data, "GAMESS-UK", "GAMESS-UK6.0", "*.out.gz")) +
             glob(os.path.join(data, "GAMESS-UK", "GAMESS-UK7.0", "*.out.gz")),
             
             glob(os.path.join(data, "Jaguar", "Jaguar4.2", "*.bz2")) +
             glob(os.path.join(data, "Jaguar", "Jaguar6.0", "*.bz2")) +
             glob(os.path.join(data, "Jaguar", "Jaguar6.5", "*.bz2")) +
             glob(os.path.join(data, "Jaguar", "Jaguar7.0", "*.bz2")) +
             glob(os.path.join(data, "Jaguar", "basicJaguar7.0", "*.out")),

             glob(os.path.join(data, "Molpro", "basicMolpro2006", "*.out")) +
             glob(os.path.join(data, "Molpro", "Molpro2006", "*.bz2")),
             ]

def normalisefilename(filename):
    """Replace all non-alphanumeric symbols by _

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

def flatten(seq):
    """Converts a list of lists [of lists] to a single flattened list."""
    # Taken from the web.
    res = []
    for item in seq:
        if (isinstance(item, (tuple, list))):
            res.extend(flatten(item))
        else:
            res.append(item)
    return res

def main():
    # Print a warning if you haven't downloaded all of the regression test files
    # or an error if all of the regression test files are not included in filenames[]
    regfile = open(os.path.join("..", "data", "regressionfiles.txt"), "r")
    regfilenames = [os.sep.join(x.strip().split("/")) for x in regfile.readlines()]
    regfile.close()
    missing = 0
    for x in regfilenames:
        if not os.path.isfile(os.path.join("..", "data", x)):
            missing += 1
        elif os.path.join("..", "data", x) not in flatten(filenames):
            print "\nERROR: The regression file %s is present, but not included in " \
                  "the 'filenames' variable.\n\nPlease add a new glob statement." % x
            sys.exit(1)       
    if missing > 0:
        print "\nWARNING: You are missing %d regression file(s).\n" \
              "         Run wget.sh in the ../data directory to update.\n" % missing
        try:
            raw_input("(Press ENTER to continue or CTRL+C to exit)")
        except KeyboardInterrupt:
            print "\n"
            sys.exit(0)
    
    failures = errors = total = 0
    for i in range(len(names)):
        print "Are the %s files ccopened and parsed correctly?" % names[i]
        for filename in filenames[i]:
            total += 1
            print "  %s..."  % filename,
            try:
                a  = ccopen(filename)
            except:
                errors += 1
                print "ccopen error"
            else:
                if type(a) == type(dummyfiles[i]):
                    try:
                        a.logger.setLevel(logging.ERROR)
                        data = a.parse()
                    except KeyboardInterrupt:
                        sys.exit(1)
                    except:
                        print "parse error"
                        errors += 1
                    else:    
                        fnname = "test" + normalisefilename("_".join(filename.split(os.sep)[2:]))
                        if fnname in globals(): # If there is a test that matches...
                            try:
                                eval(fnname)(data) # Run the test
                            except AssertionError:
                                print "test failed"
                                failures += 1
                            else:
                                print "parsed and tested"
                        else:
                            print "parsed"
                else:
                    print "ccopen failed"
                    failures += 1
        print

    print "Total: %d   Failed: %d  Errors: %d" % (total, failures, errors)

if __name__=="__main__":
    if len(sys.argv)==2 and sys.argv[1]=="test":
        import doctest
        doctest.testmod()
    else:
        main()
