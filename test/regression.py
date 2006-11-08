"""
A combined test framework for regression, ccopen and parsing which is
designed to make it easy to add new tests or datafiles.

To run the doctest, just use "python regression.py test".
"""

import os
import sys
import logging

from glob import glob

from cclib.parser import ccopen
from cclib.parser import Gaussian, GAMESS, GAMESSUK, Jaguar, ADF

# Regression tests

def testGaussian_basicGaussian03_dvb_gopt_out(logfile):
    """Example regression test for Gaussian/basicGaussian03/dvb_gopt.out

    Note: the name of the test must match the full path to the datafile
    exactly, except that all periods are replaced by underscores, and path
    separators are also replaced by underscores.
    """
    assert len(logfile.homos)==1

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

# Edit the following variable definitions to add new parsers
# or new datafiles

data = os.path.join("..","data")
names = [ "Gaussian", "GAMESS", "ADF", "GAMESS UK", "Jaguar" ]
dummyfiles = [ Gaussian(""), GAMESS(""), ADF(""), GAMESSUK(""), Jaguar("") ]

filenames = [glob(os.path.join(data, "Gaussian", "basicGaussian03", "*.out")) +  
             glob(os.path.join(data, "Gaussian", "basicGaussian03", "*.log")) +
             glob(os.path.join(data, "Gaussian", "Gaussian03", "*.bz2")) +
             glob(os.path.join(data, "Gaussian", "Gaussian98", "*.bz2")) +
             glob(os.path.join(data, "Gaussian", "Gaussian98", "*.gz")),
             
             glob(os.path.join(data, "GAMESS", "basicGAMESS-US", "*.out")) +
             glob(os.path.join(data, "GAMESS", "basicPCGAMESS", "*.out")) +
             glob(os.path.join(data, "GAMESS", "GAMESS-US", "*.bz2")) +
             glob(os.path.join(data, "GAMESS", "PCGAMESS", "*.*.bz2")),
             
             glob(os.path.join(data, "ADF", "basicADF2004.01", "*.adfout")) +
             glob(os.path.join(data, "ADF", "ADF2004.01", "*.gz")) +
             glob(os.path.join(data, "ADF", "ADF2005.01", "*.zip")),
             
             glob(os.path.join(data, "GAMESS-UK", "basicGAMESS-UK", "*.out")) +
             glob(os.path.join(data, "GAMESS-UK", "GAMESS-UK6.0", "*.out.gz")) +
             glob(os.path.join(data, "GAMESS-UK", "GAMESS-UK7.0", "*.out.gz")),
             
             glob(os.path.join(data, "Jaguar", "basicJaguar4.2", "*.out")) +
             glob(os.path.join(data, "Jaguar", "basicJaguar6.0", "*.out")) +
             glob(os.path.join(data, "Jaguar", "basicJaguar6.5", "*.out")),
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
    """Taken from the web."""
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
                        a.parse()
                    except KeyboardInterrupt:
                        sys.exit(1)
                    except:
                        print "parse error"
                        errors += 1
                    else:    
                        fnname = "test" + normalisefilename("_".join(filename.split(os.sep)[2:]))
                        if fnname in globals(): # If there is a test that matches...
                            try:
                                eval(fnname)(a) # Run the test
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
