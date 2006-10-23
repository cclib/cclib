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

def testGaussian_Gaussian03_Mo4OSibdt2_opt_log(logfile):
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

def testADF_ADF2004_01_Fe_ox3_final_out(logfile):
    """
    Make sure HOMOS are correct
    """
    assert logfile.homos[0]==59 and logfile.homos[1]==54

# Edit the following variable definitions to add new parsers
# or new datafiles

data = os.path.join("..","data")
names = [ "Gaussian", "GAMESS", "ADF", "GAMESS UK", "Jaguar" ]
dummyfiles = [ Gaussian(""), GAMESS(""), ADF(""), GAMESSUK(""), Jaguar("") ]

filenames = [glob(os.path.join(data, "Gaussian", "basicGaussian03", "*.out")) +  
             glob(os.path.join(data, "Gaussian", "basicGaussian03", "*.log")) +
             glob(os.path.join(data, "Gaussian", "Gaussian03", "*.log")) +
             glob(os.path.join(data, "Gaussian", "Gaussian98", "*.out")),
             
             glob(os.path.join(data, "GAMESS", "basicGAMESS-US", "*.out")) +
             glob(os.path.join(data, "GAMESS", "basicPCGAMESS", "*.out")) +
             glob(os.path.join(data, "GAMESS", "GAMESS-US", "*.log")) +
             glob(os.path.join(data, "GAMESS", "PCGAMESS", "*.log")),
             
             glob(os.path.join(data, "ADF", "basicADF2004.01", "*.adfout")) +
             glob(os.path.join(data, "ADF", "ADF2004.01", "*out")),
             
             glob(os.path.join(data, "GAMESS-UK", "basicGAMESS-UK", "*.out")),
             
             glob(os.path.join(data, "Jaguar", "basicJaguar4.2", "*.out")) +
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

def main():
    print
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
