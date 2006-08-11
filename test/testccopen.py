import os
import logging

from glob import glob

from cclib.parser import ccopen
from cclib.parser import Gaussian, GAMESS, GAMESSUK, Jaguar, ADF

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
             glob(os.path.join(data, "ADF", "ADF2004.01", "*.adfout")),
             glob(os.path.join(data, "GAMESS-UK", "basicGAMESS-UK", "*.out")),
             glob(os.path.join(data, "Jaguar", "basicJaguar4.2", "*.out")) +
             glob(os.path.join(data, "Jaguar", "basicJaguar6.5", "*.out")),
             ]

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
        else:
            if type(a) == type(dummyfiles[i]):
                try:
                    a.logger.setLevel(logging.ERROR)
                    a.parse()
                    print "yes"
                except:
                    print "no"
                    errors += 1
            else:
                print "file type not guessed correctly"
                failures += 1
    print

print "Total: %d   Failed: %d  Errors: %d" % (total, failures, errors)
