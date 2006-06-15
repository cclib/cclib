import os
from cclib.parser import guesstype

# The following files broke the parser in the past
_brokeparser = ['Gaussian/wildGaussian98/oo-LAN.out',
                'GAMESS/wildGAMESS-US/bv6010sup1molecule_4.inp.cml_.log',
                'GAMESS/wildPCGAMESS/gamess.log' ]
for filename in _brokeparser:
    print "Attempting to parse %s." % filename
    newname = os.path.join("..","data",*filename.split("/")) # OS independent filenames
    t = guesstype(newname)
    t.parse()

print "Completed successfully"
