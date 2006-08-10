import os
from cclib.parser import ccopen

def makefilename(shortfilename):
    """Add ../data/ to the filenames in an OS-independent way."""
    return os.path.join("..","data",*shortfilename.split("/"))

# The following files broke the parser in the past
_brokeparser = ['Gaussian/Gaussian98/oo-LAN.out',
                'GAMESS/GAMESS-US/bv6010sup1molecule_4.inp.cml_.log',
                'GAMESS/PCGAMESS/aomixgamess.log',
                'ADF/ADF2004.01/Mo5Obdt2-opt-c2v.adfout']
for filename in _brokeparser:
    print "\n\nAttempting to parse %s." % filename
    newname = makefilename(filename)
    t = ccopen(newname)
    t.parse()

# The following file had no atomcoords as it did not contain any
# "Input orientation" sections, only "Standard orientation" sections
print "\n\nAttempting to assert that %s has atomcoords." % filename
filename = 'Gaussian/Gaussian03/Mo4OSibdt2-opt.log'
t = ccopen(makefilename(filename))
t.parse()
assert hasattr(t,"atomcoords")

print "\n\nCompleted successfully"
