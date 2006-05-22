import os
from cclib.parser import ADF

os.chdir(os.path.join("..","data","ADF"))

os.chdir("basicADF2004.01")

for file in ["dvb_gopt.adfout"]:
    t = ADF(file)
    t.parse()
    print len(t.atomcoords)
