import os
from cclib.parser import G03

t = G03(os.path.join("..","data","Gaussian","dvb_gopt.out"))
t.parse()

for x in ['scftargets','geotargets','scfvalues','geovalues']:
    print x
    t.__getattribute__(x)+1

