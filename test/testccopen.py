import os
from cclib.parser import ccopen

os.chdir(os.path.join("..","data","Gaussian"))

os.chdir("basicGaussian03")

for filename in ["dvb_gopt.out"]:
    print ccopen(filename)
