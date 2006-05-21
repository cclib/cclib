import os
from cclib.parser import guesstype

os.chdir(os.path.join("..","data","Gaussian"))

os.chdir("basicGaussian03")

for filename in ["dvb_gopt.out"]:
    print guesstype(filename)
