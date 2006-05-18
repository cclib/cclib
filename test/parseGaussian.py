import os
from cclib.parser import Gaussian

os.chdir(os.path.join("..","data","Gaussian"))

os.chdir("basicGaussian03")

for file in ["dvb_gopt.out","dvb_sp.out","dvb_ir.out","dvb_raman.out",
             "dvb_un_sp.out"]:
	t = Gaussian(file)
	t.parse()


