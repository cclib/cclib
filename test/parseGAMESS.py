import os
from cclib.parser import GAMESS

os.chdir(os.path.join("..","data","GAMESS"))

for file in ["ex.out","WinGAMESS.log","exam01.out"]:
	t = GAMESS(file)
	t.parse()

os.chdir("basicPCGAMESS")

for file in ["dvb_gopt_a.out","dvb_gopt_b.out","dvb_sp.out","dvb_ir.out","dvb_raman.out",
             "dvb_un_sp.out"]:
	t = GAMESS(file)
	t.parse()


