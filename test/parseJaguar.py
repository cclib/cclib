import os
from cclib.parser import Jaguar

os.chdir(os.path.join("..","data","Jaguar","basicJaguar"))

files = [ ["eg01","dvb_gopt.out"],
          ["eg02","dvb_sp.out"],
          ["eg03","dvb_ir.out"],
          ["eg06","dvb_un_sp.out"] ]

for f in files:
    t = Jaguar(os.path.join(f[0],f[1]))
    t.parse()
    if f[0]!="eg03":
        print t.scfvalues
