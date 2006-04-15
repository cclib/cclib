import os
from cclib.parser import Jaguar

os.chdir(os.path.join("..","data","Jaguar","basicJaguar"))

os.chdir("eg01")

for file in ["dvb_gopt.out"]:
	t = Jaguar(file)
	t.parse()

print t.moenergies[0,:]
print t.homos[0]
print t.moenergies[0,t.homos[0]]
