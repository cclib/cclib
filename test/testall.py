import os
from cclib.parser import GAMESS,G03

os.chdir(os.path.join("..","data"))

testfiles = [G03(os.path.join("Gaussian","basicGaussian03","dvb_gopt.out")),
             GAMESS(os.path.join("GAMESS","basicPCGAMESS","dvb_gopt_a.out"))]

for testfile in testfiles:
	testfile.logger.setLevel(0)
	testfile.parse()

attribs = ['natom','homos','nbasis']
for attrib in attribs:
	print attrib,
	for testfile in testfiles:
		print testfile.__getattribute__(attrib),
	print

print "Energy of optimised molecule",
for testfile in testfiles:
	print testfile.scfenergies[-1],
print
print "Energy of HOMO",
for testfile in testfiles:
	print testfile.moenergies[0,testfile.homos[0]],
print
