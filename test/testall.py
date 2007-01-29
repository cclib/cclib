import unittest
import os
from cclib.parser import ADF, GAMESS, Gaussian, Jaguar, GAMESSUK

def getfile(parser,*location):
    """Returns a parsed logfile."""
    if parser.__name__ in ["GAMESS", "ADF", "Jaguar", "Gaussian"]:
        fullpath = ("..","data",parser.__name__) + location
    elif parser.__name__=="GAMESSUK":
        fullpath = ("..","data","GAMESS-UK") + location
    logfile = parser(os.path.join(*fullpath))
    logfile.logger.setLevel(0)
    logfile.parse()
    return logfile

def visualtests():
    """These are not formal tests -- but they should be eyeballed."""
    logfiles = [ getfile(Gaussian,"basicGaussian03","dvb_gopt.out"),
                 getfile(GAMESS,"basicPCGAMESS","dvb_gopt_a.out"),
                 getfile(GAMESS,"basicGAMESS-US","dvb_gopt_a.out"),
                 getfile(ADF,"basicADF2004.01","dvb_gopt.adfout"),
                 getfile(Jaguar,"basicJaguar4.2", "dvb_gopt.out"),
                 getfile(Jaguar,"basicJaguar6.5", "dvb_gopt.out"),
                 ]

    print "\n\nMO energies of optimised dvb"
    print "    ","".join(["%8s" % x for x in ['Gaussian','PCGAMESS','GAMESS-US','ADF','Jaguar']])
    print "HOMO", "  ".join(["%+2.4f" % x.moenergies[0][x.homos[0]] for x in logfiles])
    print "LUMO", "  ".join(["%+2.4f" % x.moenergies[0][x.homos[0]+1] for x in logfiles])
    print "H-L  ", "  ".join(["%2.4f" % (x.moenergies[0][x.homos[0]+1]-x.moenergies[0][x.homos[0]],) for x in logfiles])

def importName(modulename, name):
    """Import from a module whose name is determined at run-time.

    Taken from Python Cookbook 2nd ed O'Reilly Recipe 16.3
    """
    try:
        module = __import__(modulename, globals(), locals(), [name])
    except ImportError:
        return None
    return getattr(module, name)
    

if __name__=="__main__":
    perpackage = {}
    errors = []
    for module in [ "testGeoOpt", "testSP", "testSPun", "testBasis", "testvib", "testMP" ]:
        try:
            names = importName(module, "names") # i.e. from testGeoOpt import names
        except: # Parsing failed
            errors.append("ERROR: no tests run for %s as parsing failed." % module)
        else:
            tests = importName(module, "tests") # i.e. from testGeoOpt import tests
            for name,test in zip(names,tests):
                print "\n**** Testing %s (%s) ****" % (name, module)
                myunittest = unittest.makeSuite(test)
                a = unittest.TextTestRunner(verbosity=2).run(myunittest)
                l = perpackage.setdefault(name, [0, 0, 0])
                l[0] += a.testsRun
                l[1] += len(a.errors)
                l[2] += len(a.failures)

    print "\n\n********* SUMMARY PER PACKAGE ****************"
    names = perpackage.keys()
    names.sort()
    total = [0, 0, 0]
    print " "*14, "\t".join(["Total", "Passed", "Failed", "Errors"])
    for name in names:
        l = perpackage[name]
        print name.ljust(15), "%3d\t%3d\t%3d\t%3d" % (l[0], l[0]-l[1]-l[2], l[2], l[1])
        for i in range(3):
            total[i] += l[i]

    print "\n\n********* SUMMARY OF EVERYTHING **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total[0],total[0]-(total[1]+total[2]),total[2],total[1])

    if errors:
        print "\n".join(errors)

    print "\n\n*** Visual tests ***"
    visualtests()    
