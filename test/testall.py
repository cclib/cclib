import unittest
import os
from cclib.parser import G03,GAMESS,ADF,Jaguar

def getfile(parser,*location):
    """Returns a parsed logfile."""
    if parser.__name__ in ['GAMESS','ADF','Jaguar']:
        fullpath = ("..","data",parser.__name__) + location
    elif parser.__name__=="G03":
        fullpath = ("..","data","Gaussian") + location
    logfile = parser(os.path.join(*fullpath))
    logfile.logger.setLevel(0)
    logfile.parse()
    return logfile

def visualtests():
    """These are not formal tests -- but they should be eyeballed."""
    logfiles = [ getfile(G03,"basicGaussian03","dvb_gopt.out"),
                 getfile(GAMESS,"basicPCGAMESS","dvb_gopt_a.out"),
                 getfile(GAMESS,"basicGAMESS-US","dvb_gopt_a.out"),
                 getfile(ADF,"basicADF2004.01","dvb_gopt.adfout"),
                 getfile(Jaguar,"basicJaguar","eg01","dvb_gopt.out")]

    print "\n\nMO energies of optimised dvb"
    print "    ","".join(["%8s" % x for x in ['Gaussian','PCGAMESS','GAMESS-US','ADF','Jaguar']])
    print "HOMO", "  ".join(["%+2.4f" % x.moenergies[0,x.homos[0]] for x in logfiles])
    print "LUMO", "  ".join(["%+2.4f" % x.moenergies[0,x.homos[0]+1] for x in logfiles])
    print "H-L  ", "  ".join(["%2.4f" % (x.moenergies[0,x.homos[0]+1]-x.moenergies[0,x.homos[0]],) for x in logfiles])

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
    total = errors = failures = 0
    for module in [ "testGeoOpt", "testSP", "testSPun" ]:
        names = importName(module, "names") # i.e. from testGeoOpt import names
        tests = importName(module, "tests") # i.e. from testGeoOpt import tests
        for name,test in zip(names,tests):
            print "\n**** Testing %s (%s) ****" % (name, module)
            myunittest = unittest.makeSuite(test)
            a = unittest.TextTestRunner(verbosity=2).run(myunittest)
            total += a.testsRun
            errors += len(a.errors)
            failures += len(a.failures)

    print "\n\n********* SUMMARY OF EVERYTHING **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)

    print "\n\n*** Visual tests ***"
    visualtests()
