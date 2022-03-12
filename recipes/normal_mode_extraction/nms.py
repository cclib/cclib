import numpy as np
import pybel
import openbabel as ob
from data import element_radii

mDynetoMet = 1.0e-5 * 1.0e-3 * 1.0e10
# boltzmann constant (in )
Kb = 1.38064852e-23
# meters to Angstrom
MtoA = 1.0e10

class nmsgenerator():
    # xyz = initial min structure
    # nmo = normal mode displacements
    # fcc = force constants
    # spc = atomic species list
    # T   = temperature of displacement

    try:
        et = ob.OBElementTable
    except:
        et = {}  # OpenBabel 3 doesn't use this class anymore

    def __init__(self,xyz,nmo,fcc,spc,T,minfc = 1.0E-3):
        self.xyz = xyz
        self.nmo = nmo
        self.labels = spc
        self.fcc = np.array([i if i > minfc else minfc for i in fcc])
        self.chg = np.array([self.__getCharge__(i) for i in spc])
        self.radii = [element_radii[i] for i in self.chg]
        self.T = T
        self.Na = xyz.shape[0]
        self.Nf = nmo.shape[0]
        self.e0 = self.__get_energy__(xyz)

    # returns atomic charge (number)
    def __getCharge__(self,t):
        if isinstance(t, (int, np.int32)):
            return t
        return self.et.GetAtomicNum(t)

    # Checks for small atomic distances
    def __check_atomic_distances__(self,rxyz):
        for i in range(0,self.Na):
            for j in range(i+1,self.Na):
                Rij = np.linalg.norm(rxyz[i]-rxyz[j])
                if Rij < 0.6 * (self.radii[i] * self.radii[j]):
                    return True
        return False

    # get energy
    def __get_energy__(self, rxyz):
        xyzstr = '{}\n\n'.format(len(rxyz))
        for l, x in zip(self.labels, rxyz):
            xyzstr += '{} {} {} {}\n'.format(l, *x)
        mol = pybel.readstring('xyz', xyzstr)
        ff = pybel._forcefields['uff']
        if not ff.Setup(mol.OBMol):
            return None
        return ff.Energy()

    # Generate a structure
    def __genrandomstruct__(self):
        rdt = np.random.random(self.Nf+1)
        rdt[0] = 0.0
        norm = np.random.random(1)[0]
        rdt = norm*np.sort(rdt)
        rdt[self.Nf] = norm

        oxyz = self.xyz.copy()

        for i in range(self.Nf):
            # convert force constants from mDyne to atomic
            Ki = mDynetoMet * self.fcc[i]
            ci = rdt[i+1]-rdt[i]
            Sn = -1.0 if np.random.binomial(1,0.5,1) else 1.0
            Ri = Sn * MtoA * np.sqrt((3.0 * ci * float(self.Na) * Kb * self.T)/Ki)
            oxyz = oxyz + (Ri * self.nmo[i]).reshape(-1, 3)
        return oxyz

    # Call this to return a random structure
    def get_random_structure(self, emax=100):
        gs = True
        it = 0
        while gs:
            rxyz = self.__genrandomstruct__()
            #gs = self.__check_atomic_distances__(rxyz)
            e1 = self.__get_energy__(rxyz)
            if not e1:
                continue
            gs = e1 - self.e0 - emax > 0
            it += 1
            #if it % 100 == 0:
            #    print "Failed after ", it
        return e1, rxyz
