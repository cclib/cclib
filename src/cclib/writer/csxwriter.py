#the generation of a CSX file from a quantum chemistry package output file
#using cclib as a parser

from . import filewriter

class CSX(filewriter.Writer):

    def __init__(self, ccdata, splitfiles=False,
            firstgeom=False, lastgeom=True, allgeom=False,
            *args, **kwargs):
        
    # Call the __init__ method of the superclass
        super(CSX, self).__init__(ccdata, *args, **kwargs)

        self.do_firstgeom = firstgeom
        self.do_lastgeom = lastgeom
        self.do_allgeom = allgeom

        self.generate_csx()

    def generate_csx(self):
        import os
        import sys
        import math
        #import glob
        #import openbabel
        import chemElement
        import csx1_api as api
        #from cclib.parser import ccopen
        #from cclib.writer import ccwrite


        #inputfile = sys.argv[1]
        #myFile = ccopen(inputfile)
        fileName, fileExt = os.path.splitext(self.jobfilename)
        #fileName = str(self.jobfilename)
        csxfile = open( fileName + '.csx', 'w')
        data = self.ccdata
        #cmlfile = open( fileName + '.cml', 'w')
        #xyz = ccwrite(data, 'cml', cmlfile)
        atomNum = data.natom
        if ('basisname' in data.metadata):
            basisName = data.metadata['basisname']
        else:
            basisName = 'none'
        molCharge = data.charge
        molMulti = data.mult
        wfnRestricted = False
        hasOrb = True if (hasattr(data, 'mocoeffs')) else False
        hasOrbSym = True if (hasattr(data, 'mosyms')) else False
        hasFreq = True if (hasattr(data, 'vibfreqs')) else False
        hasProp = True if (hasattr(data, 'moments')) else False
        hasNMR = True if (hasattr(data, 'nmriso')) else False
        hasElec = True if (hasattr(data, 'etoscs')) else False
        molSpin = data.metadata['spintype'] if ('spintype' in data.metadata) else 'RHF'
        if molMulti == 1 :
            wfnRestricted = True
        else:
            if molSpin == 'ROHF':
                wfnRestricted = True
            else:
                molSpin = 'UHF'
        calcType = data.metadata['theory']
        molEE = data.scfenergies[-1]
        #Wavefunction
        if hasOrb:
            if wfnRestricted :
                orbEne = data.moenergies
                orbEString = ' '.join(str(x) for x in orbEne[0])
                orbNum = data.nmo
                orbOcc = []
                if hasOrbSym:
                    orbSym = data.mosyms
                    orbSymString = ' '.join( x for x in orbSym[0])
                for iorb in range (orbNum):
                    if molSpin == "ROHF":
                        if iorb < int(data.homos[0]):
                            elecNum = 2
                        elif iorb == int(data.homos[0]):
                            elecNum = 1
                        else:
                            elecNum = 0
                    else:
                        elecNum = 0 if iorb > int(data.homos) else 2
                    orbOcc.append(elecNum)
                orbOccString = ' '.join(str(x) for x in sorted(orbOcc,reverse=True))
                wfn1 = api.waveFunctionType(orbitalCount=orbNum, \
                        orbitalOccupancies=orbOccString)
                if hasOrbSym:
                    wfn1.set_orbitalSymmetry(orbSymString)
                orbe1 = api.stringArrayType(unit='cs:eV')
                orbe1.set_valueOf_(orbEString)
                orbs1 = api.orbitalsType()
                for orbArray in data.mocoeffs:
                    for iorb in range(orbNum):
                        orbCaString = ' '.join(str(x) for x in orbArray[iorb])
                        orb1 = api.stringArrayType(id=iorb+1)
                        orb1.set_valueOf_(orbCaString)
                        orbs1.add_orbital(orb1)
                wfn1.set_orbitals(orbs1)
                wfn1.set_orbitalEnergies(orbe1)
            else:
                orbNum = data.nmo
                orbCaEne = data.moenergies[0][:]
                orbCaEString = ' '.join(str(x) for x in orbCaEne)
                orbCbEne = data.moenergies[1][:]
                orbCbEString = ' '.join(str(x) for x in orbCbEne)
                orbCaOcc = []
                orbCbOcc = []
                for iorb in range (orbNum):
                    elecCa = 1 if orbCaEne[iorb] < 0.0 else 0
                    orbCaOcc.append(elecCa)
                    elecCb = 1 if orbCbEne[iorb] < 0.0 else 0
                    orbCbOcc.append(elecCb)
                orbCaOccString = ' '.join(str(x) for x in sorted(orbCaOcc,reverse=True))
                orbCbOccString = ' '.join(str(x) for x in sorted(orbCbOcc,reverse=True))
                orbCaSym = data.mosyms[0][:]
                orbCaSymString = ' '.join( x for x in orbCaSym)
                orbCbSym = data.mosyms[1][:]
                orbCbSymString = ' '.join( x for x in orbCbSym)
                wfn1 = api.waveFunctionType(orbitalCount=orbNum)
                orbe1 = api.stringArrayType(unit='cs:eV')
                orbe1.set_valueOf_(orbCaEString)
                orbs1 = api.orbitalsType()
                alphaOrb = data.mocoeffs[0][:]
                for iorb in range(orbNum):
                    orbCaString = ' '.join(str(x) for x in alphaOrb[iorb])
                    orb1 = api.stringArrayType(id=iorb+1)
                    orb1.set_valueOf_(orbCaString)
                    orbs1.add_orbital(orb1)
                wfn1.set_alphaOrbitals(orbs1)
                wfn1.set_alphaOrbitalEnergies(orbe1)
                wfn1.set_alphaOrbitalOccupancies(orbCaOccString)
                wfn1.set_alphaOrbitalSymmetry(orbCaSymString)
                orbe2 = api.stringArrayType(unit='cs:eV')
                orbe2.set_valueOf_(orbCbEString)
                orbs2 = api.orbitalsType()
                betaOrb = data.mocoeffs[1][:]
                for iorb in range(orbNum):
                    orbCbString = ' '.join(str(x) for x in betaOrb[iorb])
                    orb2 = api.stringArrayType(id=iorb+1)
                    orb2.set_valueOf_(orbCbString)
                    orbs2.add_orbital(orb2)
                wfn1.set_betaOrbitals(orbs2)
                wfn1.set_betaOrbitalEnergies(orbe2)
                wfn1.set_betaOrbitalOccupancies(orbCbOccString)
                wfn1.set_betaOrbitalSymmetry(orbCbSymString)

        #vibrational frequency
        if hasFreq:
            molFreqNum = len(data.vibfreqs)
            frqString = ' '.join(str(x) for x in data.vibfreqs)
            intString = ' '.join(str(x) for x in data.vibirs)
            vib1 = api.vibAnalysisType(vibrationCount=molFreqNum)
            freq1 = api.stringArrayType(unit="cs:cm-1")
            freq1.set_valueOf_(frqString)
            vib1.set_frequencies(freq1)
            irint1 = api.stringArrayType()
            irint1.set_valueOf_(intString)
            vib1.set_irIntensities(irint1)
            norms1 = api.normalModesType()
            normMdString = []
            for ifrq in range(molFreqNum):
                normM = []
                for iatm in range(atomNum):
                    for ixyz in range(3):
                        normM.append(data.vibdisps[ifrq][iatm][ixyz])
                normMdString.append(' '.join(str(x) for x in normM))
                norm1 = api.normalModeType(id=ifrq+1)
                norm1.set_valueOf_(normMdString[ifrq])
                norms1.add_normalMode(norm1)
            vib1.set_normalModes(norms1)

        #dipole moments information
        if len(data.moments) == 1:
            hasProp = False
        if hasProp:
        #    au2db = 2.541766
            molDipoleX = data.moments[1][0]
            molDipoleY = data.moments[1][1]
            molDipoleZ = data.moments[1][2]
            molDipoleTot = math.sqrt(molDipoleX*molDipoleX+molDipoleY*molDipoleY+molDipoleZ*molDipoleZ)
            prop1 = api.propertiesType()
            sprop1 = api.propertyType(name='dipoleMomentX',unit='cs:debye',moleculeId='m1')
            sprop1.set_valueOf_(molDipoleX)
            sprop2 = api.propertyType(name='dipoleMomentY',unit='cs:debye',moleculeId='m1')
            sprop2.set_valueOf_(molDipoleY)
            sprop3 = api.propertyType(name='dipoleMomentZ',unit='cs:debye',moleculeId='m1')
            sprop3.set_valueOf_(molDipoleZ)
            sprop4 = api.propertyType(name='dipoleMomentAverage',unit='cs:debye')
            sprop4.set_valueOf_(molDipoleTot)
            prop1.add_systemProperty(sprop1)
            prop1.add_systemProperty(sprop2)
            prop1.add_systemProperty(sprop3)
            prop1.add_systemProperty(sprop4)
        #NMR chemical shielding information
        if hasNMR:
            prop2 = api.propertiesType()
            for iatm in range(atomNum):
                aprop1 = api.propertyType(atomId='a'+str(iatm+1), moleculeId='m1', \
                        name='nmrShieldingIsotropic',unit='cs:ppm')
                aprop1.set_valueOf_(data.nmriso[iatm])
                aprop2 = api.propertyType(atomId='a'+str(iatm+1), moleculeId='m1', \
                        name='nmrShieldingAnisotropic',unit='cs:ppm')
                aprop2.set_valueOf_(data.nmranis[iatm])
                prop2.add_atomProperty(aprop1)
                prop2.add_atomProperty(aprop2)

        #Electronic transition information
        if hasElec:
            transStr = ' '.join(str(x) for x in data.etenergies)
            oscilStr = ' '.join(str(x) for x in data.etoscs)
            elec1 = api.elecSpectraType()
            trans1 = api.stringArrayType(unit="cs:cm-1")
            trans1.set_valueOf_(transStr)
            oscil1 = api.stringArrayType()
            oscil1.set_valueOf_(oscilStr)
            elec1.set_electronicTransitions(trans1)
            elec1.set_oscillatorStrength(oscil1)

        #Start to generate CSX elements
        cs1 = api.csType(version='1.0')

        #molecular publication section
        mp1 = api.mpType(title='default title', \
                abstract='default abstract', \
                publisher='default publisher', \
                status='default status', \
                category=1, \
                visibility=0, \
                tags=data.metadata['package'], \
                key=1 )
        source1 = api.sourcePackageType(name=data.metadata['package'], \
                version=data.metadata['version'])
        mp1.set_sourcePackage(source1)
        ath1 = api.authorType(creator='Wang', \
                type_='cs:corresponding', \
                organization='default organization', \
                email='name@organization.com')
        mp1.add_author(ath1)
        cs1.set_molecularPublication(mp1)

        #molecular system section
        ms1 = api.msType(systemCharge=molCharge, \
               systemMultiplicity=molMulti)
        temp1 = api.dataWithUnitsType(unit='cs:kelvin')
        temp1.set_valueOf_(0.0)
        ms1.set_systemTemperature(temp1)
        mol1 = api.moleculeType(id='m1',atomCount=atomNum)
        if hasattr(data, "atomcharges"):
            atmCharge = data.atomcharges["mulliken"]
        else:
            atmCharge = [0]*atomNum
        #obmol1 = openbabel.OBMol()
        for iatm in range(atomNum):
            #   xCoord = float(data.atomcoords[iatm,0])
            xCoord = data.atomcoords[-1,iatm,0]
            yCoord = data.atomcoords[-1,iatm,1]
            zCoord = data.atomcoords[-1,iatm,2]
            xCoord1 = api.dataWithUnitsType(unit='cs:angstrom')
            xCoord1.set_valueOf_(xCoord)
            yCoord1 = api.dataWithUnitsType(unit='cs:angstrom')
            yCoord1.set_valueOf_(yCoord)
            zCoord1 = api.dataWithUnitsType(unit='cs:angstrom')
            zCoord1.set_valueOf_(zCoord)
            atomicNum = data.atomnos[iatm]
            atm = api.atomType(id='a'+str(iatm+1), elementSymbol=chemElement.z2elm[atomicNum], \
                    atomMass=chemElement.z2mass[atomicNum], \
                    xCoord3D=xCoord1, \
                    yCoord3D=yCoord1, \
                    zCoord3D=zCoord1, \
                    basisSet='cs:'+basisName, \
                    calculatedAtomCharge=atmCharge[iatm], \
                    formalAtomCharge=0)
            mol1.add_atom(atm)
        ms1.add_molecule(mol1)
        cs1.set_molecularSystem(ms1)

        #molCalculation section
        sd_wfn_method = ['HF', 'DFT', 'MP2', 'MP3', 'MP4', 'AM1', 'PM3', 'PM6']
        md_wfn_method = ['CCD', 'CCSD', 'CCSD-T']
        mc1 = api.mcType()
        qm1 = api.qmCalcType()
        srs1 = api.srsMethodType()
        if calcType in sd_wfn_method:
            sdm1 = api.srssdMethodType()
            #SCF
            if (calcType == 'HF'):
                scf1 = api.resultType(methodology='cs:normal',spinType='cs:'+molSpin, \
                        basisSet='bse:'+basisName)
                ene1 = api.energiesType(unit='cs:eV')
                ee_ene1 = api.energyType(type_='cs:totalPotential')
                ee_ene1.set_valueOf_(float(molEE))
                ene1.add_energy(ee_ene1)
                scf1.set_energies(ene1)
                if hasOrb:
                    scf1.set_waveFunction(wfn1)
                if hasProp:
                    scf1.set_properties(prop1)
                if hasNMR:
                    scf1.set_properties(prop2)
                if hasFreq:
                    scf1.set_vibrationalAnalysis(vib1)
                if hasElec:
                    scf1.set_electronicSpectra(elec1)
                sdm1.set_abinitioScf(scf1)
            #DFT
            elif (calcType == 'DFT'):
                dft1 = api.resultType(methodology='cs:normal',spinType='cs:'+molSpin, \
                        basisSet='bse:'+basisName, \
                        dftFunctional='cs:'+data.metadata['functional'])
                ene1 = api.energiesType(unit='cs:eV')
                ee_ene1 = api.energyType(type_='cs:totalPotential')
                ee_ene1.set_valueOf_(float(molEE))
                ene1.add_energy(ee_ene1)
                dft1.set_energies(ene1)
                if hasOrb:
                    dft1.set_waveFunction(wfn1)
                if hasProp:
                    dft1.set_properties(prop1)
                if hasNMR:
                    dft1.set_properties(prop2)
                if hasFreq:
                    dft1.set_vibrationalAnalysis(vib1)
                if hasElec:
                    dft1.set_electronicSpectra(elec1)
                sdm1.set_dft(dft1)
            #MP2
            elif (calcType == 'MP2'):
                mp21 = api.resultType(methodology='cs:normal',spinType='cs:'+molSpin, \
                        basisSet='bse:'+basisName)
                ene1 = api.energiesType(unit='cs:eV')
                ee_ene1 = api.energyType(type_='cs:totalPotential')
                ee_ene1.set_valueOf_(float(data.mpenergies[-1]))
                ce_ene1 = api.energyType(type_='cs:correlation')
                ce_ene1.set_valueOf_(float(data.mpenergies[-1])-float(molEE))
                ene1.add_energy(ee_ene1)
                ene1.add_energy(ce_ene1)
                mp21.set_energies(ene1)
                if hasOrb:
                    mp21.set_waveFunction(wfn1)
                if hasProp:
                    mp21.set_properties(prop1)
                if hasNMR:
                    mp21.set_properties(prop2)
                if hasFreq:
                    mp21.set_vibrationalAnalysis(vib1)
                sdm1.set_mp2(mp21)
            #Semiempirical methods
            elif (calcType == 'AM1' or calcType == 'PM3' or calcType == 'PM6'):
                sem1 = api.resultType(methodology=calcType,spinType='cs:'+molSpin)
                ene1 = api.energiesType(unit='cs:eV')
                ee_ene1 = api.energyType(type_='cs:totalPotential')
                ee_ene1.set_valueOf_(float(molEE))
                hof_ene1 = api.energyType(type_='cs:heatofformation')
                hof_ene1.set_valueOf_(float(data.hofenergies[-1]))
                ene1.add_energy(ee_ene1)
                ene1.add_energy(hof_ene1)
                sem1.set_energies(ene1)
                if hasOrb:
                    sem1.set_waveFunction(wfn1)
                if hasProp:
                    sem1.set_properties(prop1)
                if hasNMR:
                    sem1.set_properties(prop2)
                if hasFreq:
                    sem1.set_vibrationalAnalysis(vib1)
                sdm1.set_semiEmpiricalScf(sem1)
            else:
                print ('The current CSX does not support this method')

            srs1.set_singleDeterminant(sdm1)

        if calcType in md_wfn_method:
            mdm1 = api.srsmdMethodType()
            if (calcType == 'CCSD'):
                ccsd1 = api.resultType(methodology='cs:normal',spinType='cs:'+molSpin, \
                        basisSet='bse:'+basisName)
                ene1 = api.energiesType(unit='cs:eV')
                ee_ene1 = api.energyType(type_='cs:totalPotential')
                ee_ene1.set_valueOf_(float(data.ccenergies[0]))
                ce_ene1 = api.energyType(type_='cs:correlation')
                ce_ene1.set_valueOf_(float(data.ccenergies[-1])-float(molEE))
                ene1.add_energy(ee_ene1)
                ene1.add_energy(ce_ene1)
                ccsd1.set_energies(ene1)
                if hasOrb:
                    ccsd1.set_waveFunction(wfn1)
                if hasProp:
                    ccsd1.set_properties(prop1)
                if hasNMR:
                    ccsd1.set_properties(prop2)
                if hasFreq:
                    ccsd1.set_vibrationalAnalysis(vib1)
                mdm1.set_ccsd(ccsd1)
            elif (calcType == 'CCSD-T'):
                ccsd_t1 = api.resultType(methodology='cs:normal',spinType='cs:'+molSpin, \
                        basisSet='bse:'+basisName)
                ene1 = api.energiesType(unit='cs:eV')
                ee_ene1 = api.energyType(type_='cs:totalPotential')
                ee_ene1.set_valueOf_(float(data.ccenergies[-1]))
                ce_ene1 = api.energyType(type_='cs:correlation')
                ce_ene1.set_valueOf_(float(data.ccenergies[-1])-float(molEE))
                ene1.add_energy(ee_ene1)
                ene1.add_energy(ce_ene1)
                ccsd_t1.set_energies(ene1)
                if hasOrb:
                    ccsd_t1.set_waveFunction(wfn1)
                if hasProp:
                    ccsd_t1.set_properties(prop1)
                if hasNMR:
                    ccsd_t1.set_properties(prop2)
                if hasFreq:
                    ccsd_t1.set_vibrationalAnalysis(vib1)
                mdm1.set_ccsd(ccsd_t1)
            else:
                print ('The current CSX does not support this method')
            srs1.set_multipleDeterminant(mdm1)
        qm1.set_singleReferenceState(srs1)
        mc1.set_quantumMechanics(qm1)
        cs1.set_molecularCalculation(mc1)

        csxfile.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        cs1.export(csxfile,0)
        csxfile.close()


