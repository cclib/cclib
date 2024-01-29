.. index::
    single: development; changelog

Changelog
=========

Changes in cclib-1.8.1
----------------------

**Features**

    * New parser: xTB (Andrew S. Rosen, Omri Abarbanel, #789, #1129, #1296)
    * Support parsing number of CPUs and memory used for Gaussian, ORCA, and Turbomole (#1271, #1277)
    * New attributes for the NBO parser (Weronika Zak, #1012, #1251)
    * Remove obsolete Qt4Progress (#1285, #1315)
    * Only keep the Sphinx version of the changelog (#1254, #1352)
    * Stop using Read the Docs (#1267, #1312)
    * Don't require specific versions of PySCF and Biopython optional dependencies (#1344)

**Bugfixes**

    * Fix parsing aooverlaps from GAMESS (#1278)
    * Fix ORCA excited states sometimes not being printed identically across sections (#1276)
    * Make writing ASE objects more robust to missing cclib attributes (Andrew S. Rosen, #1299, #1300)
    * Fix writing freeenergy to ASE objects (Andrew S. Rosen, #1319, #1320)

**Developer facing changes**

    * Completely migrate from ``unittest`` to ``pytest`` (#648, #649, #766, #1309, #1346)
    * Apply and enforce ruff for code formatting using pre-commit (Anselm Hahn, #536, #1122, #1303, #1306, #1307, #1310)
    * Migrate to ``pyproject.toml`` from ``setup.py`` for package installation (#833, #1180, #1311, #1327, #1328)
    * Handle nested dictionaries in setting attributes (Weronika Zak, #1257, #1258)
    * Drop 3.7 and add 3.11 to test suite (#1305)
    * Add Dependabot for keeping GitHub Actions up to date (#1329, #1337)
    * Update copyright years to 2024 with pre-commit enforcement (#1349)
    * Remove obsolete UTF-8 pragmas in Python files (#1350)
    * Fix warning about missing file for a Q-Chem regression test (#1338)
    * Continuous integration fixes for pyquante2 (#997, #1263, #1301, #1304)

Changes in cclib-1.8
--------------------

**Features**

    * [GSoC 2023] New language bindings: Julia (Victor Hugo Cano Gil, #1053, https://github.com/cclib/Cclib.jl)
    * [GSoC 2023] New parser: NBO (Weronika Zak, #1230, #1233, #1244)
    * [GSoC 2023] New parser: GAMESS *.dat (Weronika Zak, #1208, #1214, #1229)
    * New attribute: nmrcouplingtensors for NMR spin-spin couplings, initially from ORCA (#125, #1191)
    * New attribute: rotconsts for rotational constants, initially from Gaussian (Mark Payne, #1054)
    * New method: compute CM5 charges (Sai Murali Karthik Putcha, #676, #852, #1136)
    * Support parsing atomcoords for more DALTON calculation types (#641, #1237)
    * Be less restrictive in parsed data attributes required for Open Babel bridge (#641, #1236)
    * Support parsing implicit solvation parameters from Gaussian, ORCA, QChem, and Turbomole (#1171, #1217, #1232)
    * Support parsing excited state method from Gaussian, ORCA, and Turbomole (#1171, #1219)
    * Support parsing post-HF excited states from ORCA (#1168)
    * Support parsing Hessian from DALTON, GAMESS, NWChem, and Psi4 (#1199, #1200, #1201, #1202)
    * Support parsing Psi4 1.7 (Dustin Wheeler, #1185, #1210)
    * Support writing vibrational frequency attributes to MOLDEN (#1127, #1132)
    * Support writing natural orbitals instead of canonical MOs to MOLDEN (#948)
    * Support parsing Hirshfeld and CM5 charges from Gaussian, ORCA, and QChem (#1137)
    * Support parsing CCSD(T) energies from ORCA and Psi4 (#1194, #1195)
    * Support parsing NWChem 7.0 (#1131, #1133, #1188)
    * Support parsing atommasses and vibdisps from NWChem (#1131, #1198)
    * Support parsing Hartree-Fock and semiempirical excited states from ORCA (#1187, #1189)
    * Support parsing etsyms from ORCA (#1166)
    * Support parsing timing information from Gaussian, ORCA, and Turbomole (#1167)
    * Support parsing methods and basis set from ORCA (#1170)
    * Support parsing excited state using in geometry optimization from Gaussian (#1149, #1151)
    * Support parsing etmagdips from Turbomole (#1139)

**Bugfixes**

    * Fix edge case of XYZ files being interpreted as Turbomole outputs (#1207, #1238)
    * Fix units for calculating nuclear repulsion energy (migatt, #1241, #1242, #1243)
    * Fix infrastructure for reading CJSON (#1222, #1234)
    * Fix edge cases when parsing certain Turbomole outputs (#1174, #1197, #1220, #1222)
    * Fix precision in parsing Hirshfeld charges from ORCA (#1209, #1213)
    * Fix parsing Hessian from formatted checkpoint files (#1204)
    * Fix parsing etoscs for ORCA calculations with spin-orbit coupling (#1172)
    * Fix parsing mocoeffs with non-standard population printing from Gaussian (#1162, #1169)
    * Fix possible infinite loop in DDEC6 method (#1165)
    * Improve checking of mpenergies in Gaussian (#1163, #1164)
    * Fix sign of atomcharges in NWChem (#1156)
    * Update core developers in documentation (#1144)

**Developer facing changes**

    * Update copyright years to 2023 (#1245, #1246)
    * Be more flexible in extend_attribute (Weronika Zak, #1224, #1228)
    * Fix automatically updating both cclib.github.io and cclib.readthedocs.io and testing docs build (#709, #1154, #1158, #1203, #1216, #1223, #1227)
    * Initial Black and isort configuration (#1211)
    * Support development using Dev Containers (#1212)
    * Use raw string in regular expression (#1206)
    * More comprehensive testing of coupled cluster energies (#1196)
    * Continue migration from unittest to pytest (#1181, #1182, #1183, #1184, #1186)
    * Add type annotations to most functions and methods (#1179)
    * More idiomatic checking of None (#991, #1178)
    * Fix installation of pyquante2 in cclib environments (#1176)
    * Test Python 3.9 and 3.10 (#1175)
    * Increase minimum supported Python version to 3.7 (#1157, #1159, #1160, #1161)
    * Modularize atomcharges testing (#1152)
    * Update code coverage Action version (#1095)

Changes in cclib-1.7.2
----------------------

**Features**

    * Support vibfreqs, vibirs, etenergies, etsyms, etoscs and etsecs for NWChem (BenoitDemota)
    * Support temperature, pressure, enthalpy, entropy, zpve and electronic_thermal_energy for NWChem (BenoitDamota)
    * Better metadata support for point group detection
    * Updated code and test file versions to QChem 5.4 and ORCA 5.0

**Bugfixes**

    * Fixed parsing mpenergies for optimization for Turbomole (Oliver Lee)
    * Fixed ccenergies for Gaussian (Oliver Lee)
    * Fixed oscillator strengths for ORCA (Felix Plasser)
    * Fixed units of parsed MO energies for fchk

Changes in cclib-1.7.1
----------------------

**Features**

* New parser: formatted checkpoint files
* New attribute: nmrtensors for nuclear magnetics resonance chemical shielding tensors (Jonathon Vandezande)
* Support atomcharges and atomspins for APT charges in Gaussian (Elliot Farrar)
* Support scannames and scanparms for ORCA logfiles
* Support geometry optimization output and metadata in Turbomole (Oliver Lee)
* Support moments, homos, mosyms, and moenergies in Turbomole (Oliver Lee)
* Support mpenergies and ccenergies in Turbomole (Oliver Lee)
* Support excited state attributes for TD-DFT, CC2 and ADC(2) methods in Turbomole (Oliver Lee)
* Support scfenergies, grad, hessian, atommasses, etenergies and etsyms for fchk output (Javier Cerezo)
* Support zpve for QChem, GAMESS, Psi4, Jaguar, ORCA, DALTON, ADF, GAMESSUK, Molcas and Molpro
* Support walltime and cpu time metadata for QChem output (Amanda Dumi)
* Support walltime and cpu time metadata for Gaussian output (Ellior Farrar)
* Support point group metadata in DALTON
* Plumbed through gbasis and mocoeffs to pyscf bridge (Amanda Dumi)
* Added MO symmetry to Molden writer (Amanda Dumi)

**Bugfixes**

* Improved parsing and testing enthalpy and freeenergy (Felipe Schneider)
* Fixed parsing ONIOM output for Gaussian (Elliot Farrar)
* Fixed parsing of GAMESS logfiles with more than 100 SCF iterations (simonaxelrod)
* Fixed parsing of very long (10K+) ORCA logfiles (Alex Maldonado)
* Fixed parsing of Turbomole outputs that don't compute SCF energies (Oliver Lee)
* Fixed parsing natural charges in Gaussian output
* Fixed parsing vibrational analysis (last, not first) in QChem
* Fixed indices for open shell systems in QChem (Hubert Weißmann)
* Cleaned up Turbomole unit test logfiles (froessler)
* Updated documentation for grads (Cyrille Lavigne)

Changes in cclib-1.7
--------------------

**Features**

* Dropping support for Python 2
* SciPy is now a hard dependency for cclib

**Bugfixes**

* Fixed parsing of Gaussian files missing scftargets (Hubert Weißmann)
* Fixed parsing TDA excited states from QChem (srtlg)
* Fixed parsing two character elements from Turbomole

Changes in cclib-1.6.4
----------------------

**Features**

* [GSOC2020] New methods: Bader's QTAIM, Bickelhaupt, Stockholder, Hirshfeld, and DDEC6 partial charges (Minsik Cho)
* [GSOC2020] New bridge to Horton (Minsik Cho)
* [GSOC2020] Support reading cube files in volume method (Minsik Cho)
* New bridge to Atomic Simulation Environment (Felipe S. S. Schneider)
* New bridge to PySCF (Amanda Dumi)
* New attribute dispersionenergies for molecular dispersion energy corrections
* New attribute vibfconsts for vibrational force constants (Chikashi Shinagawa)
* New attribute vibrmasses for vibrational reduced masses (Chikashi Shinagawa)
* Support t1_diagnostic in metadata for most parsers

**Bugfixes**

* Fixed parsing of ORCA optimization with constraints (Jonathon Vandezande)
* Fixed parsing of too many excited states in Gaussian09 optimization (Oliver Lee)
* Fixed parsing Gaussian logfiles with NQMF / reduced number of atoms (Michael D'Addario)
* Fixed bug in QChem parser related to two letter chemical symbols (Amanda Dumi)
* Fixed Gaussian grads to align with standard orientation like other attributes (Chikashi Shinagawa)
* Fixed handling of open shell systems in modelwriter and wfxwriter (Dave Z.)

Changes in cclib-1.6.3
----------------------

**Features**

* New bridge to Psi4 (Felipe S. S. Schneider)
* New attribute zpve for zero-point vibrational energy correction (kuriba)
* New attributes for electric transition dipoles of electronic transitions (mwykes)
* Support ccenergies in ORCA
* Support mpenergies in ORCA (Alex Maldonado)
* Support grads in MOLCAS (Daniele Padula)
* Support Mulliken atomspins in Gaussian (Peter St. John)
* Support temperature, pressure, enthalpy, entropy and freenergy attributes in GAMESS (Mark Perri)
* Support fuzzy matching of attribute in ccget script
* Updated test file versions to Psi4 1.3.1, and ORCA 4.2

**Bugfixes**

* Fixed parsing of vibrational attribute for single atoms in ORCA (Felipe S. S. Schneider)
* Fixed parsing very long ORCA logfiles (Alex Maldonado)
* Fixed method code for principal moments of inertia, and mulliken charges in Gaussian (James E T Smith)
* Fixed scannames, scanparm and scanenergies in Gaussian (Dustin Wheeler)
* Fixed freeenergy in ORCA 4.2 (shijunang)
* Fixed name collisions in tests and use of periodic table in utilities (Waylon Peng)

Changes in cclib-1.6.2
----------------------

**Features**

* Molden writer now supports ghost atoms (Shiv Upadhyay)
* Handle comments in XYZ files when reading and writing
* Updated regression testing framework (Amanda Dumi, Shiv Upadhyay)
* Updated test file versions to GAMESS-US 2018 (Shiv Upadhyay)

**Bugfixes**

* Fixed parsing ORCA output with user comments in coordinates (Felix Plasser)
* Fixed parsing ORCA output with embedding potentials
* Fixed parsing ORCA output with ROCIS in version 4.1
* Fixed parsing etenergies and similar attribute in ORCA for excited states
* Fixed parsing of vibfreqs for ORCA for linear molecules
* Parsing geometry optimization in ORCA is mode robust wrt line endings

Changes in cclib-1.6.1
----------------------

**Features**

* New attribute nsocoeffs for natural spin orbital coefficients (Shiv Upadhyay)
* New attribute nsooccnos for natural spin orbital occupation numbers (Shiv Upadhyay)
* New methods: alpha and beta electron counts (Jaime Rodríguez-Guerra)
* Support coreelectrons attribute in Molcas (Kunal Sharma)
* Support etoscs for response calculations in Dalton (Peter Reinholdt)
* Support etenergies for TDDFT in GAMESS
* Support etrotats attribute in ORCA
* Support functional name in metadata for Psi4 (Alessandro Genova)
* Updated testing framework (Jaime Rodríguez-Guerra, Maxim Stolyarchuk and others)
* Updated test file version to QChem 5.1

**Bugfixes**

* Fixed parsing GAMESS output for EOM-CC output
* Fixed parsing Gaussian output for G3 jobs
* Fixed parsing ORCA output for certain invalid inputs (Felipe S. S. Schneider)
* Fixed parsing of mocoeffs in ORCA when they are glued together (Felipe S. S. Schneider)
* Fixed parsing of mocoeffs and vibfreqs in Psi4 (Alessandro Genova)
* Fixed parsing of mocoeffs in Molcas for some files (Shiv Upadhyay)
* Fixed parsing of etsecs in Dalton
* Fixed bond atom indices in CJSON output (Alessandro Genova)

Changes in cclib-1.6
--------------------

**Features**

* New parser: cclib can now parse Molcas files (Kunal Sharma)
* New parser: cclib can now parse Turbomole files (Christopher Rowley, Kunal Sharma)
* New script: ccframe writes data table files from logfiles (Felipe Schneider)
* New method: stoichiometry builds the chemical formula of a system (Jaime Rodríguez-Guerra)
* Support package version in metadata for most parsers
* Support time attribute and BOMD output in Gaussian, NWChem, ORCA and QChem
* Support grads and metadata attributes in ORCA (Jonathon Vandezande)
* Experimental support for CASSCF output in ORCA (Jonathon Vandezande)
* Added entry in metadata for successful completion of jobs
* Updated test file versions to ORCA 4.0
* Update minimum Python3 version to 3.4

**Bugfixes**

* Fixed parsing ORCA output with linear molecules (Jonathon Vandezande)
* Fixed parsing NWChem output with incomplete SCF

Changes in cclib-1.5.3
----------------------

**Features**

* New attribute transprop for electronic transitions (Jonathon Vandezande)
* Support grads attribute in Psi4 (Adam Abbott)
* Support grads attribute in Molpro (Oskar Weser)
* Support optstatus for IRCs and in Psi4 (Emmanuel LaTruelle)
* Updated test file versions to Gaussian16 (Andrew S. Rosen)
* Add ability to write XYZ coordinates for arbitrary indices

**Bugfixes**

* Fixed ccwrite script and added unit tests (Georgy Frolov)
* Fixed closed shell determination for Gaussian (Jaime Rodríguez-Guerra)
* Fixed parsing of natom for >9999 atoms in Gaussian (Jaime Rodríguez-Guerra)
* Fixed parsing of ADF jobs with no title
* Fixed parsing of charge and core electrons when using ECPs in QChem
* Fixed parsing of scfvalues for malformed output in Gaussian

Changes in cclib-1.5.2
----------------------

**Features**

* Support for writing Molden and WFX files (Sagar Gaur)
* Support for thermochemistry attributes in ORCA (Jonathon Vandezande)
* Support for chelpg atomic charges in ORCA (Richard Gowers)
* Updated test file versions to GAMESS-US 2017 (Sagar Gaur)
* Added option to print full arrays with ccget (Sagar Gaur)

**Bugfixes**

* Fixed polarizability parsing bug in DALTON (Maxim Stolyarchuk)
* Fixed IRC parsing in Gaussian for large trajectories (Dénes Berta, LaTruelle)
* Fixed coordinate parsing for heavy elements in ORCA (Jonathon Vandezande)
* Fixed parsing of large mocoeffs in fixed width format for QChem (srtlg)
* Fixed parsing of large polarizabilities in fixed width format for DALTON (Maxim Stolyarchuk)
* Fixed parsing molecular orbitals when there are more than basis set functions in QChem

Changes in cclib-1.5.1
----------------------

**Features**

* New attribute polarizabilities for static or dynamic dipole polarizability
* New attribute pressure for thermochemistry (renpj)
* Add property to detect closed shells in parsed data
* Handle RPA excited state calculation in ORCA, in addition to TDA
* Support for Python 3.6

**Bugfixes**

* Restore alias cclib.parser.ccopen for backwards compatibility
* Fixed parsing thermochemistry for single atoms in QChem
* Fixed handling of URLs (Alexey Alnatanov)
* Fixed Atom object creation in Biopython bridge (Nitish Garg)
* Fixed ccopen when working with multiple files

Changes in cclib-1.5
--------------------

**Features**

* Support for both reading and writing CJSON (Sanjeed Schamnad)
* New parser: cclib can now parse MOPAC files (Geoff Hutchison)
* New attribute time tracks coordinated for dynamics jobs (Ramon Crehuet)
* New attribute metadata holds miscellaneous information not in other attributes (bwang2453)
* Extract moments attribute for Gaussian (Geoff Hutchison)
* Extract atombasis for ADF in simple cases (Felix Plasser)
* License change to BSD 3-Clause License

**Bugfixes**

* Correct parsing of several attributes for ROHF calculations
* Fixed precision of scfvalues in ORCA
* Fixed MO parsing from older versions of Firefly (mkrompiec)

Changes in cclib-1.4.1
----------------------

**Features**

* Preliminary support for writing CJSON (Sanjeed Schamnad)
* Tentative support for BOMD trajectories in Gaussian (Ramon Crehuet)
* Support for atombasis in ADF (Felix Plasser)
* Support for nocoeffs and nooccnos in Molpro

**Bugfixes**

* Fix for non-standard basis sets in DALTON
* Fix for non-standard MO coefficient printing in GAMESS

Changes in cclib-1.4
--------------------

**Features**

* New parser: cclib can now parse DALTON files
* New parser: cclib can now parse ORCA files
* New attribute optstatus for status during geometry optimizations and scans
* Extract atommasses for GAMESS-US (Sagar Gaur)
* Extract atombasis, gbasis and mocoeffs for QChem
* Extract gbasis for ORCA (Felix Plasser)
* Handle multi-step jobs by parsing only the supersystem
* Improve parsing vibrational symmetries and displacements for Gaussian (mwykes)
* Improve support for compressed files (mwykes)
* Improve and update unit test and regression suites
* Support for Python 3.5

**Bugfixes**

* Fix StopIteration crashes for most parsers
* Fix parsing basis section for Molpro job generated by Avogadro
* Fix parsing multi-job Gaussian output with different orbitals (Geoff Hutchinson)
* Fix parsing ORCA geometry optimization with improper internal coordinates (glideht)
* Fix units in atom coordinates parsed from GAMESS-UK files (mwykes)
* Fix test for vibrational frequencies in Turbomole (mwykes)
* Fix parsing vibration symmetries for Molpro (mwykes)
* Fix parsing eigenvectors in GAMESS-US (Alexis Otero-Calvis)
* Fix duplicate parsing of symmetry labels for Gaussian (Martin Peeks)

Changes in cclib-1.3.2
----------------------

**Features**

* New attribute nooccnos for natural orbital occupation numbers
* Read data from XYZ files using Open Babel bridge
* Start basic tests for bridge functionality

**Bugfixes**

* Better handling of ONIOM logfiles in Gaussian (Clyde Fare)
* Fix IR intensity bug in Gaussian parser (Clyde Fare)
* Fix QChem parser for OpenMP output
* Fix parsing TDDFT/RPA transitions (Felix Plasser)
* Fix encoding issues for UTF-8 symbols in parsers and bridges

Changes in cclib-1.3.1
----------------------

**Features**

* New attribute nooccnos for natural orbital occupation numbers
* Read data from XYZ files using Open Babel bridge
* Start basic tests for bridge functionality

**Bugfixes**

* Better handling of ONIOM logfiles in Gaussian (Clyde Fare)
* Fix IR intensity bug in Gaussian parser (Clyde Fare)
* Fix QChem parser for OpenMP output
* Fix parsing TDDFT/RPA transitions (Felix Plasser)
* Fix encoding issues for UTF-8 symbols in parsers and bridges

Changes in cclib-1.3
--------------------

**Features**

* New parser: cclib can now parse NWChem files
* New parser: cclib can now parse Psi (versions 3 and 4) files
* New parser: cclib can now parse QChem files (by Eric Berquist)
* New method: Nuclear (currently calculates the repulsion energy)
* Handle Gaussian basis set output with GFPRINT keyword
* Attribute optdone reverted to single Boolean value by default
* Add --verbose and --future options to ccget and parsers
* Replaced PC-GAMESS test files with newer Firefly versions
* Updated test file versions to GAMESS-UK 8.0

**Bugfixes**

* Handle GAMESS-US file with LZ value analysis (Martin Rahm)
* Handle Gaussian jobs with stars in output (Russell Johnson, NIST)
* Handle ORCA singlet-only TD calculations (May A.)
* Fix parsing of Gaussian jobs with fragments and ONIOM output
* Use UTF-8 encodings for files that need them (Matt Ernst)

Changes in cclib-1.2
--------------------

**Features**

* Move project to GitHub
* Transition to Python 3 (Python 2.7 will still work)
* Add a multifile mode to ccget script
* Extract vibrational displacements for ORCA
* Extract natural atom charges for Gaussian (Fedor Zhuravlev)
* Updated test file versions to ADF2013.01, GAMESS-US 2012, Gaussian09, Molpro 2012 and ORCA 3.0.1

**Bugfixes**

* Ignore Unicode errors in logfiles
* Handle Gaussian jobs with terse output (basis set count not reported)
* Handle Gaussian jobs using IndoGuess (Scott McKechnie)
* Handle Gaussian file with irregular ONIOM gradients (Tamilmani S)
* Handle ORCA file with SCF convergence issue (Melchor Sanchez)
* Handle Gaussian file with problematic IRC output (Clyde Fare)
* Handle ORCA file with AM1 output (Julien Idé)
* Handle GAMESS-US output with irregular frequency format (Andrew Warden)

Changes in cclib-1.1
--------------------

**Features**

* Add progress info for all parsers
* Support ONIOM calculations in Gaussian (Karen Hemelsoet)
* New attribute atomcharges extracts Mulliken and Löwdin atomic charges if present
* New attribute atomspins extracts Mulliken and Löwdin atomic spin densities if present
* New thermodynamic attributes: freeenergy, temperature, enthalpy (Edward Holland)
* Extract PES information: scanenergies, scancoords, scanparm, scannames (Edward Holland)

**Bugfixes**

* Handle coupled cluster energies in Gaussian 09 (Björn Dahlgren)
* Vibrational displacement vectors missing for Gaussian 09 (Björn Dahlgren)
* Fix problem parsing vibrational frequencies in some GAMESS-US files
* Fix missing final scfenergy in ADF geometry optimisations
* Fix missing final scfenergy for ORCA where a specific number of SCF cycles has been specified
* ORCA scfenergies not parsed if COSMO solvent effects included
* Allow spin unrestricted calculations to use the fragment MO overlaps correctly for the MPA and CDA calculations
* Handle Gaussian MO energies that are printed as a row of asterisks (Jerome Kieffer)
* Add more explicit license notices, and allow LGPL versions after 2.1
* Support Firefly calculations where nmo != nbasis (Pavel Solntsev)
* Fix problem parsing vibrational frequency information in recent GAMESS (US) files (Chengju Wang)
* Apply patch from Chengju Wang to handle GAMESS calculations with more than 99 atoms
* Handle Gaussian files with more than 99 atoms having pseudopotentials (Björn Baumeier)

Changes in cclib-1.0.1
----------------------

**Features**

* New attribute atommasses - atomic masses in Dalton
* Added support for Gaussian geometry optimisations that change the number of linearly independent basis functions over the course of the calculation

**Bugfixes**

* Handle triplet PM3 calculations in Gaussian03 (Greg Magoon)
* Some Gaussian09 calculations were missing atomnos (Marius Retegan)
* Handle multiple pseudopotentials in Gaussian03 (Tiago Silva)
* Handle Gaussian calculations with >999 basis functions
* ADF versions > 2007 no longer print overlap info by default
* Handle parsing Firefly calculations that fail
* Fix parsing of ORCA calculation (Marius Retegan)

Changes in cclib-1.0
--------------------

**Features**

* Handle PBC calculations from Gaussian
* Updates to handle Gaussian09
* Support TDDFT calculations from ADF
* A number of improvements for GAMESS support
* ccopen now supports any file-like object with a read() method, so it can parse across HTTP

**Bugfixes**

* Many many additional files parsed thanks to bugs reported by users

Changes in cclib-0.9
--------------------

**Features**

* New parser: cclib can now parse ORCA files
* Added option to use setuptools instead of distutils.core for installing
* Improved handling of CI and TD-DFT data: TD-DFT data extracted from GAMESS and etsecs standardised across all parsers
* Test suite changed to include output from only the newest program versions

**Bugfixes**

* A small number of parsing errors were fixed

Changes in cclib-0.8
--------------------

**Feaures**

* New parser: cclib can now parse Molpro files
* Separation of parser and data objects: Parsed data is now returned is a ccData object that can be pickled, and converted to and from JSON
* Parsers: multiple files can be parsed with one parse command
* NumPy support: Dropped Numeric support in favour of NumPy
* API addition: 'charge' for molecular charge
* API addition: 'mult' for spin multiplicity
* API addition: 'atombasis' for indices of atom orbitals on each atom
* API addition: 'nocoeffs' for Natural Orbital (NO) coefficients
* GAMESS-US parser: added 'etoscs' (CIS calculations)
* Jaguar parser: added 'mpenergies' (LMP2 calculations)
* Jaguar parser: added 'etenergies' and 'etoscs' (CIS calculations)
* New method: Lowdin Population Analysis (LPA)
* Tests: unittests can be run from the Python interpreter, and for a single parser; the number of "passed" tests is also counted and shown

**Bugfixes**

* Several parsing errors were fixed
* Fixed some methods to work with different numbers of alpha and beta MO coefficients in mocoeffs (MPA, CSPA, OPA)

Changes in cclib-0.7
--------------------

**Feaures**

* New parser: cclib can now parse Jaguar files
* ccopen: Can handle log files which have been compressed into .zip, .bz2 or .gz files.
* API addition: 'gbasis' holds the Gaussian basis set
* API addition: 'coreelectrons' contains the number of core electrons in each atom's pseudopotential
* API addition: 'mpenergies' holds the Moller-Plesset corrected molecular electronic energies
* API addition: 'vibdisps' holds the Cartesian displacement vectors
* API change: 'mocoeffs' is now a list of rank 2 arrays, rather than a rank 3 array
* API change: 'moenergies' is now a list of rank 1 arrays, rather than rank 2 array
* GAMESS-UK parser: added 'vibramans'
* New method: Charge Decomposition Analysis (CDA) for studying electron donation, back donation, and repulsion between fragments in a molecule
* New method: Fragment Analysis for studying bonding interactions between two or more fragments in a molecule
* New method: Ability to calculate the electron density or wavefunction

**Bugfixes**

* GAMESS parser:
    - Failed to parse frequency calculation with imaginary frequencies
    - Rotations and translations now not included in frequencies
    - Failed to parse a DFT calculation
* GAMESS-UK parser:
    - 'atomnos' not being extracted
    - Rotations and translations now not included in frequencies
* bridge to Open Babel: No longer dependent on pyopenbabel

Changes in cclib-0.6.1
----------------------

**Bugfixes**

* cclib: The "import cclib.parsers" statement failed due to references to Molpro and Jaguar parsers which are not present
* Gaussian parser: Failed to parse single point calculations where the input coords are a z-matrix, and symmetry is turned off.

Changes in cclib-0.6.0
----------------------

**Feaures**

* ADF parser: If some MO eigenvalues are not present, the parser does not fail, but uses values of 99999 instead and A symmetry

**Bugfixes**

* ADF parser: The following bugs have been fixed P/D orbitals for single atoms not handled correctly Problem parsing homos in unrestricted calculations Problem skipping the Create sections in certain calculations
* Gaussian parser: The following bugs have been fixed Parser failed if standard orientation not found
* ccget: aooverlaps not included when using --list option

Changes in cclib-0.6b
---------------------

**Feaures**

* New parser: GAMESS-UK parser
* API addition: the .clean() method; the .clean() method of a parser clears all of the parsed attributes. This is useful if you need to reparse during the course of a calculation.
* Function rename: guesstype() has been renamed to ccopen()
* Speed up: Calculation of Overlap Density of States has been sped up by two orders of magnitude

**Bugfixes**

* ccopen: Minor problems fixed with identification of log files
* ccget: Passing multiple filenames now works on Windows too
* ADF parser: The following bugs have been fixed
    - Problem with parsing SFOs in certain log files
    - Handling of molecules with orbitals of E symmetry
    - Couldn't find the HOMO in log files from new versions of ADF
    - Parser used to miss attributes if SCF not converged
    - For a symmetrical molecule, mocoeffs were in the wrong order and the homo was not identified correctly if degenerate
* Gaussian parser: The following bugs have been fixed
    - SCF values was not extracting the dEnergy value
    - Was extracting Depolar P instead of Raman activity

Changes in cclib-0.5
--------------------

**Features**

* (src/scripts/ccget): Added handling of multiple filenames. It's now possible to use ccget as follows: ``ccget *.log``. This is a good way of checking out whether cclib is able to parse all of the files in a given directory. Also possible is: ``ccget homos *.log``.
* Change of license: Changed license from GPL to LGPL

**Bugfixes**

* src/cclib/parser/gamessparser.py: gamessparser was dying on GAMESS VERSION = 12 DEC 2003 gopts, as it was unable to parse the scftargets.
* src/cclib/parser/gamessparser.py: Remove assertion to catch instances where scftargets is unset. This occurs in the case of failed calculations (e.g. wrong multiplicity).
* src/cclib/parser/adfparser.py: Fixed one of the errors with the Mo5Obdt2-c2v-opt.adfout example, which had to do with the SFOs being made of more than two combinations of atoms (4, because of rotation in c2v point group). At least one error is still present with atomcoords. It looks like non-coordinate integers are being parsed as well, which makes some of the atomcoords list have more than the 3 values for x,y,z.
* src/cclib/parser/adfparser.py: Hopefully fixed the last error in Mo5Obdt2-c2v-opt. Problem was that it was adding line.split()[5:], but sometimes there was more than 3 fields left, so it was changed to [5:8]. Need to check actual parsed values to make sure it is parsed correctly.
* data/Gaussian, logfiledist, src/cclib/parser/gaussianparser.py, test/regression.py: Bug fix: Mo4OSibdt2-opt.log has no atomcoords despite being a geo-opt. This was due to the fact that the parser was extracting "Input orientation" and not "Standard orientation". It's now changed to "Standard orientation" which works for all of the files in the repository.
