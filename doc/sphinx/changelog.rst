.. index::
    single: development; changelog

Changelog
=========

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
* Handle Gaussian file with irregular ONION gradients (Tamilmani S)
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
