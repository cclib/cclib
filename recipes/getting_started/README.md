# Getting Started

cclib is a Python library that provides parsers for output files of computational chemistry packages.

cclib has two interfaces a command line interface and a python interface.

## Prerequisites
  - A quantum chemistry output file that you would like to parse. If you don't have one grab one from the test datafiles in the cclib/data folder.
  - An installed version of cclib (See https://cclib.github.io/how_to_install.html)

## Your first extraction (command line interface)

### What can cclib extract from my file?
To see all of the possible attributes that cclib can extract from your output file try running the following command
```
ccget --list [filename]
```

If all goes well, cclib should be able to identify the generating code and parse attributes from the file. The command will list all of the identified attributes. For example
```
$ ccget --list dvb_sp.out 
Attempting to read dvb_sp.out
cclib can parse the following attributes from dvb_sp.out:
  aonames
  aooverlaps
  atombasis
  atomcharges
  atomcoords
  atommasses
  atomnos
  charge
  coreelectrons
  gbasis
  homos
  metadata
  mocoeffs
  moenergies
  moments
  mult
  natom
  nbasis
  nmo
  scfenergies
  scftargets
  scfvalues

```

### Access some data
From the above list, accessing data is as simple as

```
$ ccget moenergies dvb_sp.out
Attempting to read dvb_sp.out
moenergies
[array([-272.56721939, -272.56591324, -272.23768951, -272.2376623 , -272.20762093, -272.20756651, -272.18990632,
       -272.18601509, -271.80666117, -271.80666117,  -21.9284307 ,  -20.41694711,  -19.43559572,  -18.94761395,
        -18.0598153 ,  -15.91596633,  -15.11679516,  -14.36235951,  -13.77802223,  -12.35375112,  -11.85292558,
        -11.08564616,  -10.72359868,  -10.67573385,  -10.09047138,   -9.45100383,   -9.34455289,   -8.73186135,
         -8.36614033,   -7.87407686,   -7.72482241,   -7.07153148,   -5.68083922,   -5.21321157,   -4.06521766,
          1.11969407,    2.55106735,    3.11164909,    5.05723591,    7.50985247,    9.14683497,    9.36879824,
         10.4277565 ,   10.47415191,   11.26540457,   11.28929616,   11.68609258,   12.05581367,   12.44597051,
         13.11651346,   14.40241467,   14.96816657,   15.78581427,   16.38585252,   17.10499501,   17.51156031,
         18.68279275,   19.62210255,   21.30161644,   21.71501181])]

```

## Your first extraction (Python interface)

### What can cclib extract from my file?
Loading an output file. The following code snippet will allow you to read in a data file.
```
>>> import cclib 
>>> data = cclib.io.ccread("dvb_sp.out")
```

If all goes well, cclib should be able to identify the generating code and parse attributes from the file. All of the parsed attributes should be accessible in the ccData object.
```
>>> print(data.moenergies)
[array([-272.56721939, -272.56591324, -272.23768951, -272.2376623 ,
       -272.20762093, -272.20756651, -272.18990632, -272.18601509,
       -271.80666117, -271.80666117,  -21.9284307 ,  -20.41694711,
        -19.43559572,  -18.94761395,  -18.0598153 ,  -15.91596633,
        -15.11679516,  -14.36235951,  -13.77802223,  -12.35375112,
        -11.85292558,  -11.08564616,  -10.72359868,  -10.67573385,
        -10.09047138,   -9.45100383,   -9.34455289,   -8.73186135,
         -8.36614033,   -7.87407686,   -7.72482241,   -7.07153148,
         -5.68083922,   -5.21321157,   -4.06521766,    1.11969407,
          2.55106735,    3.11164909,    5.05723591,    7.50985247,
          9.14683497,    9.36879824,   10.4277565 ,   10.47415191,
         11.26540457,   11.28929616,   11.68609258,   12.05581367,
         12.44597051,   13.11651346,   14.40241467,   14.96816657,
         15.78581427,   16.38585252,   17.10499501,   17.51156031,
         18.68279275,   19.62210255,   21.30161644,   21.71501181])]
```


## How to proceed?
More examples are availible in the `cclib/recipes` directory.
Also further information can be found in the documentation online at https://cclib.github.io/.
If you have any questions start at discussion in the google group at https://groups.google.com/g/cclib.
If you encounter any bugs (or have code to contribute!) visit the github at github.com/cclib/cclib/.
