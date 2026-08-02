# cclib_extract_last_geom.py
Extract the last geometry from any quantum chemical output file (not just geometry optimizations) using cclib. Name is the same stub, with the file extension replaced by '.xyz'.

```
usage: extract_last_geom.py [-h] [--fragment] [--trajectory] [--suffix SUFFIX]
                            outputfilename [outputfilename ...]

positional arguments:
  outputfilename

options:
  -h, --help       show this help message and exit
  --fragment       Is this a QChem Fragment calculation?
  --trajectory     Should all geometries from the QChem outputfile be saved?
  --suffix SUFFIX  output geometry format.
```

# qchem_make_opt_input_from_opt.py
Make an input file for a Q-Chem geometry optimization based on the last possible geometry from a Q-Chem geometry optimization; this effectively 'restarts' the geometry with a new filename.

The script assumes the output file being read from is called `*opt(\d*).out`, where 'opt' might be followed by a number. The script will write an input file called `*opt(\d*)+1.in`, with the previous number incremented by one.

```
usage: qchem_make_opt_input_from_opt.py [-h] [--fragment]
                                        outputfilename [outputfilename ...]

positional arguments:
  outputfilename

options:
  -h, --help      show this help message and exit
  --fragment      Is this a QChem Fragment calculation?
```

# qchem_make_freq_input_from_opt.py
Make an input file for a Q-Chem frequency calculation based on the last possible geometry from a Q-Chem geometry optimization.

The script assumes the output file being read from is called `*opt(\d*)*.out`, where 'opt' might be followed by a number. The script will write an input file called `*freq*.in`.

```
usage: qchem_make_freq_input_from_opt.py [-h] [--fragment]
                                         outputfilename [outputfilename ...]

positional arguments:
  outputfilename

options:
  -h, --help      show this help message and exit
  --fragment      Is this a QChem Fragment calculation?
```
