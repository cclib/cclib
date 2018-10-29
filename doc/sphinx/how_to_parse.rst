How to parse and write
======================

This page outlines the various ways cclib can be used to parse and write logfiles, and provides several examples to get you started.

From Python
+++++++++++

Importing cclib and parsing a file is a few lines of Python code, making it simple to access data from the output file of any supported computational chemistry program. For example:

.. code-block:: python

  import cclib

  filename = "logfile.out"
  parser = cclib.io.ccopen(filename)
  data = parser.parse()
  print("There are %i atoms and %i MOs" % (data.natom, data.nmo))

A newer command, ``ccread``, combines both the format detection and parsing steps:

.. code-block:: python

  import cclib

  filename = "logfile.out"
  data = cclib.io.ccread(filename)
  print("There are %i atoms and %i MOs" % (data.natom, data.nmo))

From command line
+++++++++++++++++

The cclib package provides four scripts to parse and write data: ``ccget``, ``ccwrite``, ``cda``, and ``ccframe``.

1. **ccget** is used to parse attribute data from output files.
2. **ccwrite** has the ability to list out all valid attribute data that can be parsed from an output format. It has the added feature of writing the output file into four different formats i.e. ``json``, ``cjson``, ``cml``, ``xyz``.
3. **cda** is used for the chemical decomposition analysis of output files.
4. **ccframe** is used to write data tables from output files.

This page describes how to use the ccget, ccwrite and ccframe scripts to obtain data from output files.

ccget
-----

The data types that can be parsed from the output file depends on the type of computation being conducted. The name of the output file used to show example usage is ``Benzeneselenol.out``.

Data type can be parsed from the output file by following this format::

    ccget <attribute> [<attribute>] <CompChemLogFile> [<CompChemLogFile>]

where ``attribute`` can be any one of the attribute names available `here`_.

.. _`here`: data_dev.html

1. Atomic Charges

    The atomic charges are obtained by using the ``atomcharges`` attribute::

        $ ccget atomcharges Benzeneselenol.out
        Attempting to read Benzeneselenol.out
        atomcharges:
        {'mulliken': array([-0.49915 ,  0.056965,  0.172161,  0.349794, -0.153072,  0.094583,
            0.016487,  0.050249,  0.002149,  0.01161 ,  0.053777, -0.173671,
            0.018118])}

2. Electronic Energies

    The molecular electronic energies after SCF (DFT) optimization of the input molecule are printed by using the ``scfenergies`` attribute::

        $ ccget scfenergies Benzeneselenol.out
        Attempting to read Benzeneselenol.out
        scfenergies:
        [-71671.43702915 -71671.4524142  -71671.4534768  -71671.45447492
        -71671.4556548  -71671.45605671 -71671.43194906 -71671.45761021
        -71671.45850275 -71671.39630296 -71671.45915119 -71671.45935854
        -71671.4594614  -71671.45947338 -71671.45948807 -71671.4594946
        -71671.4594946 ]


3. Geometry Targets

    The targets for convergence of geometry optimization can be obtained by using the ``geotargets`` attribute::

        $ ccget geotargets Benzeneselenol.out
        Attempting to read Benzeneselenol.out
        geotargets:
        [ 0.00045  0.0003   0.0018   0.0012 ]

Chaining of attributes
^^^^^^^^^^^^^^^^^^^^^^

ccget provides the user with the option to chain attributes to obtain more than one type of data with a command call. The attributes can be chained in any particular order. A few chained examples are provided below.

1. Molecular Orbitals and Multiplicity

    The number of molecular orbitals and the number of basis functions used to optimize the molecule can be obtained by running the following command::

        $ ccget nmo nbasis Benzeneselenol.out
        Attempting to read Benzeneselenol.out
        nmo:
        405
        nbasis:
        407

2. Enthalpy and Vibrational Frequency

    The enthalpy and the vibrational frequencies of the optimized molecule is conducted is obtained below::

        $ ccget enthalpy vibfreqs Benzeneselenol.out
        Attempting to read Benzeneselenol.out
        enthalpy:
        -2633.77264
        vibfreqs:
        [  129.5512   170.6681   231.4278   304.8614   407.8299   472.5026
           629.9087   679.9032   693.2509   746.7694   812.5113   850.2578
           915.8742   987.1252   988.1785  1002.8922  1038.1073  1091.4005
          1102.3417  1183.3857  1209.2727  1311.3497  1355.6441  1471.4447
          1510.1919  1611.9088  1619.0156  2391.2487  3165.1596  3171.3909
          3182.0753  3188.5786  3198.0359]

ccwrite
-------

The same Benzeneselenol.out file used in the previous examples will be used as the input file for ccwrite. When the ccwrite script is used with a valid input, it prints out the valid attributes that can be parsed from the file.

Command line format::

    ccwrite <OutputFileFormat>  <CompChemLogFile> [<CompChemLogFile>]

The valid output file formats are ``cjson``, ``cml``, and ``xyz``.

1. `Chemical markup language`_ (CML)::

    $ ccwrite cml Benzeneselenol.out
    Attempting to parse Benzeneselenol.out
    cclib can parse the following attributes from Benzeneselenol.out:
      atomcharges
      atomcoords
      atomnos
      charge
      coreelectrons
      enthalpy
      geotargets
      geovalues
      grads
      homos
      moenergies
      mosyms
      mult
      natom
      nbasis
      nmo
      optdone
      optstatus
      scfenergies
      scftargets
      temperature
      vibdisps
      vibfreqs
      vibirs
      vibsyms

.. _`chemical markup language`: http://www.xml-cml.org/

A ``Benzeneselenol.cml`` output file is generated in the same directory as the ``Benzeneselenol.out`` file:

.. code-block:: xml

    <?xml version='1.0' encoding='utf-8'?>
    <molecule id="Benzeneselenol.out" xmlns="http://www.xml-cml.org/schema">
      <atomArray>
        <atom elementType="C" id="a1" x3="-2.8947620000" y3="-0.0136420000" z3="-0.0015280000" />
        <atom elementType="C" id="a2" x3="-2.2062510000" y3="1.1938510000" z3="-0.0025210000" />
        <atom elementType="C" id="a3" x3="-0.8164260000" y3="1.2153020000" z3="-0.0022010000" />
        <atom elementType="C" id="a4" x3="-0.1033520000" y3="0.0183920000" z3="0.0031060000" />
        <atom elementType="C" id="a5" x3="-0.7906630000" y3="-1.1943840000" z3="0.0058500000" />
        <atom elementType="C" id="a6" x3="-2.1799570000" y3="-1.2059710000" z3="0.0017890000" />
        <atom elementType="H" id="a7" x3="-3.9758430000" y3="-0.0253010000" z3="-0.0029040000" />
        <atom elementType="H" id="a8" x3="-2.7502340000" y3="2.1291370000" z3="-0.0052760000" />
        <atom elementType="H" id="a9" x3="-0.2961840000" y3="2.1630180000" z3="-0.0073260000" />
        <atom elementType="H" id="a10" x3="-0.2474670000" y3="-2.1302310000" z3="0.0132260000" />
        <atom elementType="H" id="a11" x3="-2.7028960000" y3="-2.1530750000" z3="0.0036640000" />
        <atom elementType="Se" id="a12" x3="1.8210800000" y3="-0.0433780000" z3="-0.0038760000" />
        <atom elementType="H" id="a13" x3="2.0043580000" y3="1.4100070000" z3="0.1034490000" />
      </atomArray>
      <bondArray>
        <bond atomRefs2="a9 a3" order="1" />
        <bond atomRefs2="a8 a2" order="1" />
        <bond atomRefs2="a12 a4" order="1" />
        <bond atomRefs2="a12 a13" order="1" />
        <bond atomRefs2="a7 a1" order="1" />
        <bond atomRefs2="a2 a3" order="2" />
        <bond atomRefs2="a2 a1" order="1" />
        <bond atomRefs2="a3 a4" order="1" />
        <bond atomRefs2="a1 a6" order="2" />
        <bond atomRefs2="a6 a11" order="1" />
        <bond atomRefs2="a6 a5" order="1" />
        <bond atomRefs2="a4 a5" order="2" />
        <bond atomRefs2="a5 a10" order="1" />
      </bondArray>
    </molecule>

2. XYZ_

Using ``xyz`` as the ``<OutputFileFormat>`` with ``Benzeneselenol.out``, we obtain the following ``Benzeneselenol.xyz`` file::

    13
    Benzeneselenol.out: Geometry 17
    C     -2.8947620000   -0.0136420000   -0.0015280000
    C     -2.2062510000    1.1938510000   -0.0025210000
    C     -0.8164260000    1.2153020000   -0.0022010000
    C     -0.1033520000    0.0183920000    0.0031060000
    C     -0.7906630000   -1.1943840000    0.0058500000
    C     -2.1799570000   -1.2059710000    0.0017890000
    H     -3.9758430000   -0.0253010000   -0.0029040000
    H     -2.7502340000    2.1291370000   -0.0052760000
    H     -0.2961840000    2.1630180000   -0.0073260000
    H     -0.2474670000   -2.1302310000    0.0132260000
    H     -2.7028960000   -2.1530750000    0.0036640000
    Se     1.8210800000   -0.0433780000   -0.0038760000
    H      2.0043580000    1.4100070000    0.1034490000

.. _XYZ: https://en.wikipedia.org/wiki/XYZ_file_format

ccframe
-------

This script creates complete tables of data from output files in some of the formats supported by pandas_.
Since the pandas library is not a dependency of cclib, `it must be installed <https://pandas.pydata.org/pandas-docs/stable/install.html>`_ separately.

.. _pandas: https://pandas.pydata.org/

A complete data table can be parsed from many output files by following this format::

    ccframe -O <OutputDest> <CompChemLogFile> [<CompChemLogFile>...]

The argument for ``-O`` indicates the data file to be written and its extension specifies the filetype (e.g. csv, h5/hdf/hdf5, json, pickle/pkl, xlsx).
Since higher-dimensional attributes (e.g. ``atomcoords``) are handled as plain text in some file formats (such as Excel XLSX or CSV), we recommend storing JSON or HDF5 files.
