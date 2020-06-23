.. index::
    module: bridge

Bridges to other packages
=========================

The following bridges in cclib allow using the parsed data to perform further analysis using other packages.

.. index::
    single: bridge; horton

horton
------

`Horton`_ bridge in cclib supports conversion of cclib's ccData object to horton's IOData object and vice versa. This bridge is useful in performing additional population analyses using the data parsed using cclib.
Before invoking the bridge function, ccData object should be prepared first by reading in previous calculations, following the procedures introduced in `how to parse`_ section.
Then, ``ccData`` object can be passed into the bridge function ``makehorton``. Following code block is an example that performs the conversion:

.. _`Horton`: http://theochem.github.io/horton/2.1.1/
.. _`how to parse`: how_to_parse.html

.. code-block:: python

    from cclib.bridge.cclib2horton import makehorton
    from cclib.parser import ccopen

    d = ccopen(sys.argv[1]).parse()
    ht = makehorton(d)

Converted IOData object can be used to run further analyses by referring to `horton documentation`_.

.. _`horton documentation`: http://theochem.github.io/horton/2.1.1/

Becke Population Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

An example code below demonstrates how Becke charges can be calculated based on the script available on `Python interface`_ section in horton documentation.

.. _`Python interface`: https://theochem.github.io/horton/2.1.0/user_postproc_aim.html#python-interface-to-the-partitioning-code

.. code-block:: python

    from cclib.method.density import Density
    from cclib.bridge.cclib2horton import makehorton
    from cclib.parser import ccopen
    
    from horton import BeckeMolGrid, getgobasis, BeckeWPart
    from horton.matrix.dense import DenseTwoIndex

    d = ccopen(sys.argv[1]).parse()
    ht = makehorton(d)
    
    # Calculate density matrix using cclib
    dens = Density(d)
    dens.calculate()
    den = DenseTwoIndex(len(ht.orb_alpha))
    den._array = dens.density[0] + dens.density[1]

    # Create integration grid for calculating Becke charges
    grid = BeckeMolGrid(ht.coordinates, ht.numbers, ht.pseudo_numbers, mode='only')

    # Define Gaussian basis set
    gob = get_gobasis(ht.coordinates, ht.numbers, default = 'STO-3G')
    
    # Partition charges
    wpart = BeckeWPart(ht.coordinates, ht.numbers, ht.pseudo_numbers, grid, moldens, local=True)
    wpart.do_charges()
    print(wpart['charges'])

Hirshfeld Population Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An example code below demonstrates how Hirshfeld charges can be calculated based on the script available on `Python interface`_ section in horton documentation.
To calculate partial charges that require pro-atomic densities, follow the steps in `Building proatomic database`_ section in horton documentation.
Then read in the densities as below to calculate Hirshfeld or Hirshfeld-like charges:

.. _`Python interface`: https://theochem.github.io/horton/2.1.0/user_postproc_aim.html#python-interface-to-the-partitioning-code

.. code-block:: python

    from cclib.method.density import Density
    from cclib.bridge.cclib2horton import makehorton
    from cclib.parser import ccopen
    
    from horton import BeckeMolGrid, getgobasis, HirshfeldWPart
    from horton.matrix.dense import DenseTwoIndex
    from horton.part.proatomdb import ProAtomDB

    d = ccopen(sys.argv[1]).parse()
    ht = makehorton(d)
    
    # Calculate density matrix using cclib
    dens = Density(d)
    dens.calculate()
    den = DenseTwoIndex(len(ht.orb_alpha))
    den._array = dens.density[0] + dens.density[1]

    # Create integration grid
    grid = BeckeMolGrid(ht.coordinates, ht.numbers, ht.pseudo_numbers, mode='only')

    # Define Gaussian basis set
    gob = get_gobasis(ht.coordinates, ht.numbers, default = 'STO-3G')
    
    # Read in pro-atomic density database
    db = ProAtomDB.from_file('atoms.h5')

    # Partition charges
    wpart = HirshfeldWPart(ht.coordinates, ht.numbers, ht.pseudo_numbers, grid, moldens, db)
    wpart.do_charges()
    print(wpart['charges'])

..

