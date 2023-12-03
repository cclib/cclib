.. index::
    module: methods

Calculation methods
===================

The following methods in cclib allow further analysis of calculation output.

.. _`methods in the development version`: methods_dev.html

.. index::
    single: methods; C squared population analysis (CSPA)

C squared population analysis (CSPA)
------------------------------------

**CSPA** can be used to determine and interpret the electron density of a molecule. The contribution of the a-th atomic orbital to the i-th molecular orbital can be written in terms of the molecular orbital coefficients:

.. math:: \Phi_{ai} = \frac{c^2_{ai}}{\sum_k c^2_{ki}}

The CSPA class available from cclib.method performs C-squared population analysis and can be used as follows:

.. code-block:: python

    from cclib.io import ccread
    from cclib.method import CSPA

    data = ccread("mycalc.out")

    m = CSPA(data)
    m.calculate()

After the ``calculate()`` method is called, the following attributes are available:

* ``aoresults`` is a NumPy array[3] with spin, molecular orbital, and atomic/fragment orbitals as the axes (``aoresults[0][45][0]`` gives the contribution of the 1st atomic/fragment orbital to the 46th alpha/restricted molecular orbital)
* ``fragresults`` is a NumPy array[3] with spin, molecular orbital, and atoms as the axes (``atomresults[1][23][4]`` gives the contribution of the 5th atomic/fragment orbital to the 24th beta molecular orbital)
* ``fragcharges`` is a NumPy array[1] with the number of (partial) electrons in each atom (``atomcharges[2]`` gives the number of electrons on the 3rd atom)

Custom fragments
~~~~~~~~~~~~~~~~

Calling the calculate method without an argument treats each atom as a fragment in the population analysis. An optional argument can be passed - a list of lists - containing the atomic orbital numbers to be included in each fragment. Calling with this additional argument is useful if one is more interested in the contributions of certain orbitals, such as metal d, to the molecular orbitals. For example:

.. code-block:: python

    from cclib.io import ccread
    from cclib.method import CSPA

    data = ccread("mycalc.out")

    m = CSPA(data)
    m.calculate([[0, 1, 2, 3, 4], [5, 6], [7, 8, 9]]) # fragment one is made from basis functions 0 - 4
                                                      # fragment two is made from basis functions 5 & 6
                                                      # fragment three is made from basis functions 7 - 9

Custom progress
~~~~~~~~~~~~~~~

The CSPA class also can take a progress class as an argument so that the progress of the calculation can be monitored:

.. code-block:: python

    from cclib.method import CSPA
    from cclib.parser import Gaussian
    from cclib.progress import TextProgress

    import logging

    progress = TextProgress()
    p = Gaussian("mycalc.out", logging.ERROR)
    d = p.parse(progress)

    m = CSPA(d, progress, logging.ERROR)
    m.calculate()

.. index::
    single: methods; Mulliken population analysis (MPA)

Mulliken population analysis (MPA)
----------------------------------

MPA can be used to determine and interpret the electron density of a molecule. The contribution of the a-th atomic orbital to the i-th molecular orbital in this method is written in terms of the molecular orbital coefficients, c, and the overlap matrix, S:

.. math:: \Phi_{ai} = \sum_b c_{ai} c_{bi} S_{ab}

The MPA class available from cclib.method performs Mulliken population analysis and can be used as follows:

.. code-block:: python

    import sys

    from cclib.method import MPA
    from cclib.parser import ccopen

    d = ccopen(sys.argv[1]).parse()
    m = MPA(d)
    m.calculate()

After the calculate() method is called, the following attributes are available:

* ``aoresults``: a three dimensional array with spin, molecular orbital, and atomic orbitals as the axes, so that ``aoresults[0][45][0]`` gives the contribution of the 1st atomic orbital to the 46th alpha/restricted molecular orbital,
* ``fragresults``: a three dimensional array with spin, molecular orbital, and atoms as the axes, so that ``fragresults[1][23][4]`` gives the contribution of the 5th fragment orbitals to the 24th beta molecular orbital)
* ``fragcharges``: a vector with the number of (partial) electrons in each fragment, so that ``fragcharges[2]`` gives the number of electrons in the 3rd fragment.

Custom fragments
~~~~~~~~~~~~~~~~

The calculate method chooses atoms as the fragments by default, and optionally accepts a list of lists containing the atomic orbital numbers (e.g. ``[[0, 1, 2], [3, 4, 5, 6], ...]``) of arbitrary fragments. Calling it in this way is useful if one is more interested in the contributions of groups of atoms or even certain orbitals or orbital groups, such as metal d, to the molecular orbitals. In this case, fragresults and fragcharges reflect the chosen groups of atomic orbitals instead of atoms.

Custom progress
~~~~~~~~~~~~~~~

The Mulliken class also can take a progress class as an argument so that the progress of the calculation can be monitored:

.. code-block:: python

    from cclib.method import MPA
    from cclib.parser import ccopen
    from cclib.progress import TextProgress
    import logging

    progress = TextProgress()
    d = ccopen("mycalc.out", logging.ERROR).parse(progress)

    m = MPA(d, progress, logging.ERROR)
    m.calculate()

.. index::
    single: methods; Löwdin Population Analysis

Löwdin Population Analysis
--------------------------

The LPA class available from cclib.method performs Löwdin population analysis and can be used as follows:

.. code-block:: python

    import sys

    from cclib.method import LPA
    from cclib.parser import ccopen

    d = ccopen(sys.argv[1]).parse()
    m = LPA(d)
    m.calculate()

.. index::
    single: methods; Bickelhaupt Population Analysis

Bickelhaupt Population Analysis
-------------------------------

The Bickelhaupt class available from cclib.method performs Bickelhaupt population analysis that has been proposed in *Organometallics* 1996, 15, 13, 2923–2931. `doi:10.1021/om950966x <https://pubs.acs.org/doi/abs/10.1021/om950966x>`_

The contribution of the a-th atomic orbital to the i-th molecular orbital in this method is written in terms of the molecular orbital coefficients, c, and the overlap matrix, S:

.. math:: \Phi_{ai,\alpha} = \sum_b w_{ab,\alpha} c_{ai,\alpha} c_{bi,\alpha} S_{ab}

where the weights :math:`w_{ab}` that are applied on the Mulliken atomic orbital contributions are defined as following when the coefficients of the molecular orbitals are substituted into equation 11 in the original article.

.. math:: w_{ab,\alpha} = 2 \frac{\sum_k c_{ak,\alpha}^2}{\sum_i c_{ai,\alpha}^2 + \sum_j c_{bj,\alpha}^2}

In case of unrestricted calculations, :math:`\alpha` charges and :math:`\beta` charges are each determined to obtain total charge. In restricted calculations, :math:`\alpha` subscript can be ignored since the coefficients are equivalent for both spin orbitals.

The weights are introduced to replace the somewhat arbitrary partitioning of off-diagonal charges in the Mulliken population analysis, which divides the off-diagonal charges identically to both atoms. Bickelhaupt population analysis instead divides the off-diagonal elements based on the relative magnitude of diagonal elements.

.. code-block:: python

    import sys

    from cclib.method import Bickelhaupt
    from cclib.parser import ccopen

    d = ccopen(sys.argv[1]).parse()
    m = Bickelhaupt(d)
    m.calculate()

After the calculate() method is called, the following attributes are available:

* ``aoresults``: a three dimensional array with spin, molecular orbital, and atomic orbitals as the axes, so that ``aoresults[0][45][0]`` gives the contribution of the 1st atomic orbital to the 46th alpha/restricted molecular orbital,
* ``fragresults``: a three dimensional array with spin, molecular orbital, and atoms as the axes, so that ``fragresults[1][23][4]`` gives the contribution of the 5th fragment orbitals to the 24th beta molecular orbital)
* ``fragcharges``: a vector with the number of (partial) electrons in each fragment, so that ``fragcharges[2]`` gives the number of electrons in the 3rd fragment.

Custom fragments
~~~~~~~~~~~~~~~~

The calculate method chooses atoms as the fragments by default, and optionally accepts a list of lists containing the atomic orbital numbers (e.g. ``[[0, 1, 2], [3, 4, 5, 6], ...]``) of arbitrary fragments. Calling it in this way is useful if one is more interested in the contributions of groups of atoms or even certain orbitals or orbital groups, such as metal d, to the molecular orbitals. In this case, fragresults and fragcharges reflect the chosen groups of atomic orbitals instead of atoms.

Custom progress
~~~~~~~~~~~~~~~

The Bickelhaupt class also can take a progress class as an argument so that the progress of the calculation can be monitored:

.. code-block:: python

    from cclib.method import Bickelhaupt
    from cclib.parser import ccopen
    from cclib.progress import TextProgress
    import logging

    progress = TextProgress()
    d = ccopen("mycalc.out", logging.ERROR).parse(progress)

    m = Bickelhaupt(d, progress, logging.ERROR)
    m.calculate()

..
   Overlap Population Analysis
   ---------------------------

Density Matrix calculation
--------------------------

The Density class from cclib.method can be used to calculate the density matrix:

.. code-block:: python

    from cclib.parser import ccopen
    from cclib.method import Density

    parser = ccopen("myfile.out")
    data = parser.parse()

    d = Density(data)
    d.calculate()

After ``calculate()`` is called, the density attribute is available. It is simply a NumPy array with three axes. The first axis is for the spin contributions, and the second and third axes are for the density matrix, which follows the standard definition.

Mayer's Bond Orders
-------------------

This method calculates the Mayer's bond orders for a given molecule:

.. code-block:: python

    import sys

    from cclib.parser import ccopen
    from cclib.method import MBO

    parser = ccopen(sys.argv[1])
    data = parser.parse()

    d = MBO(data)
    d.calculate()

After ``calculate()`` is called, the fragresults attribute is available, which is a NumPy array of rank 3. The first axis is for contributions of each spin to the MBO, while the second and third correspond to the indices of the atoms.

Charge Decomposition Analysis
-----------------------------

The Charge Decomposition Analysis (CDA) as developed by Gernot Frenking et al. is used to study the donor-acceptor interactions of a molecule in terms of two user-specified fragments.

The CDA class available from cclib.method performs this analysis:

.. code-block:: python

    from cclib.io import ccopen
    from cclib.method import CDA

    molecule = ccopen("molecule.log")
    frag1 = ccopen("fragment1.log")
    frag2 = ccopen("fragment2.log")

    # if using CDA from an interactive session, it's best
    # to parse the files at the same time in case they aren't
    # parsed immediately---go get a drink!

    m = molecule.parse()
    f1 = frag1.parse()
    f2 = frag2.parse()

    cda = CDA(m)
    cda.calculate([f1, f2])

After ``calculate()`` finishes, there should be the donations, bdonations (back donation), and repulsions attributes to the cda instance. These attributes are simply lists of 1-dimensional NumPy arrays corresponding to the restricted or alpha/beta molecular orbitals of the entire molecule. Additionally, the CDA method involves transforming the atomic basis functions of the molecule into a basis using the molecular orbitals of the fragments so the attributes mocoeffs and fooverlaps are created and can be used in population analyses such as Mulliken or C-squared (see Fragment Analysis for more details).

There is also a script provided by cclib that performs the CDA from a command-line::

    $ cda molecule.log fragment1.log fragment2.log
    Charge decomposition analysis of molecule.log

     MO#      d       b       r
    -----------------------------
       1:  -0.000  -0.000  -0.000
       2:  -0.000   0.002   0.000
       3:  -0.001  -0.000   0.000
       4:  -0.001  -0.026  -0.006
       5:  -0.006   0.082   0.230
       6:  -0.040   0.075   0.214
       7:   0.001  -0.001   0.022
       8:   0.001  -0.001   0.022
       9:   0.054   0.342  -0.740
      10:   0.087  -0.001  -0.039
      11:   0.087  -0.001  -0.039
    ------ HOMO - LUMO gap ------
      12:   0.000   0.000   0.000
      13:   0.000   0.000   0.000
    ......

Notes
~~~~~

* Only molecular orbitals with non-zero occupancy will have a non-zero value.
* The absolute values of the calculated terms have no physical meaning and only the relative magnitudes, especially for the donation and back donation terms, are of any real value (Frenking, et al.)
* The atom coordinates in molecules and fragments must be the same, which is usually accomplished with an argument in the QM program (the NoSymm keyword in Gaussian, for instance).
* The current implementation has some subtle differences than the code from the Frenking group. The CDA class in cclib follows the formula outlined in one of Frenking's CDA papers, but contains an extra factor of 2 to give results that agree with those from the original CDA program. It also doesn't include negligible terms (on the order of 10^-6) that result from overlap between MOs on the same fragment that appears to be included in the Frenking code. Contact atenderholt (at) gmail (dot) com for discussion and more information.

.. index::
    single: methods; Bader's QTAIM

Bader's QTAIM
----------------

Bader's QTAIM charges define the border between the atoms in the molecule as a surface where each point on the surface has zero flux. In other words, the points on the surface of the Bader spaces satisfy the equation :math:`\nabla \rho (r) \cdot n(r) = 0`. In cclib, numerical calculation of QTAIM charges through the algorithm proposed in `[Henkelman2006] <http://theory.cm.utexas.edu/henkelman/code/bader/download/henkelman06_354.pdf>`_  is implemented.

Calculating the Bader's QTAIM charges in cclib follow similar steps as other population analysis methods. The following code provides an example of how QTAIM charges can be obtained.

.. code-block:: python

    from cclib.method import Bader
    from cclib.method import Volume
    from cclib.parser import ccopen
    from cclib.progress import TextProgress
    import logging

    progress = TextProgress()
    d = ccopen("mycalc.out", logging.ERROR).parse(progress)

    # Inputs for Volume object below are origin, top corner, and spacing
    # represented in Cartesian coordinates and in angstroms.
    vol = Volume((-3, -3, -3), (3, 3, 3), (.1, .1, .1))

    m = Bader(d, vol)
    m.calculate()

After the calculate() method is called, the following attributes are available:

* ``fragresults``: a three dimensional array x, y, and z coordinates from the Volume object as the axes, so that ``fragresults[1][2][3]`` gives the Bader space (in integers starting from 1) that the grid space in (0, 1, 2) position belongs to.
* ``matches``: a vector with the Bader space (integers starting from 1) that an atom is matched with.
* ``fragcharges``: a vector with the number of (partial) electrons in each atom, so that ``fragcharges[2]`` gives the number of electrons in the 3rd atom.

Since some computational chemistry packages support writing out charge densities as cube files during calculations, it is highly recommended to do so especially for larger systems.
To calculate Bader charges from a cube file, a ``Volume`` object should be prepared by reading in a cube file and should be passed into a Bader object as shown below:

.. code-block:: python

    from cclib.method import volume
    from cclib.method import Bader
    from cclib.parser import ccopen
    from cclib.progress import TextProgress
    import logging

    progress = TextProgress()
    d = ccopen("mycalc.out", logging.ERROR).parse(progress)

    # Read in from cube file
    vol = volume.read_from_cube("mycalc.cube")

    m = Bader(d, vol)
    m.calculate()


.. index::
    single: methods; DDEC6

DDEC6
----------------

DDEC6 is a Stockholder-like charge partitioning method introduced in `[Manz2016] <https://doi.org/10.1039/C6RA04656H>`_. Proatom densities should be provided for DDEC6 method. Proatom densities generated [#proatom]_ using `horton`_ can be passed on as an argument for the constructor of the DDEC6 object. The DDEC6 algorithm requires many numerical integrations so a fine grid is necessary to obtain accurate results. Calculating the electron density on a fine grid using a Volume object is slow, therefore we recommend that electron densities are imported from cube files.

Because a lot of numerical integrations are present in DDEC6 algorithm, fine grid is necessary to obtain satisfying results.

.. _`Horton`: http://theochem.github.io/horton/2.1.1/


.. code-block:: python

    from cclib.method import DDEC6
    from cclib.method import Volume
    from cclib.parser import ccopen
    from cclib.progress import TextProgress
    import logging

    progress = TextProgress()
    d = ccopen("mycalc.out", logging.ERROR).parse(progress)

    # Inputs for Volume object below are origin, top corner, and spacing
    # represented in Cartesian coordinates and in angstroms.
    vol = Volume((-3, -3, -3), (3, 3, 3), (.1, .1, .1))

    # Alternatively, read in from cube file
    vol = volume.read_from_cube("mycalc.cube")

    m = DDEC6(d, vol, '/path/to/horton_proatom_density_directory')
    m.calculate()

Third argument to the constructor of DDEC6 object points to the directory proatom densities are stored. The proatom densities can be generated by using `horton`_ . Follow the steps described in horton `documentation`_ for its ``horton-atomdb.py`` command.

.. _`horton`: http://theochem.github.io/horton/2.1.1/
.. _`documentation`: http://theochem.github.io/horton/2.1.1/user_postproc_aim.html#horton-atomdb-py-build-a-pro-atom-database

After the calculate() method is called, the following attributes are available:

* ``fragcharges``: a vector with the number of (partial) electrons in each atom, so that ``fragcharges[2]`` gives the number of electrons in the 3rd atom.
* ``refcharges``: a two dimensional array where the first index indicates whether the reference charges are from first step of DDEC6 algorithm or the second step. Second index refers to the atoms that the charges are assigned to and follows the same order as the order used by the input ``ccData`` object.

.. index::
    single: methods; Hirshfeld Population Analysis

Hirshfeld Population Analysis
-----------------------------

Hirshfeld Population Analysis is the most basic charge partitioning method among stockholder-type methods and was introduced initially in `[Hirshfeld1977] <https://doi.org/10.1007/BF00549096>`_.
In Hirshfeld method, proatom densities, which are charge densities of neutral atoms that comprise the given molecule, are used to calculate the weights that will be applied to partition charge densities on each grid point.
Proatom densities can be generated [#proatom]_ using `horton`_ and can be passed on as an argument for the constructor of the Hirshfeld object. For Hirshfeld calculations, it is recommended that electron densities are imported from cube files.

.. _`Horton`: http://theochem.github.io/horton/2.1.1/


.. code-block:: python

    from cclib.method import Hirshfeld
    from cclib.method import Volume
    from cclib.parser import ccopen
    from cclib.progress import TextProgress
    import logging

    progress = TextProgress()
    d = ccopen("mycalc.out", logging.ERROR).parse(progress)

    # Inputs for Volume object below are origin, top corner, and spacing
    # represented in Cartesian coordinates and in angstroms.
    vol = Volume((-3, -3, -3), (3, 3, 3), (.1, .1, .1))

    # Alternatively, read in from cube file
    vol = volume.read_from_cube("mycalc.cube")

    m = Hirshfeld(d, vol, '/path/to/horton_proatom_density_directory')
    m.calculate()

Third argument to the constructor of Hirshfeld object points to the directory proatom densities are stored. The proatom densities can be generated by using `horton`_ . Follow the steps described in horton `documentation`_ for its ``horton-atomdb.py`` command.

.. _`horton`: http://theochem.github.io/horton/2.1.1/
.. _`documentation`: http://theochem.github.io/horton/2.1.1/user_postproc_aim.html#horton-atomdb-py-build-a-pro-atom-database

After the calculate() method is called, the following attributes are available:

* ``fragcharges``: a vector with the number of (partial) electrons in each atom, so that ``fragcharges[2]`` gives the number of electrons in the 3rd atom.

.. index::
    single: methods; Accessing additional methods through bridge


Nuclear Properties
------------------

This method calculates various properties of nuclei based on the parsed data.

Calculate the nuclear repulsion energy using the snippet below:

.. code-block:: python

    import cclib
    from cclib.method import Nuclear

    parser = cclib.io.ccopen("test.log")
    data = parser.parse()
    nuclear = Nuclear(data)

    nre = nuclear.repulsion_energy()
    print(nre)


Accessing additional methods through bridge
-------------------------------------------

Some other population analyses methods including Hirshfeld partial charges and Iterative Stockholder charges can be calculated using bridge functions implemented in cclib. For more information, refer to `bridge`_ section of the documentation.

.. _`bridge`: bridge.html

..
   Electron Density Calculation
   ----------------------------


.. rubric:: Footnotes

.. [#proatom] To generate proatom densities from horton using Psi4, after generating Psi4 input files, add ``mv *.molden atom.default.molden`` on line 25 of ``run_psi4.sh`` before executing the script.
