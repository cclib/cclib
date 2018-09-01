.. index::
    module: data_notes

Parsed data notes
=================

This is a list of descriptions and notes for all the data attributes currently parsed by cclib, either in the official release (|release|) or development branch. In particular, this page contains technical details about the interpretation of attributes, how to produce them in the various programs and examples in some cases. For a summary and details of the current implementation by the different parsers, please see the `extracted data`_ page and its `development`_ version.

.. _`extracted data`: data.html
.. _`development`: data_dev.html

aonames
-------

This attribute contains the atomic orbital names. These are not normalised as the following examples show, although a reasonable attempt is made to get them close to each other. Users will need to know what each orbital is by knowing the basis set inside out, rather than relying on this data. Such is life, as GAMESS does not provide enough information.

* Gaussian gives names of the form::

    ['C1_1S', 'C1_2S', 'C1_2PX', 'C1_2PY', 'C1_2PZ', 'C2_1S', 'C2_2S', 'C2_2PX', 'C2_2PY', 'C2_2PZ', 'C3_1S', 'C3_2S', 'C3_2PX', 'C3_2PY', 'C3_2PZ', 'C4_1S', 'C4_2S', 'C4_2PX', 'C4_2PY', 'C4_2PZ', 'C5_1S', 'C5_2S', 'C5_2PX', 'C5_2PY', 'C5_2PZ', 'H6_1S', 'H7_1S', 'H8_1S', 'C9_1S', 'C9_2S', 'C9_2PX', 'C9_2PY', 'C9_2PZ', 'C10_1S', 'C10_2S', 'C10_2PX', 'C10_2PY', 'C10_2PZ', 'H11_1S', 'H12_1S', 'H13_1S', 'C14_1S', 'C14_2S', 'C14_2PX', 'C14_2PY', 'C14_2PZ', 'H15_1S', 'C16_1S', 'C16_2S', 'C16_2PX', 'C16_2PY', 'C16_2PZ', 'H17_1S', 'H18_1S', 'C19_1S', 'C19_2S', 'C19_2PX', 'C19_2PY', 'C19_2PZ', 'H20_1S']

* GAMESS gives names of the form::

    ['C1_1S', 'C1_2S', 'C1_3X', 'C1_3Y', 'C1_3Z', 'C2_1S', 'C2_2S', 'C2_3X', 'C2_3Y', 'C2_3Z', 'C3_1S', 'C3_2S', 'C3_3X', 'C3_3Y', 'C3_3Z', 'C4_1S', 'C4_2S', 'C4_3X', 'C4_3Y', 'C4_3Z', 'C5_1S', 'C5_2S', 'C5_3X', 'C5_3Y', 'C5_3Z', 'C6_1S', 'C6_2S', 'C6_3X', 'C6_3Y', 'C6_3Z', 'H7_1S', 'H8_1S', 'H9_1S', 'H10_1S', 'C11_1S', 'C11_2S', 'C11_3X', 'C11_3Y', 'C11_3Z', 'C12_1S', 'C12_2S', 'C12_3X', 'C12_3Y', 'C12_3Z', 'H13_1S', 'H14_1S', 'C15_1S', 'C15_2S', 'C15_3X', 'C15_3Y', 'C15_3Z', 'C16_1S', 'C16_2S', 'C16_3X', 'C16_3Y', 'C16_3Z', 'H17_1S', 'H18_1S', 'H19_1S', 'H20_1S']

And for a large basis set calculation on a single C atom:

* Gaussian::

    ['C1_1S', 'C1_2S', 'C1_3S', 'C1_4S', 'C1_5S', 'C1_6PX', 'C1_6PY', 'C1_6PZ', 'C1_7PX', 'C1_7PY', 'C1_7PZ', 'C1_8PX', 'C1_8PY', 'C1_8PZ', 'C1_9PX', 'C1_9PY', 'C1_9PZ', 'C1_10D 0', 'C1_10D+1', 'C1_10D-1', 'C1_10D+2', 'C1_10D-2', 'C1_11D 0', 'C1_11D+1', 'C1_11D-1', 'C1_11D+2', 'C1_11D-2', 'C1_12D 0', 'C1_12D+1', 'C1_12D-1', 'C1_12D+2', 'C1_12D-2', 'C1_13F 0', 'C1_13F+1', 'C1_13F-1', 'C1_13F+2', 'C1_13F-2', 'C1_13F+3', 'C1_13F-3', 'C1_14F 0', 'C1_14F+1', 'C1_14F-1', 'C1_14F+2', 'C1_14F-2', 'C1_14F+3', 'C1_14F-3', 'C1_15G 0', 'C1_15G+1', 'C1_15G-1', 'C1_15G+2', 'C1_15G-2', 'C1_15G+3', 'C1_15G-3', 'C1_15G+4', 'C1_15G-4', 'C1_16S', 'C1_17PX', 'C1_17PY', 'C1_17PZ', 'C1_18D 0', 'C1_18D+1', 'C1_18D-1', 'C1_18D+2', 'C1_18D-2', 'C1_19F 0', 'C1_19F+1', 'C1_19F-1', 'C1_19F+2', 'C1_19F-2', 'C1_19F+3', 'C1_19F-3', 'C1_20G 0', 'C1_20G+1', 'C1_20G-1', 'C1_20G+2', 'C1_20G-2', 'C1_20G+3', 'C1_20G-3', 'C1_20G+4', 'C1_20G-4']

* GAMESS::

    ['C1_1S', 'C1_2S', 'C1_3S', 'C1_4S', 'C1_5S', 'C1_6X', 'C1_6Y', 'C1_6Z', 'C1_7X', 'C1_7Y', 'C1_7Z', 'C1_8X', 'C1_8Y', 'C1_8Z', 'C1_9X', 'C1_9Y', 'C1_9Z', 'C1_10XX', 'C1_10YY', 'C1_10ZZ', 'C1_10XY', 'C1_10XZ', 'C1_10YZ', 'C1_11XX', 'C1_11YY', 'C1_11ZZ', 'C1_11XY', 'C1_11XZ', 'C1_11YZ', 'C1_12XX', 'C1_12YY', 'C1_12ZZ', 'C1_12XY', 'C1_12XZ', 'C1_12YZ', 'C1_13XXX', 'C1_13YYY', 'C1_13ZZZ', 'C1_13XXY','C1_13XXZ', 'C1_13YYX', 'C1_13YYZ', 'C1_13ZZX', 'C1_13ZZY', 'C1_13XYZ', 'C1_14XXX', 'C1_14YYY', 'C1_14ZZZ', 'C1_14XXY', 'C1_14XXZ', 'C1_14YYX', 'C1_14YYZ', 'C1_14ZZX', 'C1_14ZZY', 'C1_14XYZ', 'C1_15XXXX', 'C1_15YYYY', 'C1_15ZZZZ', 'C1_15XXXY', 'C1_15XXXZ', 'C1_15YYYX', 'C1_15YYYZ', 'C1_15ZZZX', 'C1_15ZZZY', 'C1_15XXYY', 'C1_15XXZZ', 'C1_15YYZZ', 'C1_15XXYZ', 'C1_15YYXZ', 'C1_15ZZXY', 'C1_16S', 'C1_17S', 'C1_18S', 'C1_19X', 'C1_19Y', 'C1_19Z', 'C1_20X', 'C1_20Y', 'C1_20Z', 'C1_21X', 'C1_21Y', 'C1_21Z', 'C1_22XX', 'C1_22YY', 'C1_22ZZ', 'C1_22XY', 'C1_22XZ', 'C1_22YZ', 'C1_23XX', 'C1_23YY', 'C1_23ZZ', 'C1_23XY', 'C1_23XZ', 'C1_23YZ', 'C1_24XXX', 'C1_24YYY', 'C1_24ZZZ', 'C1_24XXY', 'C1_24XXZ', 'C1_24YYX', 'C1_24YYZ', 'C1_24ZZX', 'C1_24ZZY', 'C1_24XYZ', 'C1_25S', 'C1_26X', 'C1_26Y', 'C1_26Z', 'C1_27XX', 'C1_27YY', 'C1_27ZZ', 'C1_27XY', 'C1_27XZ', 'C1_27YZ', 'C1_28XXX', 'C1_28YYY', 'C1_28ZZZ', 'C1_28XXY', 'C1_28XXZ', 'C1_28YYX', 'C1_28YYZ', 'C1_28ZZX', 'C1_28ZZY', 'C1_28XYZ', 'C1_29XXXX', 'C1_29YYYY', 'C1_29ZZZZ', 'C1_29XXXY', 'C1_29XXXZ', 'C1_29YYYX', 'C1_29YYYZ', 'C1_29ZZZX', 'C1_29ZZZY', 'C1_29XXYY', 'C1_29XXZZ', 'C1_29YYZZ', 'C1_29XXYZ', 'C1_29YYXZ', 'C1_29ZZXY']

aooverlaps
----------

This is a 2-dimensional array which holds the numerical values of the overlap between basis functions (also called atomic orbitals). It is needed for most analyses like `Mulliken`_, `C squared`_, and `Mayer's Bond Orders`_. The indices of the matrix correspond to the basis functions of interest. This matrix is symmetric, so ``aooverlaps[i,j]`` is the same as ``aooverlaps[j,i]``.

Some examples:

* ``aooverlaps[0,3]`` is the overlap between the 1st and 4th basis function
* ``aooverlaps[2,:]`` is a 1-dimensional array containing the overlap between every basis function and the 3rd basis function

**ADF**: not present by default, printed when `PRINT Smat` is in the input; do not mistake with `fooverlaps`_.

**DALTON**: no option to print as of version 2013.

**Gaussian**: ``iop(3/33=1)`` must be specified in the input file.

.. _`Mulliken`: methods.html#mulliken-population-analysis-mpa
.. _`C squared`: methods.html#c-squared-population-analysis-cspa
.. _`Mayer's Bond Orders`: methods.html#mayer-s-bond-orders

atombasis
---------

The attribute ``atombasis`` is a list, each element being a list that contains the atomic orbital indices on the respective atom. For example, ``atombasis[1]`` will contain the indices of atomic orbitals on the second atom of the molecule.

.. index::
    single: properties; atomcharges (attribute)

atomcharges
-----------

The attribute ``atomcharges`` contains the atomic partial charges as taken from the output file. Since these charges are arbitrary and depend on the details of a population analysis, this attribute is dictionary containing any number of various atomic charges. The keys in this dictionary are strings naming the population analysis, and the values are arrays of rank 1 and contain the actual charges.

Currently, cclib parses Mulliken, Löwdin, NPA and CHELPG charges, whose respective dictionary keys are ``mulliken``, ``lowdin``, ``natural`` and ``chelpg``.

In practice, these may differ somewhat from the values cclib calculates in the various `calculation methods`_.

**Molpro**: use the ``pop`` command (see http://www.molpro.net/info/2015.1/doc/manual/node515.html).

.. _`calculation methods`: methods.html

atomcoords
----------

The attribute ``atomcoords`` contains the atomic coordinates as taken from the output file. This is an array of rank 3, with a shape (n,m,3) where n is 1 for a single point calculation and >=1 for a geometry optimisation and m is the number of atoms.

**Gaussian**: for geometry optimisations, the "Standard orientation" sections are extracted.

**Molpro**: typically prints output about geometry optimisation in a separate logfile. So, both that and the initial output need to be passed to the cclib parser.

atommasses
----------

The attribute ``atommasses`` contains the masses of all atoms in unified atomic mass units, or Daltons (Da). This is an array or rank 1.

atomnos
-------

An array of integers for the atomic numbers, or the number of protons in the atom nuclei.

atomspins
---------

The attribute ``atomspins`` contains the atomic spin densities as calculated in a population analysis and taken from the output file. Since these densities are arbitrary and depend on the particular population analysis, this attribute is dictionary. In analogy to `atomcharges`_, the keys in this dictionary are strings naming the population analysis, and the values are arrays of rank 1 and contain the actual spin densities.

Currently, cclib parses Mulliken and Löwdin spin densities, whose respective dictionary keys are ``mulliken`` and ``lowdin``.

.. index::
    single: energy; ccenergies (attribute)

ccenergies
----------

A one-dimensional array holds the total molecule energies including Coupled Cluster corrections. The array's length is 1 for single point calculations and larger for optimisations. Only the highest theory level is parsed into this attribute (for example, CCSD energies as opposed to CCD energies, or CCSD(T) as opposed to CCSD energies).

charge
------

Net charge of the calculated system, in units of ``e``.

coreelectrons
-------------

The attribute ``coreelectrons`` contains the number of core electrons in each atom's pseudopotentials. It is an array of rank 1, with as many integer elements as there are atoms.

etenergies
----------

This is a rank 1 array that contains the energies of electronic transitions from a reference state to the excited states of the molecule, in ``cm<sup>-1</sup>``. There should be as many elements to this array as there are excited states calculated. Any type of excited state calculation should provide output that can be parsed into this attribute.

etoscs
------

The attribute ``etoscs`` is a rank 1 array that contains the oscillator strengths of transitions from the reference (ground) state to the excited electronic states of the of the molecule. As for `etenergies`_ and other attributes related to excited states, there should as many elements in this array as there are excited states in the calculation.

etsecs
------

The singly-excited configurations that contribute to electronic transitions are stored in ``etsecs``. It is a list (for each electronic transition from the reference ground state) of lists (for each singly-excited configuration) with three members each:

 * a tuple (moindex, alpha/beta), which indicates the MO where the transition begins
 * a tuple (moindex, alpha/beta), which indicates the MO where the transition ends
 * a float (which can be negative), the coefficient of this singly-excited configuration

In these tuples, the value of alpha/beta is 0 or 1, respectively. For a restricted calculation, this value is always 0, although some programs (GAMESS) sometimes print coefficients for both alpha and beta electrons.

The excitation coefficient is always converted to its unnormalized value by cclib - so the sum of the squared coefficients of all alpha and beta excitations should be unity. It is important to keep in mind, however, that only the square of the excitation coefficient has a physical meaning, and its sign depends on the numerical procedures used by each program.

etsyms
------

The attributes ``etsyms`` is a list containing the symmetries (strings) of the excited states found in the calculation. As for `etenergies`_ and other attributes related to excited states, there should be as many elements in this list as there are excited states in the calculation.

Note that while the symmetry descriptions start with the string ``Singlet`` or ``Triplet``, the exact format differs between programs.

fonames
-------

ADF uses symmetry-adapted fragment orbitals (SFOs) as its basis. These SFOs are generally orthonormal linear combinations of atomic orbitals. This makes it difficult to determine which individual atomic orbitals form the basis in calculations that have any symmetry. In addition, ADF allows "fragment" calculations which use the molecular orbitals of the fragments (FOs, or fragment orbitals) for building up the calculated molecular orbitals.

The difficulty in handling the basis for a molecule with symmetry and the availability of extra information in the fragment calculations makes using `aonames`_ (as specified for the other formats) inappropriate, except for certain circumstances. Therefore, an extra member called fonames is available for the adfparser.

Some examples:

``C1+C4_1S+1S`` - Orbitals from carbon 1 and carbon 4 can interact, and their ``1S`` orbitals mix in a positive manner

``C1+C4_1Px-1Px`` - Orbitals from carbon 1 and carbon 4 can interact, and their ``1Px`` orbitals mix in a negative manner

``bdt1_37A`` - Molecular orbital 37A from the fragment bdt1

**ADF**: There are no required inputfile options for fonames to be supported; however, if one wishes to have SFOs map directly to atomic basis functions, there are two requirements. First, the ``Symmetry NOSYM`` option must be given to force ADF to not linearly combine atomic orbitals into SFOs. Second, fragment calculations cannot be done (for obvious reasons). Also, it is suggested that ``Eigval 99999 99999`` be put into an ``Eprint`` block of the input file of a spin-restricted calculation so that every molecular orbital energy will be printed.

fooverlaps
----------

This is a 2-dimensional array that holds numerical values for the spacial overlap between basis functions. It is very similar to `aooverlaps`_, but differs because of the way ADF performs the calculation (see below for more details). The matrix indices correspond to the fragment orbitals; see the examples listed for `aonames`_.

**Background**

ADF uses symmetry-adapted fragment orbitals (SFOs) as its basis. These SFOs are generally orthonormal linear combinations of atomic orbitals. This makes it difficult to determine which individual atomic orbitals form the basis in calculations that have any symmetry. In addition, ADF allows "fragment" calculations which use the molecular orbitals of the fragments (FOs, or fragment orbitals) for building up the calculated molecular orbitals.

The difficulty in handling the basis for a molecule with symmetry and the availability of extra information in the fragment calculations makes using aooverlaps (as specified for the other formats) inappropriate, except for certain circumstances. Therefore, an extra member called fooverlaps is available for the ADF parser.

**ADF**: There are no required inputfile options for fooverlaps to be supported; however, if one wishes to have SFOs map directly to atomic basis functions, there are two requirements. First, the ``Symmetry NOSYM`` option must be given to force ADF to not linearly combine atomic orbitals into SFOs. Second, fragment calculations cannot be done (for obvious reasons). Also, it is suggested that ``Eigval 99999 99999`` be put into an ``Eprint`` block of the input file of a spin-restricted calculation so that every molecular orbital energy will be printed.

.. index::
    single: basis sets; gbasis (attribute)

gbasis
------

This attribute stores information about the Gaussian basis functions that were used in the calculation, per atom using the same conventions as `PyQuante <http://pyquante.sf.net>`_. Specifically, ``gbasis`` is a list of lists iterating over atoms and Gaussian basis functions. The elements (basis functions) are tuples of length 2 consisting of orbital type (e.g. 'S', 'P' or 'D') and a list (per contracted GTO) of tuples of size 2 consisting of the exponent and coefficient. Confused? Well, here's ``gbasis`` for a molecule consisting of a single C atom with a STO-3G basis:

.. code-block:: python

    [ # per atom
        [
            ('S', [
                (71.616837, 0.154329),
                (13.045096, 0.535328),
                (3.530512, 0.444635),
                ]),
            ('S', [
                (2.941249, -0.099967),
                (0.683483, 0.399513),
                (0.222290, 0.700115),
                ]),
            ('P', [
                (2.941249, 0.155916),
                (0.683483, 0.607684),
                (0.222290, 0.391957),
                ]),
        ]
    ]

For D and F functions there is an important distinction between pure (5D, 7F) or Cartesian (6D, 10F) functions. PyQuante can only handle Cartesian functions, but we should extract this information in any case, and perhaps work to extend the PyQuante basis set format to include this.

**Gaussian**: the `GFINPUT`_ keyword should normally be used (`GFPRINT`_ gives equivalent information in a different format and is supported in cclib after v1.2).

**GAMESS/GAMESS-UK**: no special keywords are required, but the basis is only available for symmetry inequivalent atoms. There does not seem to be any way to get GAMESS to say which atoms are related through symmetry. As a result, if you want to get basis set info for every atom, you need to reduce the symmetry to C1.

**Jaguar**: for more information see manual (for example at http://yfaat.ch.huji.ac.il/jaguar-help/mand.html#114223)

**ORCA**: include ``Print[ P_Basis ] 2`` in the ``output`` block

.. _`GFINPUT`: http://www.gaussian.com/g_tech/g_ur/k_gfinput.htm
.. _`GFPRINT`: http://www.gaussian.com/g_tech/g_ur/k_gfprint.htm

.. index::
    single: geometry optimisation; geotargets (attribute)

geotargets
----------

Geotargets are the target values of the criteria used to determine whether a geometry optimisation has converged. The targets are stored in an array of length ``n``, where ``n`` is the number of targets, and the actual values of these criteria are stored for every optimisation step in the attribute `geovalues`_. Note that cclib does not carry information about the meaning of these criteria, and it is up to the user to interpret the values properly for a particular program. Below we provide some details for several parsers, but it is always a good idea to refer to the source documentation.

In some special cases, the values in ``geotargets`` will be `numpy.inf`_.

**GAMESS UK**: the criteria used for geometry convergence are based on the ``TOL`` parameter, which can be set using the ``XTOLL`` directive. The fault value of this parameter and the conditions required for convergence vary among the various optimisation strategies (see the `GAMESS-UK manual section on controlling optimisation`_ for details). In ``OPTIMIZE`` mode, ``TOL`` defaults to 0.003 and the conditions are,

    - maximum change in variables below TOL,
    - average change in variables smaller than TOL * 2/3,
    - maximum gradient below TOL * 1/4,
    - average gradient below TOL * 1/6.

.. _`GAMESS-UK manual section on controlling optimisation`: http://www.cfs.dl.ac.uk/docs/html/part4/node14.html 

**Jaguar** has several geometry convergence criteria,

    * gconv1: maximum element of gradient (4.5E-04)
    * gconv2: rms of gradient elements (3.0E-04)
    * gconv5: maximum element of nuclear displacement (1.8E-03)
    * gconv6: rms of nuclear displacement elements (1.2E-03)
    * gconv7: difference between final energies from previous and current geometry optimisation iterations (5.0E-05)

Note that a value for gconv7 is not available until the second iteration, so it is set to zero in the first element of `geovalues`_.

**Molpro** has custom convergence criteria, as described in the `manual <Molpro manual convergence_>`_:

    The standard MOLPRO convergence criterion requires the maximum component of the gradient to be less then :math:`3 \cdot 10^{-4}` [a.u.] and the maximum energy change to be less than :math:`1 \cdot 10^{-6}` [H] or the maximum component of the gradient to be less then $ 3 \cdot 10^{-4}$ [a.u.] and the maximum component of the step to be less then :math:`3 \cdot 10^{-4}` [a.u.].

.. _Molpro manual convergence: https://www.molpro.net/info/2012.1/doc/manual/node592.html

**ORCA** tracks the change in energy as well as RMS and maximum gradients and displacements. As of version 3.0, an optimisation is considered converged when all the tolerances are met, and there are four exceptions:

    * the energy is within 25x the tolerance and all other criteria are met
    * the gradients are overachieved (1/3 of the tolerance) and displacements are reasonable (at most 3x the tolerance)
    * the displacements are overachieved (1/3 of the tolerance) and the gradients are reasonable (at most 3x the tolerance)
    * the energy gradients and internal coordinates are converged (bond distances, angles, dihedrals and impropers)

**Psi** normally tracks five different values, as described `in the documentation`_, but their use various depending on the strategy employed. The default strategy (QCHEM) check whether the maximum force is converged and if the maximum energy change or displacement is converged. Additionally, to aid with flat potential energy surfaces, convergence is as assumed when the root mean square force converged to 0.01 of its default target. Note that Psi print values even for targets that are not being used -- in these cases the targets are parsed as `numpy.inf`_ so that they can still be used (any value will be converged).

.. _`in the documentation`: http://sirius.chem.vt.edu/psi4manual/latest/optking.html

.. _`numpy.inf`: http://docs.scipy.org/doc/numpy-1.8.1/user/misc.html#ieee-754-floating-point-special-values

.. index::
    single: geomtry optimisation; geovalues (attribute)

geovalues
---------

These are the current values for the criteria used to determine whether a geometry has converged in the course of a geometry optimisation. It is an array of dimensions ``m x n``, where ``m`` is the number of geometry optimisation iterations and ``n`` the number of target criteria.

Note that many programs print atomic coordinates before and after a geometry optimisation, which means that there will not necessarily be ``m`` elements in `atomcoords`_.

If the optimisation has finished successfully, the values in the last row should be smaller than the values in geotargets_ (unless the convergence criteria require otherwise).

hessian
-------

An array of rank 1 that contains the elements of the `hessian <http://en.wikipedia.org/wiki/Hessian_matrix>`_ or force constant matrix. Only the lower triangular part of the 3Nx3N matrix is stored (this may change in the future, maybe also only the unweighted matrix will be parsed).

.. index::
    single: molecular orbitals; homos (attribute)

homos
-----

A 1D array that holds the indexes of the highest occupied molecular orbitals (HOMOs), which contains one element for restricted and two elements for unrestricted calculations. These indexes can be applied to other attributes describing molecular orbitals, such as `moenergies`_ and `mocoeffs`_.

.. index::
    single: molecular orbitals; mocoeffs (attribute)

metadata
--------

A dictionary containing metadata_ (data about data) for the calculation. Currently, it can contain the following possible attributes, not all of which are implemented for each parser.

* ``basis_set``: A string with the name of the basis set, if it is printed anywhere as a standard name.
* ``coord_type``: For the ``coords`` field, a string for the representation of stored coordinates. Currently, it is one of ``xyz``, ``int``/``internal``, or ``gzmat``.
* ``coords``: A list of lists with shape ``[natoms, 4]`` which contains the input coordinates (those found in the input file). The first column is the atomic symbol as a string, and the next three columns are floats. This is useful as many programs reorient coordinates for symmetry reasons.
* ``functional``: A string with the name of the density functional used.
* ``info``: A list of strings, each of which is an information or log message produced during a calculation.
* ``input_file_contents``: A string containing the entire input file, if it is echoed back during the calculation.
* ``input_file_name``: A string containing the name of the input file, with file extension. It may not contain the entire path to the file.
* ``keywords``: A list of strings corresponding to the keywords used in the input file, in the loose format used by ORCA.
* ``methods``: A list of strings containing each method used in order. Currently, the list may contain ``HF``, ``DFT``, ``LMP2``/``DF-MP2``/``MP2``, ``MP3``, ``MP4``, ``CCSD``, and/or ``CCSD(T)``/``CCSD-T``.
* ``package``: A string with the name of the quantum chemistry program used.
* ``package_version``: A string representation of the package version. It is formatted to allow comparison using relational operators.
* ``success``: A boolean for whether or not the calculation completed properly.
* ``unrestricted``: A boolean for whether or not the calculation was performed with a unrestricted wavefunction.
* ``warnings``: A list of strings, each of which is a warning produced during a calculation.

The implementation and coverage of metadata is currently inconsistent. In the future, metadata may receive its own page similar to `extracted data`_.

.. _metadata: https://en.wikipedia.org/wiki/Metadata

mocoeffs
--------

A list of rank 2 arrays containing the molecular orbital (MO) coefficients. The list is of length 1 for restricted calculations, but length 2 for unrestricted calculations. For the array(s) in the list, the first axis corresponds to molecular orbitals, and the second corresponds to basis functions.

Examples:

* ``mocoeffs[0][2,5]`` -- The coefficient of the 6th basis function of the 3rd alpha molecular orbital
* ``mocoeffs[1][:,0]`` -- An array of the 1st basis function coefficients for the every beta molecular orbital

Note: For restricted calculation, ``mocoeffs`` is still a list, but it only contains a single rank 2 array so you access the matrix with mocoeffs[0].

**GAMESS-UK** - the `FORMAT HIGH`_ directive needs to be included if you want information on all of the eigenvalues to be available. In versions before 8.0 for unrestricted calculations, ``FORMAT HIGH`` does not increase the number of orbitals for which the molecular orbital coefficents are printed, so that there may be more orbital information on the alpha orbitals compared to the beta orbitals, and as a result the extra beta molecular orbital coefficients for which information is not available will be padded out with zeros by cclib.

**Molpro** - does not print MO coefficients at all by default, and you must add in the input ``GPRINT,ORBITALS``. What's more, this prints only the occupied orbitals, and to get virtuals add also ``ORBPTIN,NVIRT``, where ``NVIRT`` is how many virtuals to print (can be a large number like 99999 to print all).

.. index::
    single: molecular orbitals; moenergies (attribute)

moenergies
----------

A list of rank 1 arrays containing the molecular orbital energies in eV. The list is of length 1 for restricted calculations, but length 2 for unrestricted calculations.

**GAMESS-UK**: similar to `mocoeffs`_, the directive `FORMAT HIGH`_ needs to be used if you want all of the eigenvalues printed.

**Jaguar**: the first ten virtual orbitals are printed by default. In order to print more, use the ``ipvirt`` keyword, with ``ipvirt=-1`` printing all virtual orbitals.

.. _`FORMAT HIGH`: http://www.cfs.dl.ac.uk/docs/html/part3/node8.html#SECTION00083000000000000000

.. index::
    single: properties; moments (attribute)

moments
-------

This attribute contains the dipole moment vector and any higher electrostatic multipole moments for the whole molecule. It comprises a list of one dimensional arrays,

* the first is the reference point used in the multipole expansion, which is normally the center of mass,
* the second is the dipole moment vector, in Debyes (:math:`\mathbf{\mathrm{D}}`),
* the third array contains the raw molecular quadrupole moments in lexicographical order, that is the XX, XY, XZ, YY, YZ and ZZ moments, in Buckinghams (:math:`\mathbf{\mathrm{B}}`),
* any further arrays contain the raw molecular multipole moments of higher rank, in lexicographical order and in units of :math:`\mathbf{\mathrm{D}} \cdot Å^{L-1} = 10^{-10} \mathrm{esu} \cdot Å^L`

Note that by default cclib will provide the last moments printed, if several are printed in the course of a geometry optimisation or other job type involving several more than one geometry. For post-Hartree-Fock calculations, such as MP2 or coupled cluster, the uncorrelated moments are reported if none are printed for the final wavefunction.

.. index::
    single: molecular orbitals; mosyms (attribute)

mosyms
------

For unrestricted calculations, this is a list of two lists containing alpha and beta symmetries (i.e. ``[[alpha_syms],[beta_syms]]``) containing strings for the orbital symmetries, arranged in order of energy. In a restricted calculation, there is only one nested list (``[[syms]]``).

The symmetry labels are normalised and cclib reports standard symmetry names:

    ======= ======= ======= ==========  ==================          ======
    cclib   ADF     GAMESS  GAMESS-UK   Gaussian                    Jaguar
    ======= ======= ======= ==========  ==================          ======
    A       A       A       a           A                           A
    A1      A1      A1      a1          A1                          A1
    Ag      A.g     AG      ag          AG                          Ag
    A'      AA      A'      a'          A'                          Ap
    A"      AAA     A' '    a" or a' '  A"                          App
    A1'     AA1     A1'     a1'         A1'                         A1p
    A1"     AAA1    A1"     a1"         A1"                         A1pp
    sigma   Sigma                       SG
    pi      Pi                          PI
    phi     Phi                         PHI (inferred)
    delta   Delta                       DLTA but DLTU/DLTG
    sigma.g Sigma.g                     SGG
    ======= ======= ======= ==========  ==================          ======

* ADF - the full list can be found [http://www.scm.com/Doc/Doc2005.01/ADF/ADFUsersGuide/page339.html here].
* GAMESS-UK - to get the list, 'grep "data yr" input.m' if you have access to the source. Note that for E, it's split into "e1+" and "e1-" for instance.
* Jaguar - to get the list, look at the examples in schrodinger/jaguar-whatever/samples if you have access to Jaguar. Note that for E, it's written as E1pp/Ap, for instance.
* NWChem - if molecular symmetry is turned off or set to C1, symmetry adaption for orbitals is also deactivated, and can be explicitly turned on with `adapt on` in the SCF block

Developers:

* The use of a function with doctests for each of these cases is recommended, to make sure that the conversion is robust. There is a prototype called normalisesym() in logfileparser.py which should be overwritten in the subclasses if necessary (there is a unittest to make sure that this has been done).
* The character tables `here <http://www.mpip-mainz.mpg.de/~gelessus/group.html>`_ may be useful in determining the correspondence between the labels used by the comp chem package and the commonly-used symbols.

.. index::
    single: energy; mpenergies (attribute)

mpenergies
----------

The attribute ``mpenergies`` holds the total molecule energies including Møller-Plesset correlation energy corrections in a two-dimensional array. The array's shape is (n,L), where ``n`` is 1 for single point calculations and larger for optimisations, and ``L`` is the order at which the correction is truncated. The order of elements is ascending, so a single point MP5 calculation will yield mpenergies as :math:`E_{MP2}, E_{MP3}, E_{MP4}, E_{MP5}`.

**ADF**: does not perform such calculations.

**GAMESS**: second-order corrections (MP2) are available in GAMESS-US, and MP2 through MP3 calculations in PC-GAMESS (use ``mplevl=n`` in the ``$contrl`` section).

**GAMESS-UK**: MP2 through MP3 corrections are available.

**Gaussian**: MP2 through MP5 energies are available using the ``MP`` keyword. For MP4 corrections, the energy with the most substitutions is used (SDTQ by default).

**Jaguar**: the LMP2 is available.

mult
----

The attribute ``mult`` is an integer and represents the spin multiplicity of the calculated system, which in turn is the total spin plus one.

natom
-----

``Natom`` is an integer, the number of atoms treated in the calculation.

.. index::
    single: basis sets; nbasis (attribute)

nbasis
------

An integer representing the number of basis functions used in the calculation.

.. index::
    single: basis sets; nmo (attribute)

nmo
---

The number of molecular orbitals in the calculation. It is an integer and is typically equal to `nbasis`_, but may be less than this if a linear dependency was identified between the basis functions.

Commands to get information on all orbitals:

**GAMESS-UK**: only usually prints information on the 5 lowest virtual orbitals. "FORMAT HIGH" should make it do this for all of the orbitals, although GAMESS-UK 7.0 has a bug that means that this only works for restricted calculations.

**Jaguar**: the first ten virtual orbitals are printed by default; in order to print more of them, use the ``ipvirt`` keyword in the input file, with ``ipvirt=-1`` printing all virtual orbitals (see the `manual <Jaguar manual nmo_>`_ for more information).

.. _Jaguar manual nmo: http://www.pdc.kth.se/doc/jaguar4.1/html/manual/mang.html#644675

optdone
-------

Flags whether a geometry optimisation has completed. Currently this attribute is a single Boolean value, which is set to True when the final `atomcoords`_ represent a converged geometry optimisation. In the future, ``optdone`` will be a list that indexes which elements of `atomcoords`_ represent converged geometries. This functionality can be used starting from version 1.3, from the command line by passing the ``--future`` option to ``ccget``,

.. code-block:: bash

    $ ccget optdone data/Gaussian/basicGaussian09/dvb_gopt.out
    Attempting to parse data/Gaussian/basicGaussian09/dvb_gopt.out
    optdone:
    True

    $ ccget --future optdone data/Gaussian/basicGaussian09/dvb_gopt.out
    Attempting to parse data/Gaussian/basicGaussian09/dvb_gopt.out
    optdone:
    [4]

or by providing the corresponding argument to ``ccopen``,

.. code-block:: python

    from cclib.parser import ccopen
    parser = ccopen("filename", optdone_as_list=True) # could also do future=True instead of optdone_as_list
    data = parser.parse()

scfenergies
-----------

An array containing the converged SCF energies of the calculation, in eV. For an optimisation log file, there will be as many elements in this array as there were optimisation steps.

**Molpro**: typically prints output about geometry optimisation in a separate logfile. So, both that and the initial output need to be passed to the cclib parser.

scftargets
----------

Target thresholds for determining whether the current SCF run has converged, stored in a ``n x m`` array, where ``n`` is the number of geometry optimisation steps (1 for a single point calculation) and ``m`` is the number of criteria. The criteria vary between programs, and depending on the program they may be constant for the whole of a geometry optimisation or they may change between optimisation steps. A more detailed description for each program follows.

**ADF**: There are two convergence criteria which are controlled by ``SCFcnv`` in the `CONVERGE subkey of the SCF block`_.

* The maximum element of the commutator of the Fock matrix and P-matrix needs to be below ``SCFcnv``.
* The norm of the same matrix needs to be below ``10*SCFcnv``.

This hard target is normally used for single point calculations and the last step of geometry optimisations, and it defaults to 1.0E-6. There is also a soft target ``scfconv2`` that defaults to 1.0E-3, which can be switched on and is used by ADF automatically in some cases such as the first step in a geometry optimization.

For intermediate steps in a geometry optimisation the situation is more complicated and depends on the gradient and the integration accuracy. A post on the ADF user's forum revealed that it is calculated as follows:

.. math:: \mathrm{new\,criteria} = max( \mathrm{SCFcnv}, \, min(\mathrm{old\,criteria}, \, \mathrm{grdmax}/30, 10^{-\mathrm{accint}})) ),

where ``old criteria`` is the initial value or from the previous geometry cycle, ``grdmax`` is the maximum gradient from the last geometry step and ``accint`` is the current integration accuracy.

.. _`CONVERGE subkey of the SCF block`: http://www.scm.com/Doc/Doc2014/ADF/ADFUsersGuide/page235.html#keyscheme%20INTEGRATION

**GAMESS**: Two criteria are, the maximum and root-mean-square (RMS) density matrix change, are used with a default starting value of 5.0E-05. It seems these values can change over the course of a geometry optimisation. ROHF calculations use SQCDF instead of the standard RMS change.

**GAMESS-UK**: According to `the manual <GAMESS-UK manual convergence_>`_, convergence is determined by convergence of density matrix elements. The default value for SCF is 1E-5, but it appears to be 1E-7 for geoopts.

.. _`GAMESS-UK manual convergence`: http://www.cfs.dl.ac.uk/docs/html/part4/node6.html

**Gaussian**: normally three criteria are used.

* The RMS change in the density matrix elements, with a default of 1.0E-4 (1.0E-8 for geo opts).
* Maximum change in the density matrix elements, with a default of 1.0E-2 (1.0E-6 for geo opts).
* The change in energy, with a default threshold of 5.0E-05 (1.0E-06 for geo opts).

**Jaguar 4.2**: The targets in Jaguar 4.2 (based on the manual) depend on whether the job is a geometry optimisation or not. For geometry optimisations and hyper/polarisability calculation, the RMS change in the density matrix elements is used as a criterion (controlled by the ``dconv`` keyword), with a default of 5.0E6.
The energy convergence criterion (keyword ``econv``) is ignored for geometry optimisation calculations but is used for SCF calculations, and the default in this case is 5.0E5, except for hyper/polarisability calcualtions where it is 1.0E6.

scfvalues
---------

The attribute ``scfvalues`` is a list of arrays of dimension ``n x m`` (one element for each step in a geometry optimisation), where ``n`` is the number of SCF cycles required for convergence and ``m`` is the number of SCF convergence target criteria. For some packages, you may need to include a directive to make sure that SCF convergence information is printed to the log file

**Gaussian**: requires the `route section`_ to start with #P

.. _`route section`: http://www.gaussian.com/g_tech/g_ur/k_route.htm

**GAMESS-UK**: convergence information is printed only for the first optimisation step by default, but can be forced at all steps by adding ``IPRINT SCF`` to the input file.

vibdisps
--------

The attribute ``vibdisps`` stores the Cartesian displacement vectors from the output of a vibrational frequency calculation. It is a rank 3 array having dimensions ``M x N x 3``, where ``M`` is the number of normal modes and ``N`` is the number of atoms. ``M`` is typically ``3N-6`` (``3N-5`` for linear molecules).
