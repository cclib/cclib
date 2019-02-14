===========
Development
===========

Basic instructions
==================

The default cclib files distributed with a release, as described in `How to install`_, do not include any unit tests and logfiles necessary to run those tests. This section covers how to download the full source along with all test data and scripts, and how to use these for development and testing.

.. _`How to install`: how_to_install.html

Cloning cclib from GitHub
~~~~~~~~~~~~~~~~~~~~~~~~~

cclib is hosted by the fantastic people at `GitHub`_ (previously at `Sourceforge`_) in a `git`_ repository. You can download a `zipped archive`_ of the current development version (called `master`) for installation and testing or browse the available `releases`_. In order to contribute any changes, however, you will need to create a local copy of the repository:

.. code-block:: bash

    git clone https://github.com/cclib/cclib.git

.. _`GitHub`: https://github.com
.. _`Sourceforge`: https://sourceforge.net
.. _`git`: https://git-scm.com
.. _`zipped archive`: https://github.com/cclib/cclib/archive/master.zip
.. _`releases`: https://github.com/cclib/cclib/releases

Guidelines
~~~~~~~~~~~~~~~~

We follow a typical GitHub collaborative model, relying on `forks and pull requests`_. In short, the development process consists of:

* `Creating your own fork`_ of cclib in order to develop
* `Creating a pull request`_ to contribute your changes
* Reviewing and merging open pull requests (by someone else)
* Using `issues`_ to plan and prioritize future work

.. _`creating your own fork`: https://help.github.com/articles/fork-a-repo
.. _`creating a pull request`: https://help.github.com/articles/creating-a-pull-request
.. _`forks and pull requests`: https://help.github.com/articles/using-pull-requests
.. _`issues`: https://github.com/cclib/cclib/issues

Here are some general guidelines for developers who are contributing code:

* Run and review the unit tests (see below) before submitting a pull request.
* There should normally not be more failed tests than before your changes.
* For larger changes or features that take some time to implement, `using branches`_ is recommended.

.. _`using branches`: https://help.github.com/articles/branching-out

Releasing a new version
~~~~~~~~~~~~~~~~~~~~~~~

The release cycle of cclib is irregular, with new versions being created as deemed necessary after significant changes or new features. We roughly follow semantic versioning with respect to the `parsed attributes`_.

When creating a new release on GitHub, the typical procedure might include the following steps:

* Update the `CHANGELOG`_, `ANNOUNCE`_ and any other files that might change content with the new version
* Make sure that `setup.py`_ has the right version number, as well as __version__ in `__init__.py`_ and any other relevant files
* Update the download and install instructions in the documentation, if appropriate
* Create a branch for the release, so that development can continue
* Run all tests for a final time and fix any remaining issues
* Tag the release (make sure to use an annotated tag using ``git -a``) and upload it (``git push --tags``)
* Run `manifest.py`_ to update the MANIFEST file
* Create the source distributions (``python setup.py sdist --formats=gztar,zip``) and Windows binary installers (``python setup.py bdist_wininst``)
* Create a release on GitHub using the created tag (see `Creating releases`_) and upload the source distributions and Windows binaries
* Email the users and developers mailing list with the message in `ANNOUNCE`_
* Update the Python package index (https://pypi.python.org/pypi/cclib), normally done by ``python setup.py register``
* For significant releases, if appropriate, send an email to the `CCL list`_ and any mailing lists for computational chemistry packages supported by cclib

.. _`parsed attributes`: data.html

.. _`ANNOUNCE`: https://github.com/cclib/cclib/blob/master/ANNOUNCE
.. _`CHANGELOG`: https://github.com/cclib/cclib/blob/master/CHANGELOG
.. _`setup.py`: https://github.com/cclib/cclib/blob/master/setup.py
.. _`__init__.py`: https://github.com/cclib/cclib/blob/master/src/cclib/__init__.py
.. _`manifest.py`: https://github.com/cclib/cclib/blob/master/manifest.py

.. _`Creating releases`: https://help.github.com/articles/creating-releases

.. _`CCL list`: http://www.ccl.net

Testing
=======

.. index::
    single: testing; unit tests

The `test directory`_, which is not included in the default download, contains the test scripts that keep cclib reliable, and keep the developers sane. With any new commit or pull request to cclib on GitHub the tests are triggered and run with `Travis CI`_, for both the current production version |release| (|travis_prod|) as well as master (|travis_master|).

The input files for tests, which are logfiles from computational chemistry programs, are located in the `data directory`_. These are a central part of cclib, and any progress should always be supported by corresponding tests. When a user opens an issue or reports a bug, it is prudent to write a test that reproduces the bug as well as fixing it. This ensures it will remain fixed in the future. Likewise, extending the coverage of data attributes to more programs should proceed in parallel with the growth of unit tests.

.. _`Travis CI`: https://travis-ci.org/cclib/cclib

.. |travis_prod| image:: https://travis-ci.org/cclib/cclib.svg?branch=v1.6.1
.. |travis_master| image:: https://travis-ci.org/cclib/cclib.svg?branch=master

.. _`data directory`: https://github.com/cclib/cclib/tree/master/data
.. _`test directory`: https://github.com/cclib/cclib/tree/master/test

.. index::
    single: testing; unit tests

Unit tests
~~~~~~~~~~

Unit tests check that the parsers work correctly for typical calculation types on small molecules, usually water or 1,4-divinylbenzene (dvb) with :math:`C_{\mathrm{2h}}` symmetry. The corresponding logfiles stored in folders like ``data/NWChem/basicNWChem6.0`` are intended to test logfiles for an approximate major version of a program, and are standardized for all supported programs to the extent possible. They are located alongside the code in the repository, but are not normally distributed with the source. Attributes are considered supported only if they are checked by at least one test, and the `table of attribute coverage`_ is generated automatically using this criterion.

The job types currently included as unit tests:

* restricted and unrestricted single point energies for dvb (RHF/STO-3G **and** B3LYP/STO-3G)
* geometry optimization and scan for dvb (RHF/STO-3G and/or B3LYP/STO-3G)
* frequency calculation with IR and Raman intensities for dvb (RHF/STO-3G or B3LYP/STO-3G)
* single point energy for carbon atom using a large basis set such as aug-cc-pCVQZ
* Møller–Plesset and coupled cluster energies for water (STO-3G or 6-31G basis set)
* static polarizabilities for tryptophan (RHF/STO-3G)

.. _`table of attribute coverage`: data_dev.html#details-of-current-implementation

Adding a new program version
----------------------------

There are a few conventions when adding a new supported program version to the unit tests:
* Two different recent versions are typically used in the unit tests. If there already are two, move the older version(s) the regression suite (see below).
* When adding files for the new version, first copy the corresponding files for the last version already in cclib. Afterwards, check in files from the new program version as changes to the copied files. This procedure makes it easy to look at the differences introduced with the new version in git clients.

.. index::
    single: testing; regressions

Regression tests
~~~~~~~~~~~~~~~~

Regression tests ensure that bugs, once fixed, stay fixed. These are real-life files that at some point broke a cclib parser, and are stored in folders like ``data/regression/Jaguar/Jaguar6.4``. The files associated with regression tests are not stored stored together with the source code as they are often quite large. A separate repository on GitHub, `cclib-data`_, is used to track these files, and we do not distribute them with any releases.

For every bug found in the parsers, there should be a corresponding regression test that tests this bug stays fixed. The process is automated by `regression.py`_, which runs through all of our test data, both the basic data and regression files, opens them, tries to parse, and runs any relevant regression tests defined for that file. New regression tests are added by creating a function ``testMyFileName_out`` according to the examples at the start of `regression.py`_.

Using both the unit and regression tests, the line-by-line `test coverage`_ shows which parts of cclib are touched by at least one test. When adding new features and tests, the Travis CI `testing script`_ can be run locally to generate the HTML coverage pages and ensure that the tests exercise the feature code.

.. _`cclib-data`: https://github.com/cclib/cclib-data
.. _`regression.py`: https://github.com/cclib/cclib/blob/master/test/regression.py

.. _`test coverage`: coverage/index.html
.. _`testing script`: https://github.com/cclib/cclib/blob/master/travis/run_pytest.sh

Websites related to cclib
=========================

* The official `cclib organization on github`_
* The `cclib project page on Sourceforge`_ (inactive now)
* The `cclib page for Travis CI`_
* The `cclib entry on PyPI`_
* The `cclib entry on Ohloh`_

.. _`cclib organization on github`: https://github.com/cclib
.. _`cclib project page on Sourceforge`: http://sourceforge.net/projects/cclib/
.. _`cclib entry on PyPI`: http://www.python.org/pypi/cclib
.. _`cclib page for Travis CI`: https://travis-ci.org/cclib/cclib
.. _`cclib entry on Ohloh`: https://www.ohloh.net/p/cclib

Developers
==========

Besides input from a number of people `listed in the repository`_, the following developers have contributed code to cclib (in alphabetical order):

* `Eric Berquist`_
* `Karol M. Langner`_
* `Noel O'Boyle`_
* Christopher Rowley
* Adam Tenderholt

.. _`listed in the repository`: https://github.com/cclib/cclib/blob/master/THANKS

.. _`Eric Berquist`: https://github.com/berquist
.. _`Karol M. Langner`: https://github.com/langner
.. _`Noel O'Boyle`: https://www.redbrick.dcu.ie/~noel/
