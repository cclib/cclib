=========================
Development (sourceforge)
=========================

This page contains old development content, from when cclib was hosted at Sourceforge, and which has become obsolete for the current version of cclib (|release|).

Websites related to cclib
=========================

* The official `cclib repository on github`_
* The `cclib project page on Sourceforge`_ (inactive now)
* The `cclib entry on Freecode`_
* The `cclib entry on Ohloh`_
* The `cclib entry on Cheeseshop`_

.. _`cclib repository on github`: https://github.com/cclib/cclib
.. _`cclib project page on Sourceforge`: http://sourceforge.net/projects/cclib/
.. _`cclib entry on Freecode`: http://freecode.com/projects/cclib
.. _`cclib entry on Ohloh`: https://www.ohloh.net/p/cclib
.. _`cclib entry on Cheeseshop`: http://www.python.org/pypi/cclib

Instructions for developers
===========================

Checking cclib out of subversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cclib is hosted by the fantastic people at Sourceforge. The development version of cclib is stored in a Subversion (SVN) repository at Sourceforge (you can read some background material on this at [1]). The following command will download the current development version into a folder called cclib:

.. code-block:: bash

    svn checkout https://cclib.svn.sourceforge.net/svnroot/cclib/trunk cclib

To install, follow the guidelines described on the Install page.

Releasing a version
~~~~~~~~~~~~~~~~~~~

* Update the CHANGELOG and ANNOUNCE
* Make sure that setup.py has the right version number, as well as __version__ in __init__.py.
* Run manifest.py to update the MANIFEST if necessary.
* Do a final merge of the trunk to release branch
* Create the source distributions 

.. code-block:: bash

    python setup.py sdist --formats=gztar,zip upload
    python setup.py bdist_wininst # rename from cclib-0.x.win32.exe to cclib-0.x-py2.3.exe
    python2.4 setup.py bdist_wininst                       # rename to cclib-0.x-py2.4.exe

* Run the tests for a final time after removing cclib (rm -rf $PYTHONDIR/Lib/site-packages/cclib) and reinstalling from the source distribution.
* Tag the release.
* Use releaseforge to make a new release
* Update the download instructions on the wiki
* Wait 24 hours (for the sourceforge mirrors to get a copy)
* Email the users and developers mailing list with the ANNOUNCE
* Create a news item and copy and paste ANNOUNCE into it
* Update http://www.freshmeat.net (click "Add release" - summarise the changelog in changes, and give the link to the changelog in the appropriate box)
* Update the Python cheeseshop 

.. code-block:: bash

    python setup.py register

* For a major release, if appropriate, send an email to the CCL list, the GAMESS users list and the ADF users list. 

Creating a logfile release
~~~~~~~~~~~~~~~~~~~~~~~~~~

Although all of the regression log files are stored on the web server, for the convenience of users we make periodic releases of all of these as one large .tar file (no point in .gz as they are already zipped). To create a release: 

.. code-block:: bash

    cd data
    tar -c -f logfiles.tar -T regressionfiles.txt

* Create a new release in package logfiles called logfiles-rnnn, where nnn is the current SVN revision.
* You can use a variation of the following text in the release notes: "A release of test logfiles for use with rnnn of the cclib code." 

Copying a logfile to sourceforge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use scp or winscp to connect to "USERNAME,cclib@web.sf.net". The logfiles are in htdocs/data. For interactive shell access, you will need to first type "ssh -t USERNAME,cclib@shell.sf.net create". After about 10 seconds, your own private shell will be created at SourceFroge which you can then log into. 

Release branches
~~~~~~~~~~~~~~~~

Historically, branches were managed by hand. Now I've now started using svnmerge.py, which keeps track of what's already been merged, and what's already been marked as 'not to merge'.

* Create a release branch for cclib-1.0.1 

.. code-block:: bash

    svn copy https://cclib.svn.sourceforge.net/svnroot/cclib/trunk https://cclib.svn.sourceforge.net/svnroot/cclib/branches/cclib-1.0.1

* Check it out 

.. code-block:: bash

    svn checkout https://svn.sourceforge.net/svnroot/branches/cclib-1.0.1 branchcclib101

* Initialise merge tracking 

.. code-block:: bash

    cd branchcclib101
    python svnmerge.py init ../trunk
    svn commit -F svnmerge-commit-message.txt

* Merge some stuff 

.. code-block:: bash

    python ..\svnmerge.py avail -S /trunk  # or add "--log"
    python ..\svnmerge.py merge -S /trunk # Merge all changes
    svn commit -F svnmerge-commit-message.txt

Source code upload policy
~~~~~~~~~~~~~~~~~~~~~~~~~

As a sort of guide for developers who are commiting source code revisions to the SVN repository, we recommend the following:

* Run the tests before commiting (at least testall.py)
* If tests that previously passed now no longer do (when we have a more complete and stable release, this will read "If any tests fail"), and you don't have time to fix things before commiting, commits your changes to a branch as follows: 

.. code-block:: bash

    svn copy https://svn.sourceforge.net/svnroot/cclib/trunk https://..../branches/brokenadfparser
           -m "Informative log message about why you're branching"
    # change directory into the 'top' of your working copy where setup.py is
    svn switch https://svn.sourceforge.net/svnroot/cclib/branches/brokenadfparser
    # (note that this preserves the local modifications, only now these modifications
       are to the branch instead of the trunk)
    svn commit # (commits the local modifications to the branch)

* As soon as the tests that previously passed are passed again, merge the changes and remove the branch (this should be within a revision or two of the branch): 

.. code-block:: bash

    svn log --stop-on-copy # (to find the revision when the branching took place)
    svn switch https://...../mytrunk
    svn merge --dry-run -r123:HEAD https://..../mybranch
    svn merge -r123:HEAD https://..../mybranch
    svn commit -m "Merging 123:HEAD of mybranch into trunk"
    svn remove https://..../mybranch

Ensuring source code quality
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To ensure source code quality, or at least consistency, we use the following tools/metrics:

* `Cheesecake index`

Testing
=======

.. index::
    single: testing; unit tests

Unit tests
~~~~~~~~~~

In order to check whether our parsers extract information in the correct format, with the correct units, we have unit tests that parse a series of basic data files (see below) of the same calculation undertaken with different programs. Running testall.py in the test directory runs the whole test suite, but it is also possible to individually run the tests for GeoOpts (testGeoOpt.py), Single Point calculations, and so on.

Note that no change should be commited to the repository if it increases the number of failed tests (unless you are adding new tests, of course). 

.. index::
    single: testing; regressions

Regression tests
~~~~~~~~~~~~~~~~

Regression tests ensure that bugs, once fixed, stay fixed. That is, for every bug found in our parsers, we should add a regression test and then fix the bug. This process is simplified by regression.py in the test directory.

regression.py runs through all of our test data, both the basic data and the real life log files, uses 'ccopen()' to guess its type (checks for mistakes) and open it, parses it (catches any errors), and runs any relevant regression tests (catches any failures).

New regression tests are added by creating a function testMyFileName_out following the examples at the start of *regression.py*.

Test data
~~~~~~~~~

.. index::
    single: testing; test data

The test directory (not included in the release version at the moment) contains all of the tests that help keep cclib working, and keep us sane. In general we use two types of data files for testing:

1. 'basic' data files are stored in folders like "basicJaguar6.4" and are b3lyp/sto-3g calculations on 1,4-divinylbenzene (dvb) with C2h symmetry. These jobs (a geometry optimisation, a single point calculation (one restricted and another unrestricted), frequency calculation, a TD-DFT calculation, and any variants of these which break the parser such as symmetry/nosymmetry) are run for each parser. These data files are stored in SVN and may be included in future releases.
2. real-life parser-breaking files are stored in folders like "Jaguar6.4". These data files are *not* stored in SVN as they are often massive but are stored on the web server and downloaded using a shell script (contained in the data directory). These files are also available as a download from the File Release page on Sourceforge (only updated every so often). 

Doc tests
~~~~~~~~~

Doc tests are a nice Python invention for unit testing individual functions. To run the doctests in a particular file, you need to run the script. For example, "python gaussianparser.py" runs the doctests contained in gaussianparser.py. To run all of the doctests at once, you need to install a testing tool like nose, and then use the following command (note that many errors may be due to missing libraries like BioPython):

.. code-block: bash

    > "C:\Program Files\Python24\Scripts\nosetests.exe" cclib --with-doctest -e test* -v
    ERROR
    ERROR
    Doctest: cclib.bridge.cclib2openbabel.makeopenbabel ... ok
    ERROR
    ERROR
    Doctest: cclib.parser.adfparser.ADF.normalisesym ... ok
    Doctest: cclib.parser.gamessparser.GAMESS.normalise_aonames ... ok
    Doctest: cclib.parser.gamessparser.GAMESS.normalisesym ... ok
    Doctest: cclib.parser.gamessukparser.GAMESSUK.normalisesym ... ok
    Doctest: cclib.parser.gaussianparser.Gaussian.normalisesym ... ok
    Doctest: cclib.parser.jaguarparser.Jaguar.normalisesym ... ok
    Doctest: cclib.parser.logfileparser.Logfile.float ... ok
    Doctest: cclib.parser.utils.PeriodicTable ... ok
    Doctest: cclib.parser.utils.convertor ... ok
    ERROR
    ERROR
    ......

Other useful pages
==================

* The `methods in the development version`_
* The `parsed data in the development version`_
* The `progress page`_, which describes where we are and what we are trying to do

.. _`methods in the development version`: methods_dev.html
.. _`parsed data in the development version`: data_dev.html
.. _`progress page`: progress.html
