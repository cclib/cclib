How to install
==============

This page describes how to download, install and use the basic functionality of cclib.

Requirements
------------

Before you install cclib, you need to make sure that you have the following:
 * Python (at least version 3.4 is recommended, although 2.7 is still tested)
 * NumPy (at least version 1.5 is recommended)

Python is an open-source programming language available from https://www.python.org. It is available for Windows as well as being included in most Linux distributions. In Debian/Ubuntu it is installed as follows (as root):

.. code-block:: bash

    apt-get install python python-dev

NumPy (Numerical Python) adds a fast array facility and linear algebra routines to Python. It is available from https://www.numpy.org. Windows users should use the most recent NumPy installation for the Python version they have (e.g. numpy-1.0.3.1.win32-py2.4.exe for Python 2.4). Linux users are recommended to find a binary package for their distribution rather than compiling it themselves. In Debian/Ubuntu it is installed as follows (as root):

.. code-block:: bash

    apt-get install python-numpy

To test whether Python is on the ``PATH``, open a command prompt window and type:

.. code-block:: bash

    python

If Python is not on the ``PATH`` and you use Windows, add the full path to the directory containing it to the end of the ``PATH`` variable under Control Panel/System/Advanced Settings/Environment Variables. If you use Linux and Python is not on the ``PATH``, put/edit the appropriate line in your .bashrc or similar startup file.

To test that NumPy is working, try importing it at the Python prompt. You should see something similar to the following::

    $ python
    Python 3.7.0 (default, Jul 15 2018, 10:44:58)
    [GCC 8.1.1 20180531] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import numpy as np
    >>> np.__version__
    '1.15.0'
    >>>

(To exit, press CTRL+D in Linux or CTRL+Z,Enter in Windows)

Installing using pip
--------------------

pip_ is a command-line tool that often comes pre-packaged with Python distributions which can download packages from the `Python Package Index`_ (PyPI). To see if it's installed, on Linux or macOS try::

    $ which pip

and on Windows::

    cmd> where.exe pip

The version of `cclib uploaded to PyPI`_ can then be installed globally using::

    python -m pip install cclib

or to your home directory using::

    python -m pip install --user cclib

.. _pip: https://pip.pypa.io/en/stable/
.. _`Python Package Index`: https://pypi.org/
.. _`cclib uploaded to PyPI`: https://pypi.python.org/pypi/cclib

Installing using a system package manager
-----------------------------------------

If you're using `Debian GNU/Linux`_, `Ubuntu`_, or a similar distribution, there are official `cclib packages`_ that you can install in any package manager (should as synaptic) or with one simple command:

.. code-block:: bash

    aptitude install cclib

There are in fact two packages, `python-cclib`_ containing the Python module, and `cclib`_ which installs just the user scripts. If you also need to also install the unittests and logfiles we distribute, you will need to install the `cclib-data`_ package from the non-free repositories (due to license issues). Because of distribution release cycles, package manager versions of cclib may be out of date compared to the PyPI version.

.. _`Debian GNU/Linux`: https://www.debian.org
.. _`Ubuntu`: https://www.ubuntu.com
.. _`cclib packages`: https://packages.debian.org/src:cclib
.. _`python-cclib`: https://packages.debian.org/stretch/python-cclib
.. _`cclib`: https://packages.debian.org/stretch/cclib
.. _`cclib-data`: https://packages.debian.org/stretch/cclib-data

Manual download and install
---------------------------

The source code of the newest release of cclib (version |release|) is distributed as:
 * A .zip file: https://github.com/cclib/cclib/releases/download/v1.6/cclib-1.6.1.zip
 * A .tar.gz file: https://github.com/cclib/cclib/releases/download/v1.6/cclib-1.6.1.tar.gz
 * Windows binary installers (see the `newest release page`_)

On Windows, if you choose to download the .exe files instead, you can install simply by double-clicking on the file. To uninstall, use the "Add and Remove Programs" menu in the Control Panel.

None of these files include the tests and logfiles used for testing. In order to download all tests, we also provide source archives on the `newest release page`_.

If you are using the .zip or .tar.gz files, extract the contents of the file at an appropriate location, which we will call INSTALLDIR. Open a command prompt and change directory to INSTALLDIR. Next, run the following commands to install cclib:

.. code-block:: bash

    python setup.py build
    python setup.py install # (as root)

or, if pip_ is available::

    python -m pip install .

To test, trying importing '''cclib''' at the Python prompt. You should see something similar to the following::

    $ python
    Python 3.7.0 (default, Jul 15 2018, 10:44:58)
    [GCC 8.1.1 20180531] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import cclib
    >>> cclib.__version__
    '1.6.1'
    >>>

.. _`newest release page`: https://github.com/cclib/cclib/releases/tag/v1.3.1

What next?
----------

 * Read the list and specifications of the `parsed data`_ and related `data notes`_
 * Test the program using the test data files included in the full source distribution
 * Run the unit and regression tests in the test directory (``testall.py`` and ``regression.py``)
 * Send any questions to the cclib-users mailing list at https://lists.sourceforge.net/lists/listinfo/cclib-users.
 * Write some computational chemistry algorithms using information parsed from cclib and donate the code to the project

.. _`parsed data`: data.html
.. _`data notes`: data_notes.html
