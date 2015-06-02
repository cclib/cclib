#!/usr/bin/env python3
#
# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""cclib: parsers and algorithms for computational chemistry

cclib is a Python library that provides parsers for computational
chemistry log files. It also provides a platform to implement
algorithms in a package-independent manner.
"""

import sys


# Chosen from http://www.python.org/pypi?:action=list_classifiers
classifiers = """Development Status :: 5 - Production/Stable
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)
Natural Language :: English
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering :: Chemistry
Topic :: Software Development :: Libraries :: Python Modules"""


def setup_cclib():

    # Import from setuptools only if requested.
    if 'egg' in sys.argv:
        sys.argv.pop(sys.argv.index('egg'))
        from setuptools import setup

    from distutils.core import setup

    # The list of packages to be installed.
    cclib_packages = [
        'cclib',
        'cclib.parser',
        'cclib.progress',
        'cclib.method',
        'cclib.bridge',
        'cclib.writer'
    ]

    doclines = __doc__.split("\n")

    setup(
        name = "cclib",
        version = "1.3.1",
        url = "http://cclib.github.io/",
        author = "cclib development team",
        author_email = "cclib-users@lists.sourceforge.net",
        maintainer = "cclib development team",
        maintainer_email = "cclib-users@lists.sourceforge.net",
        license = "LGPL",
        description = doclines[0],
        long_description = "\n".join(doclines[2:]),      
        classifiers = classifiers.split("\n"),
        platforms = ["Any."],
        packages = cclib_packages,
        package_dir = { 'cclib':'src/cclib' },
        scripts = ["src/scripts/ccget", "src/scripts/ccwrite", "src/scripts/cda"],
    )


if __name__ == '__main__':

    setup_cclib()
