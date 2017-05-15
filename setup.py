# -*- coding: utf-8 -*-
#
# Copyright (c) 2016, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

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
License :: OSI Approved :: BSD License
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
        'cclib.bridge',
        'cclib.io',
        'cclib.method',
        'cclib.parser',
        'cclib.progress',
    ]

    doclines = __doc__.split("\n")

    setup(
        name = "cclib",
        version = "1.5.1",
        url = "http://cclib.github.io/",
        author = "cclib development team",
        author_email = "cclib-users@lists.sourceforge.net",
        maintainer = "cclib development team",
        maintainer_email = "cclib-users@lists.sourceforge.net",
        license = "BSD 3-Clause License",
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
