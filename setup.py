# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""cclib: parsers and algorithms for computational chemistry

cclib is a Python library that provides parsers for computational
chemistry log files. It also provides a platform to implement
algorithms in a package-independent manner.
"""

from __future__ import absolute_import, with_statement

import setuptools


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

    doclines = __doc__.split("\n")

    setuptools.setup(
        name="cclib",
        version="1.5.3",
        url="http://cclib.github.io/",
        author="cclib development team",
        author_email="cclib-users@lists.sourceforge.net",
        maintainer="cclib development team",
        maintainer_email="cclib-users@lists.sourceforge.net",
        license="BSD 3-Clause License",
        description=doclines[0],
        long_description="\n".join(doclines[2:]),
        classifiers=classifiers.split("\n"),
        platforms=["Any."],
        packages=setuptools.find_packages('src'),
        package_dir={'': 'src'},
        entry_points={
            'console_scripts': [
                'ccget=scripts.ccget:ccget',
                'ccwrite=scripts.ccwrite:main',
                'cda=scripts.cda:main'
            ]
        }

    )


if __name__ == '__main__':
    setup_cclib()
