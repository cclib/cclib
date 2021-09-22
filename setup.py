# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""cclib: parsers and algorithms for computational chemistry

cclib is a Python library that provides parsers for computational
chemistry log files. It also provides a platform to implement
algorithms in a package-independent manner.
"""

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
        version="1.7.1",
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
        packages=setuptools.find_packages(exclude=['*test*']),
        entry_points={
            'console_scripts': [
                'ccframe=cclib.scripts.ccframe:main',
                'ccget=cclib.scripts.ccget:ccget',
                'ccwrite=cclib.scripts.ccwrite:main',
                'cda=cclib.scripts.cda:main'
            ]
        },
        install_requires=[
            "packaging>=19.0",
            "numpy",
            "periodictable",
            "scipy>=1.2.0",
        ],

    )


if __name__ == '__main__':
    setup_cclib()
