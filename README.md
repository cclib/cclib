# cclib

**IMPORTANT for upcoming 2.0 release** We are preparing for the 2.0 release now that 1.8.1 is done.
Although most of the new features are on the unstable `main` branch, we will now be making some breaking changes to the default `master` branch.
See https://github.com/cclib/cclib/issues/1395 for more information.

- If you choose to follow `main`, we reserve the right to rewrite history until the final `v2.0` tag is created, after which `main` will replace `master` as the default branch.
- We do not expect to make any further tagged or versioned releases on the `master` branch.
- This message will disappear when the final release of 2.0, after any alphas/release candidates/etc. is made.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8280878.svg)](https://doi.org/10.5281/zenodo.8280878)
[![PyPI version](http://img.shields.io/pypi/v/cclib.svg?style=flat)](https://pypi.python.org/pypi/cclib)
[![GitHub release](https://img.shields.io/github/release/cclib/cclib.svg?style=flat)](https://github.com/cclib/cclib/releases)
[![build status](https://github.com/cclib/cclib/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/cclib/cclib/actions/workflows/ci.yml)
[![license](http://img.shields.io/badge/license-BSD-blue.svg?style=flat)](https://github.com/cclib/cclib/blob/master/LICENSE)

<img src="./logo.png" alt="cclib logo" width="100" />

cclib is a Python library that provides parsers for output files of computational chemistry packages. It also provides a platform for computational chemists to implement algorithms in a platform-independent way.

For more information, go to [https://cclib.github.io](https://cclib.github.io). There is a mailing list for questions at https://groups.google.com/g/cclib.

## Citing cclib

If you use cclib, please cite [the following article](https://doi.org/10.1063/5.0216778):

```bib
@article{10.1063/5.0216778,
    author = {Berquist, Eric and Dumi, Amanda and Upadhyay, Shiv and Abarbanel, Omri D. and Cho, Minsik and Gaur, Sagar and Cano Gil, Victor Hugo and Hutchison, Geoffrey R. and Lee, Oliver S. and Rosen, Andrew S. and Schamnad, Sanjeed and Schneider, Felipe S. S. and Steinmann, Casper and Stolyarchuk, Maxim and Vandezande, Jonathon E. and Zak, Weronika and Langner, Karol M.},
    title = {cclib 2.0: An updated architecture for interoperable computational chemistry},
    journal = {The Journal of Chemical Physics},
    volume = {161},
    number = {4},
    pages = {042501},
    year = {2024},
    month = {07},
    issn = {0021-9606},
    doi = {10.1063/5.0216778},
    url = {https://doi.org/10.1063/5.0216778},
}
```

For the original paper on cclib (before version 2.0), please refer to [the following article](https://doi.org/10.1002/jcc.20823):

```bib
@article{https://doi.org/10.1002/jcc.20823,
author = {O'boyle, Noel M. and Tenderholt, Adam L. and Langner, Karol M.},
title = {cclib: A library for package-independent computational chemistry algorithms},
journal = {Journal of Computational Chemistry},
volume = {29},
number = {5},
pages = {839-845},
keywords = {computational chemistry, algorithms, Python},
doi = {https://doi.org/10.1002/jcc.20823},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.20823},
year = {2008}
}
```
