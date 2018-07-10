# Documentation for cclib

This directory contains the source of the current official website and documentation for [cclib](https://github.com/cclib/cclib), available on GitHub pages at http://cclib.github.io.

## How to update documentation

The website is generated using [Sphinx](http://sphinx-doc.org/) with some [custom adjustments](https://github.com/cclib/sphinx_rtd_theme/tree/cclib). The [reStructuredText](http://sphinx-doc.org/rest.html) sources are in the `sphinx` subdirectory, and executing `make` places the built website in `sphinx/_build/html`. Some of the content is generated automatically from the cclib code using autodoc or custom Python scripts, and this should be handled by the Makefile.
