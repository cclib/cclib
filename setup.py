"""cclib -- parsers and algorithms for computational chemistry

cclib is a Python library that provides parsers for computational
chemistry log files. It also provides a platform to implement
algorithms in a package-independent manner.
"""

from distutils.core import setup

# Chosen from http://www.python.org/pypi?:action=list_classifiers
classifiers = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)
Programming Language :: Python
Topic :: Software Development :: Libraries :: Python Modules
Operating System :: OS Independent
Topic :: Scientific/Engineering :: Chemistry
"""

doclines = __doc__.split("\n")

setup(name="cclib",
      version="0.6",
      author="cclib development team",
      author_email="cclib-users@lists.sourceforge.net",
      url="http://cclib.sf.net",
      description=doclines[0],
      long_description = "\n".join(doclines[2:]),      
      classifiers=filter(None, classifiers.split("\n")),
      license="LGPL",
      platforms=["Any."],
      scripts=["src/scripts/ccget"],
      package_dir = {'cclib':'src/cclib'},
      packages=['cclib','cclib.parser','cclib.progress','cclib.method','cclib.bridge'])
