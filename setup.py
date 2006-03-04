from distutils.core import setup
setup(name="cclib",
      version="0.1.2",
      author="Noel O'Boyle",
      author_email="baoilleach@users.sf.net",
      url="http://gausssum.sf.net",
      scripts=['src/gausssum/gausssum.py'],
      package_dir = {'cclib':'src/cclib'},
      packages=['cclib'])
