from distutils.core import setup
setup(name="cclib",
      version="0.1",
      author="Noel O'Boyle",
      author_email="baoilleach@users.sf.net",
      url="http://cclib.sourceforge.net",
      package_dir = {'cclib':'src/cclib'},
      packages=['cclib','cclib.parser','cclib.progress','cclib.method'])
