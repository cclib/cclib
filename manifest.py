# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Create MANIFEST"""

__revision__ = "$Revision$"

import glob
import os


files = ['THANKS','README','INSTALL','ANNOUNCE','CHANGELOG','LICENSE',
         'setup.py']

source = os.path.join('src','cclib')
files.append(os.path.join(source,"__init__.py"))
files.append(os.path.join("src","scripts","ccget"))
files.append(os.path.join("src","scripts","cda"))

folders = ['method','bridge','parser','progress']
for folder in folders:
    files.extend(glob.glob(os.path.join(source,folder,'*.py')))

# Include data files
files.append(os.path.join("data", "regressionfiles.txt"))
files.append(os.path.join("data", "wget.sh"))
folders = glob.glob(os.path.join('data', '*'))
for folder in folders:
    basicfolders = glob.glob(os.path.join(folder, 'basic*'))
    for basicfolder in basicfolders:
        files.extend(glob.glob(os.path.join(basicfolder, "*")))

# Include test scripts
files.extend(glob.glob(os.path.join("test", "regression.py")))
files.extend(glob.glob(os.path.join("test", "test*.py")))
files.append(os.path.join("test", "__init__.py"))
for name in ['bettertest', 'methods']:
    files.append(os.path.join("test", "%s.py" % name))
files.append(os.path.join("test", "testdata"))

for f in files:
    if not os.path.isfile(f):
        print "%s does not exist" % f

print("\n".join(files), file=open("MANIFEST", "w"))
