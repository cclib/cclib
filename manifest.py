"""Create MANIFEST"""

import glob
import os

files = ['THANKS','README','INSTALL','ANNOUNCE','CHANGELOG','LICENSE',
         'setup.py']

source = os.path.join('src','cclib')
files.append(os.path.join(source,"__init__.py"))
files.append(os.path.join("src","scripts","ccget"))

folders = ['method','bridge','parser','progress']
for folder in folders:
    files.extend(glob.glob(os.path.join(source,folder,'*.py')))

for file in files:
    if not os.path.isfile(file):
        print "%s does not exist" % file

print >> open("MANIFEST","w"), "\n".join(files)
