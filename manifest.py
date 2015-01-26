# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Create a MANIFEST file for distributing soruces with distutils."""

from __future__ import print_function

import glob
import os


files = ['ANNOUNCE', 'CHANGELOG', 'INSTALL', 'LICENSE', 'README','THANKS',]
files += ['setup.py']

source = os.path.join('src', 'cclib')
files.append(os.path.join(source, "__init__.py"))
files.append(os.path.join("src", "scripts", "ccget"))
files.append(os.path.join("src", "scripts", "ccwrite"))
files.append(os.path.join("src", "scripts", "cda"))

folders = ['bridge', 'method', 'parser', 'progress', 'writer']
for folder in folders:
    files.extend(glob.glob(os.path.join(source,folder,'*.py')))

for f in files:
    if not os.path.isfile(f):
        print("%s does not exist" % f)

print("\n".join(files), file=open("MANIFEST", "w"))
