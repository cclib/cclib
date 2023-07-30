# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Create a MANIFEST file for distributing soruces with distutils."""

import glob
import os


files = ['ANNOUNCE', 'CHANGELOG', 'LICENSE', 'README.md', 'THANKS']
files += ['setup.py']

source = 'cclib'
files.append(os.path.join(source, "__init__.py"))

folders = ['bridge', 'io', 'method', 'parser', 'progress']
for folder in folders:
    files.extend(glob.glob(os.path.join(source, folder, '*.py')))

for f in files:
    if not os.path.isfile(f):
        print(f"{f} does not exist")

print("\n".join(files), file=open("MANIFEST", "w"))
