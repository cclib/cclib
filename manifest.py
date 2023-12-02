# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Create a MANIFEST file for distributing soruces with distutils."""

import glob
import os
from pathlib import Path

files = ["ANNOUNCE", "CHANGELOG", "LICENSE", "README.md", "THANKS"]
files += ["pyproject.toml", "setup.py"]

source = "cclib"
files.extend([os.path.join(source, fname) for fname in ["__init__.py", "py.typed"]])

folders = ["bridge", "io", "method", "parser", "progress"]
for folder in folders:
    files.extend(glob.glob(os.path.join(source, folder, "*.py")))

for f in files:
    if not os.path.isfile(f):
        print(f"{f} does not exist")

Path("MANIFEST").write_text("\n".join(files), encoding="utf-8")
