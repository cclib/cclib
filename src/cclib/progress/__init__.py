# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import sys

if 'PyQt4' in list(sys.modules.keys()):
    from cclib.progress.qt4progress import Qt4Progress

from cclib.progress.textprogress import TextProgress

