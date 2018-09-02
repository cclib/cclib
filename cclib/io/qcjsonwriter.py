# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A writer for MolSSI quantum chemical JSON (QCJSON) files.
"""

from .cjsonwriter import CJSON as CJSONWriter


class QCJSONWriter(CJSONWriter):
    def __init__(self, source, *args, **kwargs):
        super(QCJSONWriter, self).__init__(source, *args, **kwargs)
