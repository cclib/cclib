#
# Copyright (c) 2018-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A writer for MolSSI quantum chemical JSON (QCJSON) files."""

from .cjsonwriter import CJSON as CJSONWriter


class QCJSONWriter(CJSONWriter):
    def __init__(self, source, *args, **kwargs):
        super().__init__(source, *args, **kwargs)
