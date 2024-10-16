#
# Copyright (c) 2021-2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A reader for MolSSI quantum chemical JSON (QCSchema) files."""

from cclib.io.cjsonreader import CJSON as CJSONReader


class QCSchemaReader(CJSONReader):
    def __init__(self, source, *args, **kwargs):
        super().__init__(source, *args, **kwargs)
