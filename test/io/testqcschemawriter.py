# -*- coding: utf-8 -*-
#
# Copyright (c) 2021, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for QCSchema writer."""

import os
import unittest
from pathlib import Path

import qcschema

import cclib

__filedir__ = os.path.dirname(__file__)
__filepath__ = Path(os.path.realpath(__filedir__))
__datadir__ = __filepath__.joinpath("..", "..").resolve()


class QCSchemaWriterTest(unittest.TestCase):
    def test_validate_output_hf(self):
        fpath = __datadir__ / "data" / "QChem" / "basicQChem5.1" / "dvb_sp.out"
        data = cclib.io.ccread(str(fpath))
        writer = cclib.io.qcschemawriter.QCSchemaWriter(data)
        qcschema.validate(writer.as_dict(), schema_type="output")
