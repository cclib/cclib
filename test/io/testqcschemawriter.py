#
# Copyright (c) 2021-2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for the QCSchema writer."""

import os
import unittest
from pathlib import Path

import cclib

import qcschema

__filedir__ = os.path.dirname(__file__)
__filepath__ = Path(os.path.realpath(__filedir__))
__datadir__ = __filepath__.joinpath("..", "..").resolve()


def validate_output(outputpath):
    data = cclib.io.ccread(str(outputpath))
    writer = cclib.io.qcschemawriter.QCSchemaWriter(data)
    qcschema.validate(writer.as_dict(), schema_type="output")


class QCSchemaWriterTest(unittest.TestCase):
    def test_validate_output_b3lyp_energy(self):
        fpath = __datadir__ / "data" / "QChem" / "basicQChem5.1" / "dvb_sp.out"
        validate_output(fpath)

    def test_validate_output_b3lyp_gradient(self):
        # This is not quite right, because a geometry optimization is a
        # procedure (repeated use of the gradient driver) and not a driver
        # itself. In the future, when requesting QCSChema output, a geometry
        # optimization should split into multiple force outputs (might depend
        # on #657).
        fpath = __datadir__ / "data" / "QChem" / "basicQChem5.1" / "dvb_gopt.out"
        validate_output(fpath)

    def test_validate_output_b3lyp_hessian(self):
        fpath = __datadir__ / "data" / "QChem" / "basicQChem5.1" / "dvb_ir.out"
        validate_output(fpath)
