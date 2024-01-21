# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

from pathlib import Path

import cclib

__filedir__ = Path(__file__).parent
__filepath__ = Path(__filedir__).resolve()
__regdir__ = (__filepath__ / ".." / "data" / "regression").resolve()


class XYZRegressionTests:
    def test_xyz_not_turbomole(self):
        """Ensure XYZ file isn't misrecognized as a Turbomole file.

        From https://github.com/cclib/cclib/issues/1207.
        """
        fpath = __regdir__ / "io" / "xyz" / "1207.xyz"
        # ccopen assumes something is a parsable logfile and doesn't try any
        # fallback mechanisms.
        logfile = cclib.io.ccopen(str(fpath))
        assert logfile is None
        # ccread tries alternative file reading mechanisms.
        data = cclib.io.ccread(str(fpath))
        assert set(data.getattributes().keys()) == {"atomcoords", "atommasses", "atomnos", "natom"}
