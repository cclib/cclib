# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

from pathlib import Path

from cclib.method import volume

import numpy as np
import pytest

__filedir__ = Path(__file__).parent
__filepath__ = Path(__filedir__).resolve()
__regdir__ = (__filepath__ / ".." / "data" / "regression").resolve()


class VolumeRegressionTests:
    """Regression tests for the cclib.method.volume.Volume class.

    The h5cube spec is located at
    https://h5cube-spec.readthedocs.io/en/latest/cubeformat.html.
    """

    def test_orca_cube_ammonia(self):
        """Test that ORCA-generated h5cube-formatted cube files parse.

        Previously, the line for DSET_IDS was being read as data.

        The example cube file is from
        https://github.com/cclib/cclib/pull/1688#issuecomment-3431984904.
        """
        fpath = __regdir__ / "method" / "volume" / "NH3.mo4a.cube"
        vol = volume.read_from_cube(fpath)
        assert vol.data.shape == (5, 5, 5)
        assert vol.data[0, 0, 0] == pytest.approx(-2.97850e-10)
        assert vol.data[4, 4, 4] == pytest.approx(-6.46340e-11)

    def test_orca_cube_hexazine(self):
        """Test that ORCA-generated h5cube-formatted cube files parse.

        Previously, the line for DSET_IDS was being read as data.

        The example cube file is from
        https://discuss.avogadro.cc/t/orca-generated-cube-file-not-displaying-orbital/6346.
        """
        fpath = __regdir__ / "method" / "volume" / "hexazine.mo20a.cube"
        vol = volume.read_from_cube(fpath)
        assert vol.data.shape == (40, 40, 40)
        np.testing.assert_allclose(
            vol.data[0, 0, :6],
            np.asarray(
                [0.00000e00, 0.00000e00, 0.00000e00, -2.69165e-11, -1.44655e-10, -2.49139e-10]
            ),
        )
