# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the CM5 method in cclib"""

from cclib.method import CM5
from cclib.parser import QChem

import numpy as np

from ..test_data import getdatafile


class CM5Test:
    """Tests for Charge Model 5 (CM5) calculations."""

    def testcm5restricted(self) -> None:
        """Check that our computed CM5 charges match those parsed from a
        restricted calculation.
        """
        data, _ = getdatafile(QChem, "basicQChem5.4", ["dvb_sp.out"])
        res = CM5(data).charges()
        ref = data.atomcharges["cm5"]
        np.testing.assert_allclose(res, ref, atol=1.0e-6)

    def testcm5unrestricted(self) -> None:
        """Check that our computed CM5 charges match those parsed from an
        unrestricted calculation.
        """
        data, _ = getdatafile(QChem, "basicQChem5.4", ["dvb_sp_un.out"])
        res = CM5(data).charges()
        ref = data.atomcharges["cm5"]
        np.testing.assert_allclose(res, ref, atol=1.0e-6)
