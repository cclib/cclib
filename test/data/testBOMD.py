# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test Born-Oppenheimer molecular dynamics (BOMD) logfiles in cclib."""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from cclib.parser.data import ccData


class GenericBOMDTest:
    """Generic Born-Oppenheimer molecular dynamics unittest"""

    # To calculate the initial set of forces/velocities, programs
    # first converge the energy at the input geometry, so there is
    # always one more geometry/energy than MD step.
    nsteps = 35
    nenergies = 36

    def testdimscfenergies(self, data: "ccData") -> None:
        """Are the number of parsed energies consistent with the number of MD
        steps?
        """
        assert data.scfenergies.shape == (self.nenergies,)

    def testdimatomcoords(self, data: "ccData") -> None:
        """Are the number of parsed geometries consistent with the number of
        MD steps?
        """
        assert data.atomcoords.shape == (self.nenergies, 20, 3)

    def testdimtime(self, data: "ccData") -> None:
        """Are the number of time points consistent with the number of MD
        steps?
        """
        assert data.time.shape == (self.nsteps,)


class GaussianBOMDTest(GenericBOMDTest):
    """Customized Born-Oppenheimer molecular dynamics unittest"""

    # This may have something to do with our corrections for
    # extrapolation step rejection.
    nenergies = 35
