# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test NMR logfiles in cclib."""

import numpy
import pytest


class GenericNMRTest:
    """Generic NMR unittest"""

    def testkeys(self, data) -> None:
        """Check dictionary keys are ints."""
        key_types = {type(key) for key in data.nmrtensors.keys()}
        assert len(key_types) == 1 and int in key_types

    def testsize(self, data) -> None:
        """Check to make sure there are the correct number of tensors parsed"""
        assert len(data.nmrtensors) == data.natom
        assert len(data.nmrtensors[0]) == 4
        assert data.nmrtensors[0]["total"].shape == (3, 3)

    def testisotropic(self, data) -> None:
        """Check the total isotropic value matches the computed value."""
        tensor = data.nmrtensors[0]
        total = 0.0
        for t_type in ("diamagnetic", "paramagnetic"):
            eigvals = numpy.linalg.eigvals(tensor[t_type])
            total += numpy.mean(eigvals)

        assert total == pytest.approx(tensor["isotropic"], abs=3)


class GenericNMRCouplingTest:
    """Generic NMR spin-spin coupling unittest"""

    def testkeys(self, data) -> None:
        """Check dictionary keys are ints."""
        assert all(
            {
                all((isinstance(key[0], int), isinstance(key[1], int)))
                for key in data.nmrcouplingtensors.keys()
            }
        )

        assert all(
            [
                all((isinstance(isotopekey[0], int), isinstance(isotopekey[1], int)))
                for isotopes in data.nmrcouplingtensors.values()
                for isotopekey in isotopes.keys()
            ]
        )

    def testsize(self, data) -> None:
        """Check to make sure there are the correct number of tensors parsed"""
        assert len(data.nmrcouplingtensors) == 139
        tensor = list(list(data.nmrcouplingtensors.values())[0].values())[0]
        assert len(tensor) == 7
        assert tensor["total"].shape == (3, 3)

    def testisotropic(self, data) -> None:
        """Check the total isotropic value matches the computed value."""
        tensor = list(list(data.nmrcouplingtensors.values())[0].values())[0]
        # Check the total isotropic value matches the computed value.
        total = 0.0
        for t_type in (
            "diamagnetic",
            "paramagnetic",
            "fermi",
            "spin-dipolar",
            "spin-dipolar-fermi",
        ):
            eigvals = numpy.linalg.eigvals(tensor[t_type])
            total += numpy.mean(eigvals)

        assert total == pytest.approx(tensor["isotropic"], abs=3)
