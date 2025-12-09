# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test NMR logfiles in cclib."""

import numpy
import pytest
from skip import skipForParser


class GenericNMRTest:
    """Generic NMR unittest"""
    # How many different tensor types do we expect?
    num_tensors = 4

    def testkeys(self, data) -> None:
        """Check dictionary keys are ints."""
        key_types = {type(key) for key in data.nmrtensors.keys()}
        assert len(key_types) == 1 and int in key_types

    def testsize(self, data) -> None:
        """Check to make sure there are the correct number of tensors parsed"""
        assert len(data.nmrtensors) == data.natom
        assert len(data.nmrtensors[0]) == self.num_tensors
        assert data.nmrtensors[0]["total"].shape == (3, 3)

    def testisotropic(self, data) -> None:
        """Check the total isotropic value matches the computed value."""
        tensor = data.nmrtensors[0]
        total = 0.0
        other_tensors = {
            key: value for key, value in tensor.items() if key not in ["total", "isotropic"]
        }
        if len(other_tensors) > 0:
            for t_type in other_tensors:
                eigvals = numpy.linalg.eigvals(tensor[t_type])
                total += numpy.mean(eigvals)

            assert total == pytest.approx(tensor["isotropic"], abs=3)

        eigvals = numpy.linalg.eigvals(tensor["total"])
        total = numpy.mean(eigvals)
        assert total == pytest.approx(tensor["isotropic"], abs=3)

    def testtotalshift(self, data):
        """Test the total chemical shift matches our expected value."""
        # TODO: A bit crude, but at least it will alert to unexpected changes in structure.
        assert sum(
            [data.nmrtensors[atom]["isotropic"] for atom in data.nmrtensors]
        ) == pytest.approx(1455, abs=5)


class GaussianNMRTest(GenericNMRTest):
    num_tensors = 2


class PySCFNMRTest(GenericNMRTest):
    num_tensors = 2


class TurbomoleNMRTest(GenericNMRTest):
    num_tensors = 5


class GenericNMRCouplingTest:
    """Generic NMR spin-spin coupling unittest"""
    # The total molecular mass, varies depending on how the program handles isotope masses.
    tot_mass = 140.111
    
    # How many different tensor types do we expect?
    num_tensors = 7

    # How many atom pairs do we expect?
    num_couplings = 190

    @skipForParser("Turbomole", "atommasses is not available")
    def testmass(self, data) -> None:
        """Check we are using 13C"""
        assert sum(data.atommasses) == pytest.approx(self.tot_mass, 1e-2)

    def testisotopes(self, data) -> None:
        """Check we are using 13C"""
        for atom1, atom2 in data.nmrcouplingtensors:
            for isotope1, isotope2 in data.nmrcouplingtensors[(atom1, atom2)]:
                if data.atomnos[atom1] == 6:
                    assert isotope1 == 13

                if data.atomnos[atom2] == 6:
                    assert isotope2 == 13

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
        assert len(data.nmrcouplingtensors) == self.num_couplings
        tensor = list(list(data.nmrcouplingtensors.values())[0].values())[0]
        assert len(tensor) == self.num_tensors
    
    @skipForParser("Gaussian", "no coupling tensors are available")
    def testtotal(self, data) -> None:
        """Does the total tensor have the correct shape?"""
        tensor = list(list(data.nmrcouplingtensors.values())[0].values())[0]
        assert tensor["total"].shape == (3, 3)

    @skipForParser("Gaussian", "no coupling tensors are available")
    @skipForParser("PySCF", "only the total tensor is available")
    @skipForParser("Turbomole", "only the total tensor is available")
    def testtensors(self, data) -> None:
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

    @skipForParser("Gaussian", "no coupling tensors are available")
    def testtotaltensor(self, data) -> None:
        """Check the total isotropic value matches the computed value."""
        tensor = list(list(data.nmrcouplingtensors.values())[0].values())[0]
        # Check the total isotropic value matches the computed value.
        eigvals = numpy.linalg.eigvals(tensor["total"])
        total = numpy.mean(eigvals)

        assert total == pytest.approx(tensor["isotropic"], abs=3)

    def testtotalcoupling(self, data):
        """Test the total chemical shift matches our expected value."""
        assert sum(
            [
                list(data.nmrcouplingtensors[atom].values())[0]["isotropic"]
                for atom in data.nmrcouplingtensors
            ]
        ) == pytest.approx(1825, abs=5)


class OrcaNMRCouplingTest(GenericNMRCouplingTest):
    """Orca NMR spin-spin coupling unittest"""
    # In Orca, the NMR isotopes do not impact the main calculation isotopes.
    tot_mass = 130.190
    num_tensors = 7


class GaussianNMRCouplingTest(GenericNMRCouplingTest):
    """Gaussian NMR spin-spin coupling unittest"""
    num_tensors = 1


class PySCFNMRCouplingTest(GenericNMRCouplingTest):
    """PySCF NMR spin-spin coupling unittest"""
    num_tensors = 2


class TurbomoleNMRCouplingTest(GenericNMRCouplingTest):
    """Turbomole NMR spin-spin coupling unittest"""
    num_tensors = 2
