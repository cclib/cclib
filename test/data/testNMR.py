# -*- coding: utf-8 -*-
#
# Copyright (c) 2021, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test NMR logfiles in cclib."""

import os
import numpy
import unittest
import pytest

from skip import skipForParser

__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericNMRTest(unittest.TestCase):
    """Generic NMR unittest"""
    
    def testkeys(self):
        """Check dictionary keys are ints."""
        key_types = set([type(key) for key in self.data.nmrtensors.keys()])
        assert len(key_types) == 1 and int in key_types

    def testsize(self):
        """Check to make sure there are the correct number of tensors parsed"""
        assert len(self.data.nmrtensors) == self.data.natom
        assert len(self.data.nmrtensors[0]) == 4
        assert self.data.nmrtensors[0]["total"].shape == (3, 3)
        
    def testisotropic(self):
        """Check the total isotropic value matches the computed value."""
        tensor = self.data.nmrtensors[0]
        total = 0.0
        for t_type in ("diamagnetic", "paramagnetic"):
            total += sum(numpy.linalg.eigvals(tensor[t_type])) / len(numpy.linalg.eigvals(tensor[t_type]))
        
        assert total == pytest.approx(tensor['isotropic'], abs = 3)


class GenericNMRCouplingTest(unittest.TestCase):
    """Generic NMR spin-spin coupling unittest"""
    
    def testkeys(self):
        """Check dictionary keys are ints."""
        assert all(set([all((type(key[0]) == int, type(key[1]) == int)) for key in self.data.nmrcouplingtensors.keys()]))
        
        assert all(
            [
                all((type(isotopekey[0]) == int, type(isotopekey[1]) == int))
                for isotopes in self.data.nmrcouplingtensors.values()
                for isotopekey in isotopes.keys()
            ]
        )

    def testsize(self):
        """Check to make sure there are the correct number of tensors parsed"""
        assert len(self.data.nmrcouplingtensors) == 139
        tensor = list(list(self.data.nmrcouplingtensors.values())[0].values())[0]
        assert len(tensor) == 7
        assert tensor['total'].shape == (3, 3)
    
    def testisotropic(self):
        """Check the total isotropic value matches the computed value."""
        tensor = list(list(self.data.nmrcouplingtensors.values())[0].values())[0]
        # Check the total isotropic value matches the computed value.
        total = 0.0
        for t_type in ("diamagnetic", "paramagnetic", "fermi", "spin-dipolar", "spin-dipolar-fermi"):
            total += sum(numpy.linalg.eigvals(tensor[t_type])) / len(numpy.linalg.eigvals(tensor[t_type]))
        
        assert total == pytest.approx(tensor['isotropic'], abs = 3)

if __name__ == "__main__":
    import sys

    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite

    suite = DataSuite(['NMR'])
    suite.testall()
