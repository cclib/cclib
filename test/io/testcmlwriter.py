# Copyright (c) 2025-2026, the cclib development team

"""Unit tests for the CML writer."""

import numpy as np

import cclib
from cclib.parser.data import ccData


class CMLWriterTest:
    def test_geometry_indices(self) -> None:
        """The selected geometry is used for CML coordinates."""
        data = ccData({
            "natom": 1,
            "atomnos": np.array([1]),
            "atomcoords": np.array([[[0.0, 0.0, 0.0]], [[9.0, 0.0, 0.0]]]),
            "metadata": {},
        })
        first = cclib.io.cmlwriter.CML(data, indices=[0]).generate_repr()
        last = cclib.io.cmlwriter.CML(data, indices=[1]).generate_repr()
        assert 'x3="0.0000000000"' in first
        assert 'x3="9.0000000000"' in last
