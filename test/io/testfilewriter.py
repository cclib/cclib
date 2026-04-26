# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for general file writer."""

import os

import cclib

import pytest

__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")


class FileWriterTest:
    def test_init(self):
        """Does the class initialize properly?"""

        # You cannot instantiate a class with abstract methods.
        fpath = os.path.join(__datadir__, "data/ADF/basicADF2007.01/dvb_gopt.adfout")
        data = cclib.io.ccread(fpath)
        with pytest.raises(TypeError):
            cclib.io.filewriter.Writer(data)

    def test_fix_indices(self) -> None:
        """Check the implementation of cleaning up geometry indices."""
        fix_indices = cclib.io.filewriter._fix_indices

        assert fix_indices(None, None) == set()
        assert fix_indices(None, 10) == set()
        assert fix_indices(4, None) == {4}
        assert fix_indices(4, 10) == {4}
        assert fix_indices([4], None) == {4}
        assert fix_indices([4, 5], None) == {4, 5}
        assert fix_indices([4, 4], None) == {4}
        assert fix_indices([-2, -1], 10) == {8, 9}
        assert fix_indices([3, 4], 10) == {3, 4}
        with pytest.raises(RuntimeError):
            fix_indices([-1], None)
