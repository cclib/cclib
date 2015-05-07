import os
import unittest

from cclib.parser import ccData
from cclib.bridge import cclib2openbabel


class OpenbabelTest(unittest.TestCase):
    """Tests for the cclib2openbabel bridge in cclib."""

    def setUp(self):
        self.path = os.path.abspath(os.path.dirname(__file__))

    def test_xyz_uracyl(self):
        """Try to load an XYZ file with uracyl through Openbabel"""
        data = cclib2openbabel.readfile(self.path + "/uracil.xyz", "XYZ")
        assert data.natom == 12


if __name__ == "__main__":

    unittest.main()
