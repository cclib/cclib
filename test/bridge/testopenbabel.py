from cclib.parser import ccData
from cclib.bridge import cclib2openbabel


def test_uracil():
    data = cclib2openbabel.readfile("uracil.xyz", "XYZ")
    assert data.natom == 12


if __name__ == "__main__":
    test_uracil()
