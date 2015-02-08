from cclib.parser import ccData
from cclib.bridge import cclib2openbabel


def test_uracil():
    data = ccData(cclib2openbabel.readfile("../data/misc/uracil.xyz", "XYZ"))
    assert data.natom == 12


if __name__ == "__main__":
    test_uracil()