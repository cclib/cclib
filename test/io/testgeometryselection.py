# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Geometry selection in single-geometry exports."""

import json
import xml.etree.ElementTree as ET

from cclib.io import cjsonwriter, cmlwriter, filewriter
from cclib.io.xyzwriter import XYZ
from cclib.parser.data import ccData

import numpy as np
import pytest


@pytest.fixture
def geometries():
    # The first H2 geometry is bonded; the second has dissociated.
    return ccData(
        {
            "natom": 2,
            "atomnos": [1, 1],
            "charge": 0,
            "mult": 1,
            "atomcoords": [[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]], [[0.0, 0.0, 0.0], [9.0, 0.0, 0.0]]],
            "metadata": {},
        }
    )


@pytest.fixture(params=[False, True], ids=["without-openbabel", "with-openbabel"])
def babel(request, monkeypatch):
    if request.param and not filewriter._has_openbabel:
        pytest.skip("Open Babel not installed")
    for module in (filewriter, cmlwriter, cjsonwriter):
        monkeypatch.setattr(module, "_has_openbabel", request.param)
    return request.param


def exported_coordinates(writer):
    output = writer.generate_repr()
    if isinstance(writer, cmlwriter.CML):
        atoms = ET.fromstring(output).findall(".//{*}atom")
        return np.array([[float(a.attrib[k]) for k in ("x3", "y3", "z3")] for a in atoms])
    return np.array(json.loads(output)["atoms"]["coords"]["3d"]).reshape(-1, 3)


@pytest.mark.parametrize("writer_class", [cmlwriter.CML, cjsonwriter.CJSON])
@pytest.mark.parametrize(
    "indices, expected",
    [
        (None, 1),
        ([], 1),
        (0, 0),
        ([0], 0),
        ({0}, 0),
        ([0, 0], 0),
        (1, 1),
        ([-1], 1),
        (-2, 0),
        ([0, -2], 0),
    ],
)
def test_selected_geometry(writer_class, indices, expected, geometries, babel):
    original = geometries.atomcoords.copy()
    writer = writer_class(geometries, indices=indices)
    np.testing.assert_allclose(exported_coordinates(writer), original[expected])
    assert writer.ccdata is geometries
    np.testing.assert_array_equal(geometries.atomcoords, original)
    if babel:
        assert writer.obmol.NumBonds() == (1 if expected == 0 else 0)
        assert writer.obmol.GetAtom(2).GetX() == pytest.approx(original[expected, 1, 0])
        if writer_class is cmlwriter.CML:
            bonds = ET.fromstring(writer.generate_repr()).findall(".//{*}bond")
            assert len(bonds) == (1 if expected == 0 else 0)
        else:
            bonds = writer.as_dict()["bonds"]["connections"]["index"]
            assert bonds == ([0, 1] if expected == 0 else [])


@pytest.mark.parametrize("writer_class", [cmlwriter.CML, cjsonwriter.CJSON])
@pytest.mark.parametrize("indices", [2, -3, [99]])
def test_invalid_geometry(writer_class, indices, geometries, babel):
    with pytest.raises(IndexError, match="Geometry index out of range"):
        writer_class(geometries, indices=indices)


@pytest.mark.parametrize("writer_class", [cmlwriter.CML, cjsonwriter.CJSON])
def test_multiple_geometries_rejected(writer_class, geometries, babel):
    with pytest.raises(ValueError, match="supports only one geometry"):
        writer_class(geometries, indices=[0, 1])


@pytest.mark.parametrize("writer_class", [cmlwriter.CML, cjsonwriter.CJSON])
def test_missing_coordinates(writer_class, geometries, babel):
    del geometries.atomcoords
    # Existing atom-only exports should still work.
    writer_class(geometries).generate_repr()
    with pytest.raises(IndexError, match="Geometry index out of range"):
        writer_class(geometries, indices=0)


def test_xyz_still_accepts_multiple_geometries(geometries, babel):
    output = XYZ(geometries, indices=[0, 1]).generate_repr()
    assert "0.7400000000" in output
    assert "9.0000000000" in output


@pytest.mark.parametrize("method", ["writecml", "writejson"])
def test_public_write_api(method, geometries, babel):
    output = getattr(geometries, method)(indices=0)
    if method == "writecml":
        atoms = ET.fromstring(output).findall(".//{*}atom")
        assert float(atoms[1].attrib["x3"]) == pytest.approx(0.74)
    else:
        assert json.loads(output)["atoms"]["coords"]["3d"][3] == pytest.approx(0.74)
