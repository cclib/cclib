{ buildPythonPackage
, lib
, setuptools
, pytestCheckHook
, numpy
, scipy
, periodictable
, pyyaml
, ase
, openbabel-bindings
, h5py
, pyscf
, iodata
, psi4
, pandas
, biopython
, pyquante
}:

buildPythonPackage rec {
  pname = "cclib";
  version = "1.8.1";

  src = lib.cleanSource ../.;

  pyproject = true;

  build-system = [
    setuptools
  ];

  postPatch = ''
    # Upstream uses versioningit to set the version
    sed -i "/versioningit>=/d" pyproject.toml
    sed -i '/^name =.*/a version = "${version}"' pyproject.toml
    sed -i "/dynamic =/d" pyproject.toml
    echo '__version__ = "${version}"' > cclib/_version.py
  '';

  dependencies = [
    numpy
    scipy
    periodictable
    pyyaml
    ase
    openbabel-bindings
    h5py
    pyscf
    iodata
    psi4
    pandas
    biopython
    pyquante
  ];

  checkInputs = [ pytestCheckHook ];

  pythonImportsCheck = [ "cclib" ];

  disabledTests = [
    "test_ccread_url"
    "test_multi_url_io"
    "test_url_io"
    "test_url_seek"
  ];

  meta = with lib; {
    description = "Parsers and algorithms for computational chemistry logfiles";
    license = licenses.bsd3;
    maintainers = with maintainers; [ berquist sheepforce ];
  };
}
