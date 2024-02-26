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
  version = "1.8";

  src = lib.cleanSource ../.;

  pyproject = true;

  nativeBuildInputs = [
    setuptools
  ];

  propagatedBuildInputs = [
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

  disabledTests = [
    "test_url_io"
    "test_multi_url_io"
    "test_url_seek"
  ];

  meta = with lib; {
    description = "Parsers and algorithms for computational chemistry logfiles ";
    license = licenses.bsd3;
    maintainers = with maintainers; [ sheepforc3 ];
  };
}
