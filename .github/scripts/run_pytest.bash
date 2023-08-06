#!/usr/bin/env bash

# run_pytest.bash: Run pytest on cclib with coverage checking. Requires
# `pytest` and `pytest-cov`.

set -euxo pipefail

python -m pytest -v --capture=no --cov=cclib --cov-report=term --cov-report=xml:coverage-unit.xml --terse test -k "not test_method"
pushd data
bash ./regression_download.sh
popd
python -m pytest -v --capture=no --cov=cclib --cov-report=term --cov-report=xml:coverage-regression.xml --cov-append -k test_regression test/regression.py
python -m pytest -v -s --cov=cclib --cov-report=term --cov-report=xml:coverage-regression.xml --cov-append test/regression_io.py
