#!/usr/bin/env bash

# run_pytest.bash: Run pytest on cclib with coverage checking. Requires
# `pytest` and `pytest-cov`.

set -euxo pipefail

PYTEST_ADDOPTS='-v -s --cov=cclib --cov-report=term --cov-report=xml:coverage-unit.xml -m "not is_method"' python -m pytest
pushd data
bash ./regression_download.sh
popd
PYTEST_ADDOPTS='-v -s --cov=cclib --cov-report=term --cov-report=xml:coverage-regression.xml --cov-append' python -m pytest test/regression.py
PYTEST_ADDOPTS='-v -s --cov=cclib --cov-report=term --cov-report=xml:coverage-regression.xml --cov-append ' python -m pytest test/regression_io.py
