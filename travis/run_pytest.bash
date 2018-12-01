#!/usr/bin/env bash

# run_pytest.bash: Run pytest on cclib with coverage checking. Requires
# `pytest` and `pytest-cov`.

set -euxo pipefail

python -m pytest -v --capture=no --cov=cclib --cov-report=term --cov-report=html test
pushd data
bash ./regression_download.sh
popd
python -m pytest -v --capture=no --cov=cclib --cov-report=term --cov-report=html --cov-append -k test_regression test/regression.py
