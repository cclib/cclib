#!/usr/bin/env bash

# run_pytest.sh: Run pytest on cclib with coverage checking. Requires `pytest`
# and `pytest-cov`.

set -eu

python -m pytest -v --capture=no --cov=cclib --cov-report=term --cov-report=html test &&
cd data && bash ./regression_download.sh && cd .. &&
python -m pytest -v --capture=no --cov=cclib --cov-report=term --cov-report=html --cov-append -k test_regression test/regression.py
