#!/usr/bin/env bash

# run_pytest.bash: Run pytest on cclib with coverage checking. Requires
# `pytest` and `pytest-cov`.

set -euxo pipefail

python -m pytest -v --capture=no --cov=cclib --cov-report=term --cov-report=xml:coverage-method.xml --terse test -k "test_method"
