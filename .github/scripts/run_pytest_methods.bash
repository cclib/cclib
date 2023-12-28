#!/usr/bin/env bash

# run_pytest.bash: Run pytest on cclib with coverage checking. Requires
# `pytest` and `pytest-cov`.

set -euxo pipefail

PYTEST_ADDOPTS='-v -s --cov=cclib --cov-report=term --cov-report=xml:coverage-method.xml -m "is_method"' python -m pytest
