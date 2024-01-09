#!/usr/bin/env bash

# run_pytest.bash: Run pytest on cclib with coverage checking. Requires
# `pytest` and `pytest-cov`.

set -exo pipefail

PYTEST_NWORKERS="${1}"

if [[ -n "${PYTEST_NWORKERS}" ]]; then
    PYTEST_PARALLELISM="-n ${PYTEST_NWORKERS} --dist worksteal"
else
    PYTEST_PARALLELISM=""
fi

PYTEST_ADDOPTS="-v -s ${PYTEST_PARALLELISM} --cov=cclib --cov-report=term --cov-report=xml:coverage-method.xml -m 'is_method'" python -m pytest
