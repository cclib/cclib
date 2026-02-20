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

if [[ "${PYTEST_ENABLE_COVERAGE}" == "true" ]]; then
    PYTEST_COVERAGE="--cov=cclib --cov-report=term"
else
    PYTEST_COVERAGE=""
fi

PYTEST_ADDOPTS="-v -s ${PYTEST_PARALLELISM} ${PYTEST_COVERAGE} --cov-report=xml:coverage-method.xml -m 'is_method'" python -m pytest
PYTEST_ADDOPTS="-v -s ${PYTEST_PARALLELISM} ${PYTEST_COVERAGE} --cov-report=xml:coverage-method.xml --cov-append -m 'is_method'" python -m pytest test/regression_method.py
