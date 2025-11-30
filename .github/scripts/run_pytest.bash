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

PYTEST_ADDOPTS="-v -s ${PYTEST_PARALLELISM} ${PYTEST_COVERAGE} --cov-report=xml:coverage-unit.xml -m 'not is_method' -k 'not (data and pyscf)'" python -m pytest
pushd data
bash ./regression_download.sh
popd
#PYTEST_ADDOPTS="-v -s ${PYTEST_PARALLELISM} ${PYTEST_COVERAGE} --cov-report=xml:coverage-regression.xml --cov-append" python -m pytest test/regression.py
#PYTEST_ADDOPTS="-v -s ${PYTEST_PARALLELISM} ${PYTEST_COVERAGE} --cov-report=xml:coverage-regression.xml --cov-append" python -m pytest test/regression_io.py
