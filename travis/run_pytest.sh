#!/bin/sh

# run_pytest.sh: Run pytest on cclib with coverage checking. Requires `pytest`
# and `pytest-cov`.

pytest -v --doctest-modules --cov=cclib test
