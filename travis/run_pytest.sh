#!/bin/sh

# run_pytest.sh: Run pytest on cclib with coverage checking. Requires `pytest`
# and `pytest-cov`.

pytest -v --doctest-modules --capture=no --cov=cclib test
pytest -v --capture=no --cov=cclib --cov-append -k test_regression test/regression.py
