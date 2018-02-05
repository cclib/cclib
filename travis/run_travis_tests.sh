#!/bin/sh

# run_travis_tests.sh: Run the tests for Travis CI.

python -m test.test_bridge &&
python -m test.test_utils &&
python -m test.test_method &&
python -m test.test_parser &&
python -m test.test_io &&
python -m test.test_data --status --terse &&
cd data && bash regression_download.sh &&
cd .. && python -m test.regression --status --traceback
