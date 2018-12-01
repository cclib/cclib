#!/usr/bin/env bash

# run_unittest.bash: Run cclib tests using the built-in `unittest` module.

python -m test.test_bridge &&
python -m test.test_io &&
python -m test.test_method &&
python -m test.test_parser &&
python -m test.test_utils &&
python -m test.test_data --terse &&
cd data && bash regression_download.sh &&
cd .. && python -m test.regression --traceback
