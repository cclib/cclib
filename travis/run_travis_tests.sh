#!/bin/sh

# run_travis_tests.sh: Run the tests for Travis CI.

cd test
python test_bridge.py
python test_utils.py
python test_method.py
python test_parser.py
python test_io.py
python test_data.py --status --terse
cd ../data && bash regression_download.sh
cd ../test && python run_regressions.py --status --traceback
