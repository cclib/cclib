#!/bin/sh

# install_pyquante.sh: Install the original PyQuante. Only for Python
# 2.7. Used for Travis CI.

curl -L https://sourceforge.net/projects/pyquante/files/PyQuante-1.6/PyQuante-1.6.5/PyQuante-1.6.5.tar.gz > PyQuante-1.6.5.tar.gz
tar xf PyQuante-1.6.5.tar.gz
cd PyQuante-1.6.5
python setup.py install
cd ..
