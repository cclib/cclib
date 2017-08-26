#!/bin/bash
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

if [ -d regression/.git ]; then
    if [ -e regression/README.md ] && [ $(head -1 regression/README.md | cut -d " " -f 2) == "cclib-data" ]; then
        echo "Updating regression files..."
        cd regression
        if [[ -n $(git status -s) ]]; then
            echo "WARNING: stashing uncommited changes."
            git stash save --keep-index
        fi
        git pull
        cd ..
    else
        echo "The repository in regression/ does not seem to contain cclib-data."
        echo "Please remove it and run this script again."
    fi
else
    echo "Downloading repository of regression files..."
    git clone https://github.com/cclib/cclib-data.git regression
fi
