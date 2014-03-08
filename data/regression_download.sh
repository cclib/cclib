#!/bin/bash

if [ -d regression/.git ]; then
    if [ -e regression/README.md ] && [ `head -1 regression/README.md` == "cclib-data" ]; then
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

