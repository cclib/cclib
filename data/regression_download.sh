#!/bin/bash

if [ -d regression/.git ]; then
    echo "Updating regression files..."
    cd regression
    git pull
    cd ..
else
    echo "Downloading repository of regression files..."
    git clone git@github.com:langner/cclib-data.git regression
fi

