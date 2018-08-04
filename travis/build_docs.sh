#!/usr/bin/env bash

# build_docs.sh: Build the documentation as part of Travis CI.

set -o errexit

bold=$(tput bold)
normal=$(tput sgr0)

if [ -z ${GH_TOKEN+x} ]; then
    echo "${bold}\$GH_TOKEN is not set!${normal}"
    exit 1
elif [ -z ${THEME_DIR+x} ]; then
    echo "${bold}\$THEME_DIR is not set!${normal}"
    exit 1
elif [ -z ${DOCS_DIR+x} ]; then
    echo "${bold}\$DOCS_DIR is not set!${normal}"
    exit 1
fi

git clone -b cclib \
    "https://${GH_TOKEN}@github.com/cclib/sphinx_rtd_theme.git" \
    "${THEME_DIR}"/sphinx_rtd_theme

install -dm755 "${DOCS_DIR}"/_build/html
touch "${DOCS_DIR}"/_build/html/.nojekyll

pushd "${DOCS_DIR}"/..
make
popd
