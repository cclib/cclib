#!/usr/bin/env bash

# build_docs.sh: Build the documentation as part of Travis CI.

set -o errexit

bold=$(tput bold)
normal=$(tput sgr0)

if [ -z ${DOCS_BUILD_DIR+x} ]; then
    echo "${bold}\$DOCS_BUILD_DIR is not set!${normal}"
    exit 1
elif [ -z ${THEME_DIR+x} ]; then
    echo "${bold}\$THEME_DIR is not set!${normal}"
    exit 1
fi

git clone -b cclib \
    "https://git@github.com/cclib/sphinx_rtd_theme.git" \
    "${THEME_DIR}"/sphinx_rtd_theme

install -dm755 "${DOCS_BUILD_DIR}"
touch "${DOCS_BUILD_DIR}"/.nojekyll

pushd "${DOCS_BUILD_DIR}"/../..
make default
popd

TEST_COVERAGE_DIR="${TRAVIS_BUILD_DIR}/htmlcov"
if [[ -d "${TEST_COVERAGE_DIR}" ]]; then
    echo "${bold}Copying test coverage results...${normal}"
    mv "${TEST_COVERAGE_DIR}" "${DOCS_BUILD_DIR}"/coverage
fi
