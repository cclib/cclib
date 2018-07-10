#!/usr/bin/env bash

# Run this to build the HTML documentation using Sphinx, then commit
# and push it to the docs repo.

set -o errexit

GH_USER=cclib
GH_REPO_NAME=cclib.github.io
DOCS_BRANCH_NAME=gh-pages
GH_REPO_REF=github.com:$GH_USER/$GH_REPO_NAME.git

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

pushd "${SCRIPTDIR}"/../doc

# Make sure a copy of the custom theme is present and up-to-date.
mkdir -p sphinx/_themes
if [ ! -d sphinx/_themes/sphinx_rtd_theme ]; then
    pushd sphinx/_themes/
    git clone -b cclib git@github.com:cclib/sphinx_rtd_theme.git
else
    pushd sphinx/_themes/sphinx_rtd_theme
    git pull
fi
popd

# Build the HTML documentation.
make

if [ ! -d $GH_REPO_NAME ]; then
    echo "Cloning $DOCS_BRANCH_NAME branch of $GH_REPO_REF..."
    git clone -b $DOCS_BRANCH_NAME git@$GH_REPO_REF
fi
cd $GH_REPO_NAME
rm -rf ./* ./.*
echo "" > .nojekyll
echo "Copying built HTML..."
cp -a ../sphinx/_build/html/* .

echo "Adding changes..."
git add --all
echo "Committing..."
# This will return 1 if there are no changes, which should not result
# in failure.
git commit -m "Deploy documentation `date`" || ret=$?
git push --force origin $DOCS_BRANCH_NAME

popd
