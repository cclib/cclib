#!/usr/bin/env bash

# deploy_docs_travis.bash: Deploy the documentation to the website using
# Travis. Assume that the documentation has already been built.

# Adapted from https://gist.github.com/vidavidorra/548ffbcdae99d752da02

set -o errexit

# Required environment variables:
# - DOCS_BRANCH_NAME: name of the remote branch serving the book
# - DOCS_BUILD_DIR: directory where documentation HTML files live
# - DOCS_REPO_NAME: name of the remote repo
# - DOCS_REPO_OWNER: name of the remote repo's owner
# - GH_TOKEN: [secret] Personal Access Token

bold=$(tput bold)
normal=$(tput sgr0)

if [ -z ${DOCS_BRANCH_NAME+x} ]; then
    echo "${bold}\$DOCS_BRANCH_NAME is not set!${normal}"
    exit 1
elif [ -z ${DOCS_BUILD_DIR+x} ]; then
    echo "${bold}\$DOCS_BUILD_DIR is not set!${normal}"
    exit 1
elif [ -z ${DOCS_REPO_NAME+x} ]; then
    echo "${bold}\$DOCS_REPO_NAME is not set!${normal}"
    exit 1
elif [ -z ${DOCS_REPO_OWNER+x} ]; then
    echo "${bold}\$DOCS_REPO_OWNER is not set!${normal}"
    exit 1
elif [ -z ${GH_TOKEN+x} ]; then
    echo "${bold}\$GH_TOKEN is not set!${normal}"
    exit 1
fi

git config user.name "Travis CI User"
git config user.email "travis@travis-ci.org"

GH_REPO_REF="github.com/${DOCS_REPO_OWNER}/${DOCS_REPO_NAME}.git"

git clone -b "${DOCS_BRANCH_NAME}" https://git@"${GH_REPO_REF}"
pushd ./"${DOCS_REPO_NAME}"
rm -rf ./*
echo "Copying built HTML..."
cp -a "${DOCS_BUILD_DIR}"/* .

if [ -f "index.html" ]; then
    echo "" > .nojekyll
    echo "${bold}Adding changes...${normal}"
    git add --all
    echo "${bold}Committing...${normal}"
    # This will return 1 if there are no changes, which should not result in
    # failure.
    git commit \
        -m "Deploy code docs to GitHub Pages Travis build: ${TRAVIS_BUILD_NUMBER}" \
        -m "Commit: ${TRAVIS_COMMIT}" || ret=$?
    git push "https://${GH_TOKEN}@${GH_REPO_REF}"
else
    echo "" >&2
    echo "${bold}Warning: No documentation (html) files have been found!${normal}" >&2
    echo "${bold}Warning: Not going to push the documentation to GitHub!${normal}" >&2
    exit 1
fi

popd
