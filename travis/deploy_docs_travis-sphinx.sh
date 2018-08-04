#!/usr/bin/env bash

# deploy_docs_travis-sphinx.sh: Deploy the built documentation to
# GitHub Pages using Travis CI and travis-sphinx.

MESSAGE="Deploy code docs to GitHub Pages Travis build: ${TRAVIS_BUILD_NUMBER} Commit: ${TRAVIS_COMMIT}"
travis-sphinx --verbose deploy --message="${MESSAGE}"
