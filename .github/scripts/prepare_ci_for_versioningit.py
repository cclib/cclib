#!/usr/bin/env python

# Copyright (c) 2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

# prepare_ci_for_versioningit.py: The sed command that updates the Open Babel
# dependency in CI dirties the working copy, which would mess with the version
# string in the build artifacts and cause problems when attempting to upload
# to a package index.  Disable this check since this workaround is not for the
# Python version that saves build artifacts and uploads to package indices.

import os
import sys


minor_version = sys.version_info.minor
if minor_version >= 10:
    with open(os.environ["GITHUB_ENV"], "a") as github_env_file:
        github_env_file.write("VERSIONINGIT_FOR_PACKAGE_INDEX=1\n")
