# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Canonical exception classes for cclib."""


class MissingAttributeError(Exception):
    """Raised when a required ccData attribute is absent.

    This is the single canonical definition for the package.
    It replaces the local copies that were previously defined in
    each bridge and method module.
    """

    pass
