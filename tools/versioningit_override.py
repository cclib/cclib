# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from os import getenv
from typing import Any, Dict

from versioningit import VCSDescription
from versioningit.basics import DEFAULT_FORMATS

_ENVVARNAME = "VERSIONINGIT_FOR_PACKAGE_INDEX"


def cclib_format(
    *, description: VCSDescription, base_version: str, next_version: str, params: Dict[str, Any]
) -> str:
    """Override the formation of the version string created by versioningit.

    The default version strings are not suitable for uploading to (Test)PyPI,
    which does not accept 'local' PEP 440 strings.  If the environment
    variable 'VERSIONINGIT_FOR_PACKAGE_INDEX' is set to a truthy value, a
    simplified version string is formed that counts the number of commits
    since the last tag.  This is safe since we only push to TestPyPI on
    commits (merges) to the default branch, and pushing to PyPI only happens
    for tags when the base version is used and this code is never even
    executed.  Because of this and it only making sense to run with this
    setting in continuous integration, 'dirty' builds are not allowed.

    When building normally, the version strings are the defaults but do not
    include mention of the VCS used (Git).
    """
    state = description.state
    assert state in {"distance", "dirty", "distance-dirty"}

    if getenv(_ENVVARNAME, "False").lower() in ("true", "1", "t"):
        fmt_distance = "{base_version}.post{distance}"
        if state != "distance":
            raise RuntimeError("dirty state doesn't make sense when building for a package index")
    else:
        # Default but missing {vcs} before {rev}
        fmt_distance = "{base_version}.post{distance}+{rev}"
        # Default
        fmt_dirty = DEFAULT_FORMATS["dirty"]
        # Default but missing {vcs} before {rev}
        fmt_distance_dirty = "{base_version}.post{distance}+{rev}.d{build_date:%Y%m%d}"

    if state == "distance":
        fmt = fmt_distance
    elif state == "dirty":
        fmt = fmt_dirty
    elif state == "distance-dirty":
        fmt = fmt_distance_dirty

    return fmt.format_map({**description.fields, "base_version": base_version})
