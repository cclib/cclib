# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
import os
from pathlib import Path
from typing import Any, MutableMapping

import requests

# adapted from https://stackoverflow.com/a/46271487
API_TOKEN_ENVVAR = "GITHUB_API_TOKEN"
API_TOKEN = os.getenv(API_TOKEN_ENVVAR, None)
if not API_TOKEN:
    raise RuntimeError(f"GitHub API token not stored in env variable ${API_TOKEN_ENVVAR}")
URL = "https://api.github.com/graphql"


def execute_query(queryfile: Path) -> str:
    """Execute a GraphQL query against the GitHub API."""
    assert queryfile.exists()
    assert queryfile.is_file()
    query = {"query": queryfile.read_text(encoding="utf-8")}
    headers = {"Authorization": f"token {API_TOKEN}"}
    result = requests.post(url=URL, json=query, headers=headers)
    return result.text


def transform_author(item: MutableMapping[str, Any]) -> None:
    """Flatten the author field on an item.

    We have no need for the additional fields in an Actor aside from the login name.
    https://docs.github.com/en/graphql/reference/interfaces#actor
    """
    item["author"] = item["author"]["login"]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("queryfile", type=Path)

    args = parser.parse_args()
