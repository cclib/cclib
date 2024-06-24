# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Classes for hierarchical storage of parsed data."""

from typing import TYPE_CHECKING, Any, Dict, List, Mapping, Optional

from cclib.attribute_parsers import ccData

if TYPE_CHECKING:
    from cclib.combinator import Combinator
    from cclib.tree import Tree


class ccCollection:
    def __init__(
        self, combinator: Optional["Combinator"] = None, tree: Optional["Tree"] = None
    ) -> None:
        """Initialize the ccCollection object.

        Inputs:
            attributes - optional dictionary of attributes to load as data
        """
        # we should be using the combinator object to figure out what the storage structure should be.
        # For now its a list of a ccData
        # [
        #         ccData[attr1, attr2]
        # ]

        self._combinator = combinator
        self._tree = tree
        self._parsed_data = [ccData() for i in range(len(self._tree))]
        # [
        #         [attrparser1, attrparser2, attrparser3, etc]
        # ]
        # if self._combinator != None:
        #    assert len(self._combinator.job_list) == 1

    @property
    def parsed_data(self) -> List[ccData]:
        return self._parsed_data
