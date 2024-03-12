# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Classes for hierarchical storage of parsed data."""

import logging
from collections import namedtuple
from typing import Any, Dict, List, Mapping, Optional

from cclib.parser import ccData

import numpy


class ccCollection:
    def __init__(self, combinator=None, tree=None) -> None:
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
        self._parsed_data = [ccData() for i in range(self._tree.num_nodes)]
        # [
        #         [attrparser1, attrparser2, attrparser3, etc]
        # ]
        # if self._combinator != None:
        #    assert len(self._combinator.job_list) == 1
