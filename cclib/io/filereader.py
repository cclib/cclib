# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Generic file reader and related tools"""

from abc import ABC, abstractmethod


class Reader(ABC):
    """Abstract class for reader objects."""

    def __init__(self, source, *args, **kwargs):
        """Initialize the Reader object.

        This should be called by a subclass in its own __init__ method.

        Inputs:
          source - A single filename, stream [TODO], or list of filenames/streams [TODO].
        """
        if isinstance(source, str):
            self.filename = source
        else:
            raise ValueError

    def parse(self):
        """Read the raw contents of the source into the Reader."""
        # TODO This cannot currently handle streams.
        with open(self.filename) as handle:
            self.filecontents = handle.read()

        return None

    @abstractmethod
    def generate_repr(self):
        """Convert the raw contents of the source into the internal representation."""
