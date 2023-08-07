# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Generic file reader and related tools"""

from abc import ABC, abstractmethod
import typing

from cclib.parser.logfilewrapper import FileWrapper


class Reader(ABC):
    """Abstract class for reader objects."""

    def __init__(self,
        source: typing.Union[str, typing.IO, FileWrapper, typing.List[typing.Union[str, typing.IO]]],
        *args,
        **kwargs)-> None:
        """Initialize the Reader object.

        This should be called by a subclass in its own __init__ method.

        Inputs:
          source - A single filename, stream [TODO], or list of filenames/streams [TODO].
        """
        if not isinstance(source, FileWrapper):
            source = FileWrapper(source)
        
        self.inputfile = source

    def parse(self) -> None:
        """Read the raw contents of the source into the Reader."""
        self.filecontents = self.inputfile.read()

    @abstractmethod
    def generate_repr(self) -> None:
        """Convert the raw contents of the source into the internal representation."""
