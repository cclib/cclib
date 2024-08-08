# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

# TOOD: This file belongs in cclib.io, but circular dependency issues mean it can't go there just now.

import bz2
import codecs
import collections
import gzip
import inspect
import io
import logging
import pathlib
import re
import zipfile
from tempfile import NamedTemporaryFile
from typing import IO, Any, Iterable, List, Optional, Tuple
from urllib.error import URLError
from urllib.request import urlopen

from cclib.parser import utils

# Regular expression for validating URLs
URL_PATTERN = re.compile(
    r"^(?:http|ftp)s?://"  # http:// or https://
    r"(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|"  # domain...
    r"localhost|"  # localhost...
    r"\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})"  # ...or ip
    r"(?::\d+)?"  # optional port
    r"(?:/?|[/?]\S+)$",
    re.IGNORECASE,
)


def logerror(error):
    """
    Log a unicode decode/encode error to the logger and return a replacement character.
    """
    logging.getLogger("cclib").warning(str(error))

    # Return type is a tuple.
    # First item is a replacement character. Second is the position to continue from.
    return ("", error.start + 1)


codecs.register_error("logerror", logerror)


class FileHandler:
    """Wrap any supported input file type."""

    def __init__(self, *sources) -> None:
        # The total size of all our files.
        self.size = 0
        # Current 'byte' position in all our files.
        self.pos = 0
        # The current file position.
        self.file_pointer = 0

        self.filenames = []
        self.files = []
        self.sizes = []

        # First, check if we were given an unpacked list for backwards compatibility.
        expanded_sources = []
        for source in sources:
            if isinstance(source, list):
                expanded_sources.extend(source)
            else:
                expanded_sources.append(source)

        # Each source can be a lot of different things, go through and process them now.
        for source in expanded_sources:
            # 'open' the file. If the file is compressed, this function will uncompress it.
            # Likewise, appropriate decoding and error handling will be applied.
            #
            # If no file has been opened yet (source is a string-like), open it.
            filename, fileobject = self.open_log_file(source)

            # Move to the end of the file to determine how big it is.
            fileobject.seek(0, 2)
            size = fileobject.tell()
            fileobject.seek(0, 0)

            self.files.append(fileobject)
            # open_log_file returns a pathlib.Path object, cast to str for compatibility.
            self.filenames.append(str(filename))
            self.sizes.append(size)

        # Total number of bytes in all our files.
        self.size = sum(self.sizes)

        # A short buffer of previously read lines.
        # This permits primitive 'look-behind' functionality in some parsers (see Turbomole).
        # Init with 10 empty strings (empty lines).
        self.last_lines = collections.deque([""] * 10, 10)

        self.virtual_reset_position: Optional[int] = None

    @property
    def file_name(self) -> str:
        return ", ".join(self.filenames)

    def sort(self, order: list) -> None:
        """
        Sort the individual files that make up this FileHandler.

        order is an ordered list of filenames.
        """
        # Reset all our streams.
        self.reset()

        filenames = []
        files = []
        sizes = []

        # Reorder.
        for filename in order:
            index = self.filenames.index(filename)
            filenames.append(filename)
            files.append(self.files[index])
            sizes.append(self.sizes[index])

        self.filenames = filenames
        self.files = files
        self.sizes = sizes

    @classmethod
    def open_log_file(
        self, source, mode: str = "r", encoding: str = "utf-8", errors: str = "logerror"
    ) -> Tuple[str, IO]:
        """
        Open a possibly compressed file, returning both the filename of the file and an open file object.
        """
        # First, work out what source is (could be a filename, a URL, an open file etc).
        if isinstance(source, str) and URL_PATTERN.match(source):
            # This file is a URL.
            try:
                # Cache the file to a temp location.
                response = urlopen(source)
                fileobject = NamedTemporaryFile(delete=True)
                fileobject.write(response.read())
                fileobject.seek(0, 0)

                fileobject = io.TextIOWrapper(fileobject, encoding=encoding, errors=errors)
                filename = source

            except (ValueError, URLError) as error:
                # Maybe no need to raise a different exception?
                raise ValueError(
                    "Encountered an error processing the URL '{}'".format(source)
                ) from error

        elif hasattr(source, "read") or hasattr(source, "readline"):
            # This file is a file.
            # If this file supports seek, we don't need to do anything.
            # If not, we'll cache it to file.
            if not hasattr(source, "seek") or (
                hasattr(source, "seekable") and not source.seekable()
            ):
                fileobject = NamedTemporaryFile(delete=True)
                fileobject.write(source.read())
                fileobject.seek(0, 0)

                fileobject = io.TextIOWrapper(fileobject, encoding=encoding, errors=errors)

            else:
                fileobject = source
            filename = getattr(source, "name", f"stream {str(type(source))}")

        else:
            # This file is something else, assume we can open() it.
            filename = source
            fileobject = None

        filename = pathlib.Path(filename)
        extension = filename.suffix

        if extension == ".gz":
            fileobject = io.TextIOWrapper(
                gzip.GzipFile(filename, mode, fileobj=fileobject), encoding=encoding, errors=errors
            )

        elif extension == ".zip":
            fileobject = zipfile.ZipFile(fileobject if fileobject else filename, mode)
            # TODO: Need to check that we're not leaving any open file objects here...
            # TODO: We should be able to handle multiple files...
            assert len(fileobject.namelist()) == 1, "ERROR: Zip file contains more than 1 file"

            fileobject = io.TextIOWrapper(
                fileobject.open(fileobject.namelist()[0]), encoding=encoding, errors=errors
            )

        elif extension in [".bz", ".bz2"]:
            # Module 'bz2' is not always importable.
            assert bz2 is not None, "ERROR: module bz2 cannot be imported"
            fileobject = io.TextIOWrapper(
                bz2.BZ2File(fileobject if fileobject else filename, mode),
                encoding=encoding,
                errors=errors,
            )

        elif fileobject is not None:
            # Assuming that object is text file encoded in utf-8
            # If the file/stream has already been opened, we have no ability to handle decoding errors.
            pass

        else:
            # Normal text file.

            fileobject = open(filename, mode, encoding=encoding, errors=errors)

        return filename, fileobject

    def __enter__(self):
        """
        Enter context.
        """
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """
        Exit context
        """
        self.close()

    def virtual_set(self) -> None:
        """
        For use with property parsers to handle parsing
        through a block of text without changing the state of the handler.

        Sets the virtual_file_pointer to the current file_pointer
        """
        self.virtual_reset_position = self.pos

    def virtual_reset(self) -> None:
        """
        For use with property parsers to handle parsing
        through a block of text without changing the state of the handler.

        Sets the virtual_file_pointer to the current file_pointer
        """
        if self.virtual_reset_position is None:
            raise RuntimeError(
                "virtual_set() must be called before virtual_reset() and virtual_next()"
            )
        self.files[self.file_pointer].seek(self.virtual_reset_position)
        self.pos = self.virtual_reset_position
        self.virtual_reset_position = None

    def virtual_next(self) -> str:
        """
        For use with property parsers to handle parsing
        through a block of text without changing the state of the handler

        virtual_next will _not_ advance to the next file
        """
        if self.virtual_reset_position is None:
            raise RuntimeError("virtual_set() must be called before reset and virtual_next")
        try:
            line = next(self.files[self.file_pointer])
            self.pos += len(line)
            return line
        except:
            # possibly raise a warning? but we are ok just reaching the end of a file for a subparser parsing
            return

    def __next__(self) -> str:
        """
        Get the next line from this log file.
        """
        try:
            try:
                line = next(self.files[self.file_pointer])
                self.last_lines.append(line)
                self.pos += len(line)
                return line

            except StopIteration:
                self.file_pointer += 1
                return next(self)

        except IndexError:
            # raise StopIteration()
            return False

    @property
    def last_line(self) -> str:
        """
        Return the last line read by this parser.
        """
        return self.last_lines[-1]

    def __iter__(self):
        return self

    # TODO: support size parameter.
    def readline(self) -> str:
        """
        Read one line from this file.
        """
        return next(self)

    def readlines(self) -> List[str]:
        """
        Read all the lines from this file.
        """
        return list(self)

    # TODO: support size parameter.
    def read(self) -> str:
        """
        Read everything from this file.

        Be aware that this function will load the entire file into a single string.
        """
        return "".join(list(self))

    def close(self) -> None:
        """
        Close all open files.
        """
        for file in self.files:
            file.close()

    def __del__(self) -> None:
        """
        Make sure to close any open files when we go out of scope.

        Note that there is no guarantee when or if this function will get called;
        user's should ensure to close their own files once they are finished with them.
        """
        self.close()

    def seek(self, offset: int, whence: int = 0) -> None:
        if offset != 0 or whence not in (0, 2):
            raise NotImplementedError("FileHandler only supports seeking to start or end")

        if whence == 0:
            self.reset()

        elif whence == 2:
            self.finish()

    #     def seek_from_current(self, offset):
    #         """
    #         Seek forwards or backwards based on the current position.
    #         """
    #         remainder = offset
    #         while remainder != 0:
    #             if remainder > 0:
    #                 # Seek forwards please.
    #                 file_size = self.sizes[self.file_pointer]

    def reset(self) -> None:
        # Equivalent to seeking to 0 for all our files.
        for file in self.files:
            file.seek(0, 0)

        self.file_pointer = 0
        self.pos = 0

    def finish(self) -> None:
        # Equivalent to seeking to 2 for all our files.
        for file in self.files:
            file.seek(0, 2)

        self.file_pointer = len(self.files) - 1
        self.pos = self.size

    def skip_lines(self, sequence: Iterable[str], virtual: bool = False) -> List[str]:
        """Read trivial line types and check they are what they are supposed to be.

        This function will read len(sequence) lines and do certain checks on them,
        when the elements of sequence have the appropriate values. Currently the
        following elements trigger checks:
            'blank' or 'b'      - the line should be blank
            'dashes' or 'd'     - the line should contain only dashes (or spaces)
            'equals' or 'e'     - the line should contain only equal signs (or spaces)
            'stars' or 's'      - the line should contain only stars (or spaces)
        """

        expected_characters = {"-": ["dashes", "d"], "=": ["equals", "e"], "*": ["stars", "s"]}

        lines = []
        for expected in sequence:
            # Read the line we want to skip.
            if virtual:
                line = self.virtual_next()
            else:
                line = next(self)

            # Blank lines are perhaps the most common thing we want to check for.
            if expected in ["blank", "b"]:
                try:
                    assert line.strip() == ""
                except AssertionError:
                    frame, fname, lno, funcname, funcline, index = inspect.getouterframes(
                        inspect.currentframe()
                    )[1]
                    parser = fname.split("/")[-1]
                    msg = (
                        f"In {parser}, line {int(lno)}, line not blank as expected: {line.strip()}"
                    )
                    # self.logger.warning(msg)

            # All cases of heterogeneous lines can be dealt with by the same code.
            for character, keys in expected_characters.items():
                if expected in keys:
                    try:
                        assert utils.str_contains_only(line.strip(), [character, " "])
                    except AssertionError:
                        frame, fname, lno, funcname, funcline, index = inspect.getouterframes(
                            inspect.currentframe()
                        )[1]
                        parser = fname.split("/")[-1]
                        msg = f"In {parser}, line {int(lno)}, line not all {keys[0]} as expected: {line.strip()}"
                        # self.logger.warning(msg)
                        continue

            # Save the skipped line, and we will return the whole list.
            lines.append(line)

        return lines
