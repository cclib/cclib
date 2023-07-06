# TOOD: This file belongs in cclib.io, but circular dependency issues mean it can't go there just now.

import bz2
import gzip
import zipfile
import pathlib
from tempfile import NamedTemporaryFile
from urllib.request import urlopen
from urllib.error import URLError
import collections
import typing
import re
import io
import logging
import codecs


# Regular expression for validating URLs
URL_PATTERN = re.compile(

    r'^(?:http|ftp)s?://'  # http:// or https://
    r'(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|'  # domain...
    r'localhost|'  # localhost...
    r'\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})'  # ...or ip
    r'(?::\d+)?'  # optional port
    r'(?:/?|[/?]\S+)$', re.IGNORECASE

)


def logerror(error):
    """
    Log a unicode decode/encode error to the logger and return a replacement character.
    """
    logging.warning(str(error))
    
    # Return type is a tuple.
    # First item is a replacement character. Second is the position to continue from.
    return (u'', error.start +1)
    
codecs.register_error('logerror', logerror)


class FileWrapper:
    """Wrap any supported input file type."""

    def __init__(self, *sources) -> None:
        # Each source can be a lot of different things, go through and process them now.
        self.input_files = {}
        self.file_pointer = 0
        
        # First, check if we were given an unpacked list for backwards compatibility.
        expanded_sources = []
        for source in sources:
            if isinstance(source, list):
                expanded_sources.extend(source)
            
            else:
                expanded_sources.append(source)
        
        for source in expanded_sources:
            # What are you?
            if isinstance(source, str) and URL_PATTERN.match(source):
                # This file is a URL.
                try:
                    # Cache the file to a temp location.
                    response = urlopen(source)
                    tfile = NamedTemporaryFile(delete = True)
                    tfile.write(response.read())
                    
                    fileobject = tfile
                    filename = source
                    
                except (ValueError, URLError) as error:
                    # Maybe no need to raise a different exception?
                    raise ValueError(
                        "Encountered an error processing the URL '{}'".format(source)
                    ) from error
                
            elif hasattr(source, "read") or hasattr(source, "readline"):
                # This file is a file.
                # We don't need to do anything.
                fileobject = source
                # Not sure what we're supposed to do if we don't have an available file name,
                # probably ok.
                filename = getattr(source, "name", f"stream {str(type(source))}")
                
            else:
                fileobject = None
                filename = source
                
            # Now 'open' the file. If the file is compressed, this function will uncompress it.
            # Likewise, appropriate decoding and error handling will be applied.
            #
            # If no file has been opened yet (source is a string-like), open it.
            self.input_files[filename] = self.open_log_file(filename, fileobject = fileobject)
        
        # TODO: Implement this for progress updating.
        self.size = None
        
        # A short buffer of previously read lines.
        # This permits primitive 'look-behind' functionality in some parsers (see Turbomole).
        # Init with 10 empty strings (empty lines).
        self.last_lines = collections.deque([""] * 10, 10)
        
    @property
    def file_name(self):
        return ", ".join(self.input_files)
    
    @classmethod
    def open_log_file(
            self,
            filename: str,
            mode: str = "r",
            encoding: str = "utf-8",
            errors: str = "logerror",
            fileobject: typing.Optional[typing.IO] = None
        ) -> typing.IO:
        """
        Open a possibly compressed file.
        """
        filename = pathlib.Path(filename)
        extension = filename.suffix
    
        if extension == ".gz":
            fileobject = io.TextIOWrapper(gzip.GzipFile(filename, mode, fileobj = fileobject), encoding = encoding, errors = errors)
    
        elif extension == ".zip":
            fileobject = zipfile.ZipFile(fileobject if fileobject else filename, mode)
            # TODO: Need to check that we're not leaving any open file objects here...
            # TODO: We should be able to handle multiple files...
            assert len(fileobject.namelist()) == 1, "ERROR: Zip file contains more than 1 file"
            
            fileobject = io.TextIOWrapper(
                fileobject.open(fileobject.namelist()[0]),
                encoding = encoding, errors = errors
            )
    
        elif extension in ['.bz', '.bz2']:
            # Module 'bz2' is not always importable.
            assert bz2 is not None, "ERROR: module bz2 cannot be imported"
            fileobject = io.TextIOWrapper(bz2.BZ2File(fileobject if fileobject else filename, mode), encoding = encoding, errors = errors)
    
        elif fileobject is not None:
            # Assuming that object is text file encoded in utf-8
            #fileobject = io.StringIO(fileobject.decode(encoding, errors))
            fileobject = io.TextIOWrapper(fileobject, encoding = encoding, errors = errors)
            
        else:
            # Normal text file.
            
            fileobject = open(filename, mode, encoding = encoding, errors = errors)
        
        return fileobject

    def next(self):
        try:
            # TODO: Wasteful to make a list each iteration here...
            try:
                line = next(list(self.input_files.values())[self.file_pointer])
                self.last_lines.append(line)
                return line
            
            except StopIteration:
                self.file_pointer += 1
                return self.next()
            
        except IndexError:
            raise StopIteration()
    
    @property
    def last_line(self):
        return self.last_lines[-1]

    def __next__(self):
        return self.next()

    def __iter__(self):
        return self
        # This works great, but means we loose support for next() on this class.
#         for file in self.input_files:
#             for line in file:
#                 self.last_lines.append(line)        
#                 yield line

    def readline(self) -> str:
        """
        Read one line from this file.
        """
        return next(self)
    
    def read(self) -> str:
        """
        Read everything from this file.
        
        Be aware that this function will load the entire file into a single string.
        """
        return "".join(list(self))

    def close(self) -> None:
        for input_file in self.input_files:
            input_file.close()

    def seek(self, pos: int, ref: int) -> None:
        raise NotImplementedError("FileWrapper does not support seek()")
    
    def reset(self):
        # Equivalent to seeking to 0 for all our files.
        for file in self.input_files.values():
            file.seek(0,0)
            
        self.file_pointer = 0
    
#     def peek_ahead(self):
#         """
#         Acquire a copy of this file object to read lines without advancing
#         the file pointer of the original object. 
#         """