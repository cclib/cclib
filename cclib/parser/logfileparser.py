# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
"""Generic output file parser and related tools"""


import bz2
import fileinput
import gzip
import inspect
import io
import logging
import os
import random
import sys
import zipfile
from abc import ABC, abstractmethod

import numpy

from cclib.parser import utils
from cclib.parser.data import ccData
from cclib.parser.data import ccData_optdone_bool


# This seems to avoid a problem with Avogadro.
logging.logMultiprocessing = 0


class myBZ2File(bz2.BZ2File):
    """Return string instead of bytes"""
    def __next__(self):
        line = super(bz2.BZ2File, self).__next__()
        return line.decode("ascii", "replace")

    def next(self):
        line = self.__next__()
        return line


class myGzipFile(gzip.GzipFile):
    """Return string instead of bytes"""
    def __next__(self):
        super_ob = super(gzip.GzipFile, self)
        # seemingly different versions of gzip can have either next or __next__
        if hasattr(super_ob, 'next'):
            line = super_ob.next()
        else:
            line = super_ob.__next__()
        return line.decode("ascii", "replace")

    def next(self):
        line = self.__next__()
        return line


class FileWrapper:
    """Wrap a file-like object or stream with some custom tweaks"""

    def __init__(self, source, pos=0):

        self.src = source

        # Most file-like objects have seek and tell methods, but streams returned
        # by urllib.urlopen in Python2 do not, which will raise an AttributeError
        # in this code. On the other hand, in Python3 these methods do exist since
        # urllib uses the stream class in the io library, but they raise a different
        # error, namely io.UnsupportedOperation. That is why it is hard to be more
        # specific with except block here.
        try:
            self.src.seek(0, 2)
            self.size = self.src.tell()
            self.src.seek(pos, 0)

        except (AttributeError, IOError, io.UnsupportedOperation):
            # Stream returned by urllib should have size information.
            if hasattr(self.src, 'headers') and 'content-length' in self.src.headers:
                self.size = int(self.src.headers['content-length'])
            else:
                self.size = pos

        # Assume the position is what was passed to the constructor.
        self.pos = pos
        
        self.last_line = None

    def next(self):
        line = next(self.src)
        self.pos += len(line)
        self.last_line = line
        return line

    def __next__(self):
        return self.next()

    def __iter__(self):
        return self

    def close(self):
        self.src.close()

    def seek(self, pos, ref):

        # If we are seeking to end, we can emulate it usually. As explained above,
        # we cannot be too specific with the except clause due to differences
        # between Python2 and 3. Yet another reason to drop Python 2 soon!
        try:
            self.src.seek(pos, ref)
        except:
            if ref == 2:
                self.src.read()
            else:
                raise

        if ref == 0:
            self.pos = pos
        if ref == 1:
            self.pos += pos
        if ref == 2 and hasattr(self, 'size'):
            self.pos = self.size


def openlogfile(filename, object=None):
    """Return a file object given a filename or if object specified decompresses it
    if needed and wrap it up.

    Given the filename or file object of a log file or a gzipped, zipped, or bzipped
    log file, this function returns a file-like object.

    Given a list of filenames, this function returns a FileInput object,
    which can be used for seamless iteration without concatenation.
    """

    # If there is a single string argument given.
    if type(filename) in [str, str]:

        extension = os.path.splitext(filename)[1]

        if extension == ".gz":
            fileobject = myGzipFile(filename, "r", fileobj=object)

        elif extension == ".zip":
            zip = zipfile.ZipFile(object, "r") if object else zipfile.ZipFile(filename, "r")
            assert len(zip.namelist()) == 1, "ERROR: Zip file contains more than 1 file"
            fileobject = io.StringIO(zip.read(zip.namelist()[0]).decode("ascii", "ignore"))

        elif extension in ['.bz', '.bz2']:
            # Module 'bz2' is not always importable.
            assert bz2 is not None, "ERROR: module bz2 cannot be imported"
            fileobject = myBZ2File(object, "r") if object else myBZ2File(filename, "r")

        else:
            # Assuming that object is text file encoded in utf-8
            fileobject = io.StringIO(object.decode('utf-8')) if object \
                    else FileWrapper(io.open(filename, "r", errors='ignore'))

        return fileobject

    elif hasattr(filename, "__iter__"):

        # This is needed, because fileinput will assume stdin when filename is empty.
        if len(filename) == 0:
            return None

        return fileinput.input(filename, openhook=fileinput.hook_compressed)


class Logfile(ABC):
    """Abstract class for logfile objects.

    Subclasses defined by cclib:
        ADF, DALTON, GAMESS, GAMESSUK, Gaussian, Jaguar, Molpro, MOPAC,
        NWChem, ORCA, Psi, Q-Chem
    """

    def __init__(self, source, loglevel=logging.ERROR, logname="Log",
                 logstream=sys.stderr, datatype=ccData_optdone_bool, **kwds):
        """Initialise the Logfile object.

        This should be called by a subclass in its own __init__ method.

        Inputs:
            source - a logfile, list of logfiles, or stream with at least a read method
            loglevel - integer corresponding to a log level from the logging module
            logname - name of the source logfile passed to this constructor
            logstream - where to output the logging information
            datatype - class to use for gathering data attributes
        """

        # Set the filename to source if it is a string or a list of strings, which are
        # assumed to be filenames. Otherwise, assume the source is a file-like object
        # if it has a read method, and we will try to use it like a stream.
        self.isfileinput = False
        if isinstance(source, str):
            self.filename = source
            self.isstream = False
        elif isinstance(source, list) and all([isinstance(s, str) for s in source]):
            self.filename = source
            self.isstream = False
        elif isinstance(source, fileinput.FileInput):
            self.filename = source
            self.isstream = False
            self.isfileinput = True
        elif hasattr(source, "read"):
            self.filename = f"stream {str(type(source))}"
            self.isstream = True
            self.stream = source
        else:
            raise ValueError("Unexpected source type.")

        # Set up the logger.
        # Note that calling logging.getLogger() with one name always returns the same instance.
        # Presently in cclib, all parser instances of the same class use the same logger,
        #   which means that care needs to be taken not to duplicate handlers.
        self.loglevel = loglevel
        self.logname = logname
        self.logger = logging.getLogger(f"{self.logname} {self.filename}")
        self.logger.setLevel(self.loglevel)
        if len(self.logger.handlers) == 0:
            handler = logging.StreamHandler(logstream)
            handler.setFormatter(logging.Formatter("[%(name)s %(levelname)s] %(message)s"))
            self.logger.addHandler(handler)

        # Set up the metadata.
        if not hasattr(self, "metadata"):
            self.metadata = {}
            self.metadata["package"] = self.logname
            self.metadata["methods"] = []
            # Indicate if the computation has completed successfully
            self.metadata['success'] = False


        # Periodic table of elements.
        self.table = utils.PeriodicTable()

        # This is the class that will be used in the data object returned by parse(), and should
        # normally be ccData or a subclass of it.
        self.datatype = datatype

        # Change the class used if we want optdone to be a list or if the 'future' option
        # is used, which might have more consequences in the future.
        optdone_as_list = kwds.get("optdone_as_list", False) or kwds.get("future", False)
        optdone_as_list = optdone_as_list if isinstance(optdone_as_list, bool) else False
        if optdone_as_list:
            self.datatype = ccData
        # Parsing of Natural Orbitals and Natural Spin Orbtials into one attribute
        self.unified_no_nso = kwds.get("future",False)

    def __setattr__(self, name, value):

        # Send info to logger if the attribute is in the list of attributes.
        if name in ccData._attrlist and hasattr(self, "logger"):

            # Call logger.info() only if the attribute is new.
            if not hasattr(self, name):
                if type(value) in [numpy.ndarray, list]:
                    self.logger.info(f"Creating attribute {name}[]")
                else:
                    self.logger.info(f"Creating attribute {name}: {str(value)}")

        # Set the attribute.
        object.__setattr__(self, name, value)

    def parse(self, progress=None, fupdate=0.05, cupdate=0.002):
        """Parse the logfile, using the assumed extract method of the child."""

        # Check that the sub-class has an extract attribute,
        #  that is callable with the proper number of arguemnts.
        if not hasattr(self, "extract"):
            raise AttributeError(
                f"Class {self.__class__.__name__} has no extract() method."
            )
        if not callable(self.extract):
            raise AttributeError(
                f"Method {self.__class__.__name__}._extract not callable."
            )
        if len(inspect.getfullargspec(self.extract)[0]) != 3:
            raise AttributeError(
                f"Method {self.__class__.__name__}._extract takes wrong number of arguments."
            )

        # Save the current list of attributes to keep after parsing.
        # The dict of self should be the same after parsing.
        _nodelete = list(set(self.__dict__.keys()))

        # Initiate the FileInput object for the input files.
        # Remember that self.filename can be a list of files.
        if not self.isstream:
            if not self.isfileinput:
                inputfile = openlogfile(self.filename)
            else:
                inputfile = self.filename
        else:
            inputfile = FileWrapper(self.stream)

        # Intialize self.progress
        is_compressed = isinstance(inputfile, myGzipFile) or isinstance(inputfile, myBZ2File)
        if progress and not (is_compressed):
            self.progress = progress
            self.progress.initialize(inputfile.size)
            self.progress.step = 0
        self.fupdate = fupdate
        self.cupdate = cupdate

        # Maybe the sub-class has something to do before parsing.
        self.before_parsing()

        # Loop over lines in the file object and call extract().
        # This is where the actual parsing is done.
        for line in inputfile:
            self.updateprogress(inputfile, "Unsupported information", cupdate)

            # This call should check if the line begins a section of extracted data.
            # If it does, it parses some lines and sets the relevant attributes (to self).
            # Any attributes can be freely set and used across calls, however only those
            #   in data._attrlist will be moved to final data object that is returned.
            try:
                self.extract(inputfile, line)
            except StopIteration:
                self.logger.error("Unexpectedly encountered end of logfile.")
                break
            except Exception as e:
                self.logger.error("Encountered error when parsing.")
                self.logger.error(f"Last line read: {inputfile.last_line}")
                raise

        # Close input file object.
        if not self.isstream:
            inputfile.close()

        # Maybe the sub-class has something to do after parsing.
        self.after_parsing()

        # If atomcoords were not parsed, but some input coordinates were ("inputcoords").
        # This is originally from the Gaussian parser, a regression fix.
        if not hasattr(self, "atomcoords") and hasattr(self, "inputcoords"):
            self.atomcoords = numpy.array(self.inputcoords, 'd')

        # Set nmo if not set already - to nbasis.
        if not hasattr(self, "nmo") and hasattr(self, "nbasis"):
            self.nmo = self.nbasis

        # Create a default coreelectrons array, unless it's impossible
        # to determine.
        if not hasattr(self, "coreelectrons") and hasattr(self, "natom"):
            self.coreelectrons = numpy.zeros(self.natom, "i")
        if hasattr(self, "incorrect_coreelectrons"):
            self.__delattr__("coreelectrons")

        # Create the data object we want to return. This is normally ccData, but can be changed
        # by passing the datatype argument to the constructor. All supported cclib attributes
        # are copied to this object, but beware that in order to be moved an attribute must be
        # included in the data._attrlist of ccData (or whatever else).
        # There is the possibility of passing assitional argument via self.data_args, but
        # we use this sparingly in cases where we want to limit the API with options, etc.
        data = self.datatype(attributes=self.__dict__)

        # Now make sure that the cclib attributes in the data object are all the correct type,
        # including arrays and lists of arrays.
        data.arrayify()

        # Delete all temporary attributes (including cclib attributes).
        # All attributes should have been moved to a data object, which will be returned.
        for attr in list(self.__dict__.keys()):
            if not attr in _nodelete:
                self.__delattr__(attr)

        # Perform final checks on values of attributes.
        data.check_values(logger=self.logger)

        # Update self.progress as done.
        if hasattr(self, "progress"):
            self.progress.update(inputfile.size, "Done")

        return data

    def before_parsing(self):
        """Set parser-specific variables and do other initial things here."""
        pass

    def after_parsing(self):
        """Correct data or do parser-specific validation after parsing is finished."""
        pass

    def updateprogress(self, inputfile, msg, xupdate=0.05):
        """Update progress."""

        if hasattr(self, "progress") and random.random() < xupdate:
            newstep = inputfile.pos
            if newstep != self.progress.step:
                self.progress.update(newstep, msg)
                self.progress.step = newstep

    @abstractmethod
    def normalisesym(self, symlabel):
        """Standardise the symmetry labels between parsers."""

    def new_internal_job(self):
        """Delete attributes that can be problematic in multistep jobs.

        TODO: instead of this hack, parse each job in a multistep comptation
        as a different ccData object (this is for 2.x).

        Some computations are actually sequences of several jobs, and some
        attributes won't work well if parsed across jobs. There include:
            mpenergies: if different jobs go to different orders then
                        these won't be consistent and can't be converted
                        to an array easily
        """
        for name in ("mpenergies",):
            if hasattr(self, name):
                delattr(self, name)

    def set_attribute(self, name, value, check_change=True):
        """Set an attribute and perform an optional check when it already exists.

        Note that this can be used for scalars and lists alike, whenever we want
        to set a value for an attribute.
        
        Parameters
        ----------
        name: str
            The name of the attribute.
        value: str
            The value for the attribute.
        check_change: bool
            By default we want to check that the value does not change
            if the attribute already exists.
        """
        if check_change and hasattr(self, name):
            try:
                numpy.testing.assert_equal(getattr(self, name), value)
            except AssertionError:
                self.logger.warning(
                    f"Attribute {name} changed value ({getattr(self, name)} -> {value})"
                )

        setattr(self, name, value)

    def append_attribute(self, name, value):
        """Appends a value to an attribute."""

        if not hasattr(self, name):
            self.set_attribute(name, [])
        getattr(self, name).append(value)

    def extend_attribute(self, name, values):
        """Appends an iterable of values to an attribute."""
        
        if not hasattr(self, name):
            self.set_attribute(name, [])
        getattr(self, name).extend(values)

    def _assign_coreelectrons_to_element(self, element, ncore,
                                         ncore_is_total_count=False):
        """Assign core electrons to all instances of the element.

        It's usually reasonable to do this for all atoms of a given element,
        because mixed usage isn't normally allowed within elements.

        Parameters
        ----------
        element: str
          the chemical element to set coreelectrons for
        ncore: int
          the number of core electrons
        ncore_is_total_count: bool
          whether the ncore argument is the total count, in which case it is
          divided by the number of atoms of this element
        """
        atomsymbols = [self.table.element[atomno] for atomno in self.atomnos]
        indices = [i for i, el in enumerate(atomsymbols) if el == element]
        if ncore_is_total_count:
            ncore = ncore // len(indices)

        if not hasattr(self, 'coreelectrons'):
            self.coreelectrons = numpy.zeros(self.natom, 'i')
        self.coreelectrons[indices] = ncore

    def skip_lines(self, inputfile, sequence):
        """Read trivial line types and check they are what they are supposed to be.

        This function will read len(sequence) lines and do certain checks on them,
        when the elements of sequence have the appropriate values. Currently the
        following elements trigger checks:
            'blank' or 'b'      - the line should be blank
            'dashes' or 'd'     - the line should contain only dashes (or spaces)
            'equals' or 'e'     - the line should contain only equal signs (or spaces)
            'stars' or 's'      - the line should contain only stars (or spaces)
        """

        expected_characters = {
            '-': ['dashes', 'd'],
            '=': ['equals', 'e'],
            '*': ['stars', 's'],
        }

        lines = []
        for expected in sequence:

            # Read the line we want to skip.
            line = next(inputfile)

            # Blank lines are perhaps the most common thing we want to check for.
            if expected in ["blank", "b"]:
                try:
                    assert line.strip() == ""
                except AssertionError:
                    frame, fname, lno, funcname, funcline, index = inspect.getouterframes(inspect.currentframe())[1]
                    parser = fname.split('/')[-1]
                    msg = f"In {parser}, line {int(lno)}, line not blank as expected: {line.strip()}"
                    self.logger.warning(msg)

            # All cases of heterogeneous lines can be dealt with by the same code.
            for character, keys in expected_characters.items():
                if expected in keys:
                    try:
                        assert utils.str_contains_only(line.strip(), [character, ' '])
                    except AssertionError:
                        frame, fname, lno, funcname, funcline, index = inspect.getouterframes(inspect.currentframe())[1]
                        parser = fname.split('/')[-1]
                        msg = f"In {parser}, line {int(lno)}, line not all {keys[0]} as expected: {line.strip()}"
                        self.logger.warning(msg)
                        continue

            # Save the skipped line, and we will return the whole list.
            lines.append(line)

        return lines

    skip_line = lambda self, inputfile, expected: self.skip_lines(inputfile, [expected])

