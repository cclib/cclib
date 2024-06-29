# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
"""Generic output file parser and related tools"""

import inspect
import logging
import random
import sys
import typing
from abc import ABC, abstractmethod
from typing import Any, Iterable, List, Optional

from cclib.parser import utils
from cclib.parser.data import ccData
from cclib.parser.logfilewrapper import FileWrapper

import numpy

# This seems to avoid a problem with Avogadro.
logging.logMultiprocessing = 0


class StopParsing(Exception):
    """
    An exception to signal that parsing should stop.
    """


class Logfile(ABC):
    """Abstract class for logfile objects."""

    def __init__(
        self,
        source: typing.Union[
            str, typing.IO, FileWrapper, typing.List[typing.Union[str, typing.IO]]
        ],
        loglevel: int = logging.ERROR,
        logname: str = "Log",
        logstream=sys.stderr,
        datatype=ccData,
        **kwds,
    ):
        """Initialise the Logfile object.

        This should be called by a subclass in its own __init__ method.

        Inputs:
            source - a logfile, list of logfiles, or stream with at least a read method
            loglevel - integer corresponding to a log level from the logging module
            logname - name of the logging object to use for this parser
            logstream - where to output the logging information
            datatype - class to use for gathering data attributes
        """
        if not isinstance(source, FileWrapper):
            source = FileWrapper(source)
            # Probably the wrong type given.
            # raise TypeError("Source does not have an 'input_files' attribute, are you sure it inherits from FileWrapper?") from None

        self.inputfile = source
        # If our parser needs a certain file ordering, set that now.
        self.inputfile.sort(self.sort_input(self.inputfile.filenames))

        # Set up the logger.
        # Note that calling logging.getLogger() with one name always returns the same instance.
        # Presently in cclib, all parser instances of the same class use the same logger,
        #   which means that care needs to be taken not to duplicate handlers.
        self.loglevel = loglevel
        self.logname = logname
        self.logger = logging.getLogger(f"{self.logname} {self.inputfile.file_name}")
        self.logger.setLevel(self.loglevel)
        if len(self.logger.handlers) == 0:
            handler = logging.StreamHandler(logstream)
            handler.setFormatter(logging.Formatter("[%(name)s %(levelname)s] %(message)s"))
            self.logger.addHandler(handler)

        # Set up the metadata.
        if not hasattr(self, "metadata"):
            self.metadata = {"package": self.logname, "methods": []}

        # Periodic table of elements.
        self.table = utils.PeriodicTable()

        # This is the class that will be used in the data object returned by parse(), and should
        # normally be ccData or a subclass of it.
        self.datatype = datatype

        self.future = kwds.get("future", False)
        # Parsing of Natural Orbitals and Natural Spin Orbtials into one attribute
        self.unified_no_nso = self.future

    @property
    def filename(self):
        return self.inputfile.file_name

    @classmethod
    def sort_input(self, file_names: typing.List[str]) -> typing.List:
        """
        If this parser expects multiple files to appear in a certain order, return that ordering.
        """
        return file_names

    def __setattr__(self, name, value):
        # Send info to logger if the attribute is in the list of attributes.
        if name in ccData._attrlist and hasattr(self, "logger"):
            # Call logger.info() only if the attribute is new.
            if not hasattr(self, name):
                if type(value) in [numpy.ndarray, list]:
                    self.logger.info("Creating attribute %s[]", name)
                else:
                    self.logger.info("Creating attribute %s: %s", name, value)

        # Set the attribute.
        object.__setattr__(self, name, value)

    def parse(self, progress=None, fupdate=0.05, cupdate=0.002):
        """Parse the logfile, using the assumed extract method of the child."""

        # Check that the sub-class has an extract attribute,
        #  that is callable with the proper number of arguemnts.
        if not hasattr(self, "extract"):
            raise AttributeError(f"Class {self.__class__.__name__} has no extract() method.")
        if not callable(self.extract):
            raise AttributeError(f"Method {self.__class__.__name__}._extract not callable.")
        if len(inspect.getfullargspec(self.extract)[0]) != 3:
            raise AttributeError(
                f"Method {self.__class__.__name__}._extract takes wrong number of arguments."
            )

        # Save the current list of attributes to keep after parsing.
        # The dict of self should be the same after parsing.
        _nodelete = list(set(self.__dict__.keys()))

        # Intialize self.progress
        if progress:
            self.progress = progress
            self.progress.initialize(self.inputfile.size)
            self.progress.step = 0
        self.fupdate = fupdate
        self.cupdate = cupdate

        # Maybe the sub-class has something to do before parsing.
        self.before_parsing()

        # Loop over lines in the file object and call extract().
        # This is where the actual parsing is done.
        for line in self.inputfile:
            self.updateprogress(self.inputfile, "Unsupported information", cupdate)

            # This call should check if the line begins a section of extracted data.
            # If it does, it parses some lines and sets the relevant attributes (to self).
            # Any attributes can be freely set and used across calls, however only those
            #   in data._attrlist will be moved to final data object that is returned.
            try:
                self.extract(self.inputfile, line)
            except StopParsing:
                # This is fine
                break
            except StopIteration:
                self.logger.error("Unexpectedly encountered end of logfile.")
                break
            except Exception:
                self.logger.error("Encountered error when parsing.")

                # Not all input files support last_line.
                if hasattr(self.inputfile, "last_line"):
                    self.logger.error("Last line read: %s", self.inputfile.last_line)
                raise

        # Maybe the sub-class has something to do after parsing.
        self.after_parsing()

        # If atomcoords were not parsed, but some input coordinates were ("inputcoords").
        # This is originally from the Gaussian parser, a regression fix.
        if not hasattr(self, "atomcoords") and hasattr(self, "inputcoords"):
            self.atomcoords = numpy.array(self.inputcoords, "d")

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
            if attr not in _nodelete:
                self.__delattr__(attr)

        # Perform final checks on values of attributes.
        data.check_values(logger=self.logger)

        # Update self.progress as done.
        if hasattr(self, "progress"):
            self.progress.update(self.inputfile.size, "Done")

        return data

    def before_parsing(self) -> None:
        """Set parser-specific variables and do other initial things here."""
        pass

    def after_parsing(self) -> None:
        """Correct data or do parser-specific validation after parsing is finished."""
        if (
            hasattr(self, "enthalpy")
            and hasattr(self, "entropy")
            and hasattr(self, "temperature")
            and not hasattr(self, "freeenergy")
        ):
            self.set_attribute("freeenergy", self.enthalpy - self.entropy * self.temperature)

    def updateprogress(self, inputfile, msg: str, xupdate: float = 0.05) -> None:
        """Update progress."""

        if hasattr(self, "progress") and random.random() < xupdate:
            newstep = inputfile.pos
            if newstep != self.progress.step:
                self.progress.update(newstep, msg)
                self.progress.step = newstep

    @abstractmethod
    def normalisesym(self, symlabel: str) -> None:
        """Standardise the symmetry labels between parsers."""

    def new_internal_job(self) -> None:
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

    def hasattrs(self, names: Iterable[str]) -> bool:
        """Does this logfile have all the given attributes?"""
        return all(hasattr(self, name) for name in names)

    def set_attribute(self, name: str, value: Any, check_change: bool = True) -> None:
        """Set an attribute and perform an optional check when it already exists.

        Note that this can be used for scalars and lists alike, whenever we want
        to set a value for an attribute.

        Parameters
        ----------
        name: str
            The name of the attribute.
        value: any
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
                    "Attribute %s changed value (%s -> %s)", name, getattr(self, name), value
                )

        setattr(self, name, value)

    def del_attribute(self, name):
        """Safely remove/delete an attribute, even if it does not exist."""
        try:
            delattr(self, name)

        except AttributeError:
            pass

    def append_attribute(self, name: str, value: Any) -> None:
        """Appends a value to an attribute."""

        if not hasattr(self, name):
            self.set_attribute(name, [])
        getattr(self, name).append(value)

    def extend_attribute(
        self, name: str, values: Iterable[Any], index: Optional[int] = None
    ) -> None:
        """Appends an iterable of values to an attribute."""

        if not hasattr(self, name):
            self.set_attribute(name, [])

        if isinstance(getattr(self, name), list) and index is not None:
            getattr(self, name)[index].extend(values)
        else:
            getattr(self, name).extend(values)

    def _assign_coreelectrons_to_element(
        self, element: str, ncore: int, ncore_is_total_count: bool = False
    ) -> None:
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

        if not hasattr(self, "coreelectrons"):
            self.coreelectrons = numpy.zeros(self.natom, "i")
        self.coreelectrons[indices] = ncore

    def skip_lines(self, inputfile, sequence: Iterable[str]) -> List[str]:
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
            line = next(inputfile)

            # Blank lines are perhaps the most common thing we want to check for.
            if expected in ["blank", "b"]:
                try:
                    assert line.strip() == ""
                except AssertionError:
                    frame, fname, lno, funcname, funcline, index = inspect.getouterframes(
                        inspect.currentframe()
                    )[1]
                    parser = fname.split("/")[-1]
                    self.logger.warning(
                        "In %s, line %d, line not blank as expected: %s", parser, lno, line.strip()
                    )

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
                        self.logger.warning(
                            "In %s, line %d, line not all %s as expected: %s",
                            parser,
                            lno,
                            keys[0],
                            line.strip(),
                        )
                        continue

            # Save the skipped line, and we will return the whole list.
            lines.append(line)

        return lines

    @staticmethod
    def next_filled_line(inputfile):
        """Return the next line that contains something other than whitespace."""
        while True:
            line = next(inputfile)

            if line.strip() != "":
                return line

    def skip_line(self, inputfile: "FileWrapper", expected: str) -> List[str]:
        return self.skip_lines(inputfile, [expected])
