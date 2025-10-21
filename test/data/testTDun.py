# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test single point, unrestricted time-dependent logfiles in cclib"""

from typing import TYPE_CHECKING

from skip import skipForLogfile, skipForParser

if TYPE_CHECKING:
    from cclib.parser.data import ccData


# The Gaussian log files for this test are a normal restricted calculation,
# is this class misnamed?
class GenericTDunTest:
    """Generic time-dependent unrestricted HF/DFT unittest"""

    number = 24

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    def testenergiesnumber(self, data: "ccData") -> None:
        """Is the length of etenergies correct?"""
        assert len(data.etenergies) == self.number

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForLogfile(
        "Turbomole/basicTurbomole7.4/CO_cc2_TD_un",
        "Oscillator strengths are not available for triplets with Turbomole's ricc2",
    )
    def testoscsnumber(self, data: "ccData") -> None:
        """Is the length of eotscs correct?"""
        assert len(data.etoscs) == self.number

    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("PySCF", "etrotats are not yet implemented")
    @skipForLogfile(
        "Turbomole/basicTurbomole7.4/CO_cc2_TD",
        "Rotatory strengths are not currently available for ricc2",
    )
    def testrotatsnumber(self, data: "ccData") -> None:
        """Is the length of etrotats correct?"""
        assert len(data.etrotats) == self.number

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    def testsecsnumber(self, data: "ccData") -> None:
        """Is the length of etsecs correct?"""
        assert len(data.etsecs) == self.number

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    def testsymsnumber(self, data: "ccData") -> None:
        """Is the length of etsyms correct?"""
        assert len(data.etsyms) == self.number

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "Turbomole etsyms are not available for UHF")
    def testsyms(self, data: "ccData") -> None:
        """Is etsyms populated by singlets and triplets 50/50?"""
        singlets = [sym for sym in data.etsyms if "Singlet" in sym]
        triplets = [sym for sym in data.etsyms if "Triplet" in sym]
        assert len(singlets) == self.number / 2
        assert len(triplets) == self.number / 2


class CFOUREOMCCSDunTest(GenericTDunTest):
    """Customized UHF/EOMEE-CCSD unittest"""

    number = 10


class TurbomoleTDunTest(GenericTDunTest):
    """Customized time-dependent unrestricted HF/DFT unittest"""

    number = 10


class PySCFTDunTest(GenericTDunTest):
    """Customized time-dependent unrestricted HF/DFT unittest"""

    number = 10
