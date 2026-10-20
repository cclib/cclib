# Copyright (c) 2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from pathlib import Path

from cclib.io import ccread


DATA = Path(__file__).resolve().parent / ".." / ".." / "data" / "NWChem" / "basicNWChem7.0"


def test_basis_summary_multiword_description(tmp_path: Path) -> None:
    original = (DATA / "C_bigbasis.out").read_text()
    summary = " C                        aug-cc-pVQZ               20      105   6s5p4d3f2g"
    assert summary in original
    modified = original.replace(
        summary, " C                       user specified             20      105   6s5p4d3f2g"
    )
    logfile = tmp_path / "C_user_specified.out"
    logfile.write_text(modified)

    data = ccread(str(logfile))

    assert data.metadata["basis_set"] == "user specified"
    assert len(data.atombasis[0]) == 105
