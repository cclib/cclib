# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A writer for MolSSI quantum chemical JSON (QCJSON) files.
"""

import json

from cclib.parser.utils import convertor

from .cjsonwriter import CJSON as CJSONWriter
from .cjsonwriter import JSONIndentEncoder, NumpyAwareJSONEncoder


class QCJSONWriter(CJSONWriter):
    """A writer for QCSchema files."""

    def __init__(self, ccdata, *args, **kwargs):
        super(QCJSONWriter, self).__init__(ccdata, *args, **kwargs)

    def as_dict(self):
        qcschema_dict = dict()

        qcschema_dict["schema_name"] = "qcschema_output"
        qcschema_dict["schema_version"] = 1

        qcschema_dict["molecule"] = dict()
        qcschema_dict["molecule"]["schema_name"] = "qcschema_molecule"
        qcschema_dict["molecule"]["schema_version"] = 2

        # FIXME hack
        qcschema_dict["driver"] = "energy"

        metadata = self.ccdata.metadata

        qcschema_dict["success"] = metadata["success"]

        # FIXME
        # if "keywords" in metadata:
        #     qcschema_dict["keywords"] = set(metadata["keywords"])
        # else:
        qcschema_dict["keywords"] = dict()

        qcschema_dict["model"] = dict()
        # TODO methods should be an enum
        if metadata["methods"]:
            method = metadata["methods"][-1]
        else:
            # FIXME
            method = "HF"
        qcschema_dict["model"]["method"] = method

        # FIXME
        qcschema_dict["model"]["basis"] = ""

        qcschema_dict["provenance"] = dict()
        qcschema_dict["provenance"]["creator"] = metadata["package"]
        qcschema_dict["provenance"]["version"] = metadata["package_version"]
        # FIXME
        qcschema_dict["provenance"]["routine"] = ""

        qcschema_dict["properties"] = dict()
        qcschema_dict["properties"]["calcinfo_nbasis"] = self.ccdata.nbasis
        qcschema_dict["properties"]["calcinfo_nmo"] = self.ccdata.nmo
        # TODO nalpha and nbeta
        qcschema_dict["properties"]["calcinfo_natom"] = self.ccdata.natom

        qcschema_dict["properties"]["scf_total_energy"] = convertor(
            self.ccdata.scfenergies[-1], "eV", "hartree"
        )

        qcschema_dict["molecule"] = {
            "schema_name": "qcschema_molecule",
            "schema_version": 2,
            "geometry": self.ccdata.atomcoords[-1].flatten().tolist(),
            "symbols": [self.pt.element[atomno] for atomno in self.ccdata.atomnos],
        }

        if method == "HF":
            return_result = self.ccdata.scfenergies[-1]
        elif method == "CCSD":
            # FIXME for CCSD(T) etc.
            return_result = self.ccdata.ccenergies[-1]
        else:
            pass
        qcschema_dict["return_result"] = return_result

        return qcschema_dict

    def generate_repr(self):
        """Generate the QCSchema representation of the logfile data."""
        qcschema_dict = self.as_dict()
        if self.terse:
            return json.dumps(qcschema_dict, cls=NumpyAwareJSONEncoder)
        else:
            return json.dumps(qcschema_dict, cls=JSONIndentEncoder, sort_keys=True, indent=4)
