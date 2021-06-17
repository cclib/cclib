# -*- coding: utf-8 -*-
#
# Copyright (c) 2021, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A writer for MolSSI quantum chemical JSON (QCSchema) files.
"""

import json

from cclib.io.cjsonwriter import CJSON as CJSONWriter
from cclib.io.cjsonwriter import JSONIndentEncoder, NumpyAwareJSONEncoder
from cclib.parser.utils import convertor

from qcschema import validate


class QCSchemaWriter(CJSONWriter):
    """A writer for QCSchema files."""

    def __init__(self, ccdata, *args, **kwargs):
        super().__init__(ccdata, *args, **kwargs)

    def as_dict(self):
        metadata = self.ccdata.metadata

        qcschema_dict = {
            "schema_name": "qcschema_output",
            "schema_version": 1,
            "molecule": {
                "schema_name": "qcschema_molecule",
                "schema_version": 2,
                "validated": True,
            },
        }

        # TODO This should be derived from a parsed job type.
        if hasattr(self.ccdata, "hessian"):
            driver = "hessian"
        elif hasattr(self.ccdata, "grads"):
            driver = "gradient"
        else:
            driver = "energy"
        # TODO what is "properties" for?
        qcschema_dict["driver"] = driver

        qcschema_dict["success"] = metadata["success"]

        # FIXME
        # if "keywords" in metadata:
        #     qcschema_dict["keywords"] = set(metadata["keywords"])
        # else:
        qcschema_dict["keywords"] = dict()

        # TODO methods should be an enum
        if metadata["methods"]:
            method = metadata["methods"][-1]
            if method == "DFT":
                method = metadata["functional"]
        else:
            # FIXME
            method = "HF"

        qcschema_dict["model"] = {"method": method, "basis": metadata["basis_set"]}

        qcschema_dict["provenance"] = {
            "creator": metadata["package"],
            "version": metadata["package_version"],
            # FIXME
            "routine": "",
        }

        qcschema_dict["molecule"] = {
            "schema_name": "qcschema_molecule",
            "schema_version": 2,
            "geometry": self.ccdata.atomcoords[-1].flatten().tolist(),
            "symbols": [self.pt.element[atomno] for atomno in self.ccdata.atomnos],
        }

        scf_total_energy = convertor(self.ccdata.scfenergies[-1], "eV", "hartree")
        mp2_correlation_energy = None
        mp2_total_energy = None
        ccsd_correlation_energy = None
        ccsd_total_energy = None

        if method == "HF":
            return_energy = scf_total_energy
        elif metadata["methods"][-1] == "DFT":
            # TODO parse XC energy
            return_energy = scf_total_energy
        elif method == "CCSD":
            mp2_total_energy = convertor(self.ccdata.mpenergies[-1][-1], "eV", "hartree")
            mp2_correlation_energy = convertor(
                self.ccdata.mpenergies[-1][-1] - self.ccdata.scfenergies[-1], "eV", "hartree"
            )
            ccsd_total_energy = convertor(self.ccdata.ccenergies[-1], "eV", "hartree")
            ccsd_correlation_energy = convertor(
                self.ccdata.ccenergies[-1] - self.ccdata.scfenergies[-1], "eV", "hartree"
            )
            return_energy = ccsd_total_energy
        else:
            raise RuntimeError("Don't know what to do with method {}".format(method))

        qcschema_dict["properties"] = {
            "calcinfo_nbasis": self.ccdata.nbasis,
            "calcinfo_nmo": self.ccdata.nmo,
            "calcinfo_natom": self.ccdata.natom,
            "calcinfo_nalpha": self.ccdata.homos[0] + 1,
            "calcinfo_nbeta": self.ccdata.homos[-1] + 1,
            "scf_iterations": self.ccdata.scfvalues[-1].shape[0],
            # TODO check units
            "scf_dipole_moment": self.ccdata.moments[1].tolist(),
            "scf_total_energy": scf_total_energy,
            "return_energy": return_energy,
            # TODO scf_one_electron energy
            # TODO scf_two_electron_energy
            # TODO nuclear_repulsion_energy
            # TODO mp2_same_spin_correlation_energy
            # TODO mp2_opposite_spin_correlation_energy
            # TODO mp2_singles_energy
            # TODO mp2_doubles_energy
            # TODO ccsd_same_spin_correlation_energy
            # TODO ccsd_opposite_spin_correlation_energy
            # TODO ccsd_singles_energy
            # TODO ccsd_doubles_energy
            # TODO ccsd_prt_pr_correlation_energy
            # TODO ccsd_prt_pr_total_energy
            # TODO ccsd_iterations
        }
        if mp2_correlation_energy is not None:
            qcschema_dict["properties"].update(
                {
                    "mp2_correlation_energy": mp2_correlation_energy,
                    "mp2_total_energy": mp2_total_energy,
                }
            )
        if ccsd_correlation_energy is not None:
            qcschema_dict["properties"].update(
                {
                    "ccsd_correlation_energy": ccsd_correlation_energy,
                    "ccsd_total_energy": ccsd_total_energy,
                }
            )

        if driver == "energy":
            return_result = return_energy
        elif driver == "gradient":
            return_result = self.ccdata.grads[-1].flatten().tolist()
        elif driver == "hessian":
            return_result = self.ccdata.hessian.flatten().tolist()
        elif driver == "properties":
            pass
        else:
            raise RuntimeError("Don't know driver {}".format(driver))
        qcschema_dict["return_result"] = return_result

        validate(qcschema_dict, schema_type="output")

        return qcschema_dict

    def generate_repr(self):
        """Generate the QCSchema representation of the logfile data."""
        qcschema_dict = self.as_dict()
        if self.terse:
            return json.dumps(qcschema_dict, cls=NumpyAwareJSONEncoder)
        else:
            return json.dumps(qcschema_dict, cls=JSONIndentEncoder, sort_keys=True, indent=4)
