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
from cclib.parser.utils import convertor, find_package

_found_qcschema = find_package("qcschema")
if _found_qcschema:
    import qcschema


class QCSchemaWriter(CJSONWriter):
    """A writer for QCSchema files."""

    def __init__(self, ccdata, *args, **kwargs):
        super().__init__(ccdata, *args, **kwargs)

    def as_dict(self, validate=True):
        metadata = self.ccdata.metadata

        qcschema_dict = {
            "schema_name": "qcschema_output",
            "schema_version": 1,
            "molecule": {
                "geometry": self.ccdata.atomcoords[-1].flatten().tolist(),
                "molecular_charge": self.ccdata.charge,
                "molecular_multiplicity": self.ccdata.mult,
                "schema_name": "qcschema_molecule",
                "schema_version": 2,
                "symbols": [self.pt.element[atomno] for atomno in self.ccdata.atomnos],
                "validated": True,
            },
            "provenance": {
                "creator": metadata["package"],
                "version": metadata["package_version"],
                "routine": f"{self.__module__}.{self.__class__.__qualname__}",
            },
            "success": metadata["success"],
            "error": None,
            "stdout": None,
            "stderr": None,
        }

        # TODO This should be derived from a parsed job type.  It's also not
        # quite right when considering that many jobs (for example, geometry
        # optimizations) are actually procedures and not drivers.
        if hasattr(self.ccdata, "hessian"):
            driver = "hessian"
        elif hasattr(self.ccdata, "grads"):
            driver = "gradient"
        else:
            driver = "energy"
        # TODO what is "properties" for?
        qcschema_dict["driver"] = driver

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

        qcschema_dict["model"] = {"method": method.lower(), "basis": metadata["basis_set"].lower()}

        scf_total_energy = convertor(self.ccdata.scfenergies[-1], "eV", "hartree")
        mp2_correlation_energy = None
        mp2_total_energy = None
        ccsd_correlation_energy = None
        ccsd_total_energy = None

        if method == "HF":
            return_energy = scf_total_energy
        elif metadata["methods"][-1] == "DFT":
            return_energy = scf_total_energy
        elif method == "CCSD":
            if hasattr(self.ccdata, "mpenergies"):
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

        scf_dipole_moment = None
        if hasattr(self.ccdata, "moments"):
            # FIXME We currently have no easy way to tell if the moments from
            # a correlated calculation are true correlated (relaxed) density
            # moments, or just SCF density moments.
            dipole_moment = convertor(self.ccdata.moments[1], "Debye", "ebohr")
            scf_dipole_moment = dipole_moment.tolist()

        qcschema_dict["properties"] = {
            "calcinfo_nalpha": self.ccdata.homos[0] + 1,
            "calcinfo_natom": self.ccdata.natom,
            "calcinfo_nbasis": self.ccdata.nbasis,
            "calcinfo_nbeta": self.ccdata.homos[-1] + 1,
            "calcinfo_nmo": self.ccdata.nmo,
            "return_energy": return_energy,
            "scf_iterations": self.ccdata.scfvalues[-1].shape[0],
            "scf_total_energy": scf_total_energy,
            # TODO These properties aren't parsed yet.
            # ccsd_dipole_moment
            # ccsd_doubles_energy
            # ccsd_iterations
            # ccsd_opposite_spin_correlation_energy
            # ccsd_prt_pr_correlation_energy
            # ccsd_prt_pr_dipole_moment
            # ccsd_prt_pr_total_energy
            # ccsd_same_spin_correlation_energy
            # ccsd_singles_energy
            # ccsdt_correlation_energy
            # ccsdt_dipole_moment
            # ccsdt_iterations
            # ccsdt_total_energy
            # ccsdtq_correlation_energy
            # ccsdtq_dipole_moment
            # ccsdtq_iterations
            # ccsdtq_total_energy
            # mp2_dipole_moment
            # mp2_doubles_energy
            # mp2_opposite_spin_correlation_energy
            # mp2_same_spin_correlation_energy
            # mp2_singles_energy
            # nuclear_repulsion_energy
            # scf_one_electron_energy
            # scf_two_electron_energy
            # scf_vv10_energy
            # scf_xc_energy
        }
        if hasattr(self.ccdata, "dispersionenergies"):
            qcschema_dict["properties"][
                "scf_dispersion_correction_energy"
            ] = self.ccdata.dispersionenergies[-1]
        if scf_dipole_moment is not None:
            qcschema_dict["scf_dipole_moment"] = scf_dipole_moment
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
        else:
            raise RuntimeError("Don't know driver {}".format(driver))
        qcschema_dict["return_result"] = return_result

        if validate:
            qcschema.validate(qcschema_dict, schema_type="output")

        return qcschema_dict

    def generate_repr(self):
        """Generate the QCSchema representation of the logfile data."""
        qcschema_dict = self.as_dict()
        if self.terse:
            return json.dumps(qcschema_dict, cls=NumpyAwareJSONEncoder)
        else:
            return json.dumps(qcschema_dict, cls=JSONIndentEncoder, sort_keys=True, indent=4)
