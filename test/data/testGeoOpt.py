# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test geometry optimization logfiles in cclib"""

import numpy
from common import get_minimum_carbon_separation
from skip import skipForLogfile, skipForParser


class GenericGeoOptTest:
    """Generic geometry optimization unittest"""

    # In STO-3G, H has 1, C has 3.
    nbasisdict = {1: 1, 6: 5}

    # Some programs print surplus atom coordinates by default.
    extracoords = 0

    # Some programs do surplus SCF cycles by default.
    extrascfs = 0

    # Approximate B3LYP energy of dvb after SCF in STO-3G.
    b3lyp_energy = -380.90674109218116
    b3lyp_tolerance = 1.4699729516340807

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("MOPAC", "The success status is not parsed yet")
    def test_success(self, data) -> None:
        assert data.metadata["success"]

    def testnatom(self, data) -> None:
        """Is the number of atoms equal to 20?"""
        assert data.natom == 20

    def testatomnos(self, data) -> None:
        """Are the atomnos correct?"""
        # This will work only for numpy
        # self.assertEqual(data.atomnos.dtype.char, 'i')

        atomnos_types = [numpy.issubdtype(atomno, numpy.signedinteger) for atomno in data.atomnos]
        assert numpy.all(atomnos_types)

        assert data.atomnos.shape == (20,)

        count_C = sum(data.atomnos == 6)
        count_H = sum(data.atomnos == 1)
        assert count_C + count_H == 20

    def testatomcoords(self, data) -> None:
        """Are atomcoords consistent with natom and Angstroms?"""
        natom = len(data.atomcoords[0])
        ref = data.natom
        msg = f"natom is {int(ref)} but len(atomcoords[0]) is {int(natom)}"
        assert natom == ref, msg

    def testatomcoords_units(self, data) -> None:
        """Are atomcoords consistent with Angstroms?"""
        min_carbon_dist = get_minimum_carbon_separation(data)
        dev = abs(min_carbon_dist - 1.34)
        assert dev < 0.15, f"Minimum carbon dist is {min_carbon_dist:.2f} (not 1.34)"

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("xTB", "Not implemented yet")
    def testcharge_and_mult(self, data) -> None:
        """Are the charge and multiplicity correct?"""
        assert data.charge == 0
        assert data.mult == 1

    @skipForParser("MOPAC", "Not implemented.")
    @skipForParser("xTB", "not implemented yet")
    def testnbasis(self, data) -> None:
        """Is the number of basis set functions correct?"""
        count = sum([self.nbasisdict[n] for n in data.atomnos])
        assert data.nbasis == count

    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testcoreelectrons(self, data) -> None:
        """Are the coreelectrons all 0?"""
        ans = numpy.zeros(data.natom, "i")
        numpy.testing.assert_array_equal(data.coreelectrons, ans)

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("MOPAC", "Not implemented.")
    @skipForParser("xTB", "Not implemented yet")
    def testhomos(self, data) -> None:
        """Is the index of the HOMO equal to 34?"""
        ref = numpy.array([34], "i")
        msg = f"{numpy.array_repr(data.homos)} != array([34], 'i')"
        numpy.testing.assert_array_equal(data.homos, ref, msg)

    @skipForParser("MOPAC", "The scfvalues attribute is not parsed yet")
    @skipForParser("xTB", "not implemented yet")
    def testscfvaluetype(self, data) -> None:
        """Are scfvalues and its elements the right type?"""
        assert isinstance(data.scfvalues, list)
        assert isinstance(data.scfvalues[0], numpy.ndarray)

    def testscfenergy(self, data) -> None:
        """Is the SCF energy close to target?"""
        scf = data.scfenergies[-1]
        ref = self.b3lyp_energy
        tol = self.b3lyp_tolerance
        msg = f"Final SCF energy: {scf:f} not {int(ref)} +- {int(tol)}eV"
        assert abs(scf - ref) < 40, msg

    @skipForParser("xTB", "Not implemented yet")
    def testscfenergydim(self, data) -> None:
        """Is the number of SCF energies consistent with atomcoords?"""
        count_scfenergies = data.scfenergies.shape[0] - self.extrascfs
        count_atomcoords = data.atomcoords.shape[0] - self.extracoords
        assert count_scfenergies == count_atomcoords

    @skipForParser("MOPAC", "The scftargets attribute is not parsed yet")
    @skipForParser("xTB", "not implemented yet")
    def testscftargetdim(self, data) -> None:
        """Do the scf targets have the right dimensions?"""
        dim_scftargets = data.scftargets.shape
        dim_scfvalues = (len(data.scfvalues), len(data.scfvalues[0][0]))
        assert dim_scftargets == dim_scfvalues

    @skipForParser("MOPAC", "Not implemented.")
    @skipForParser("xTB", "not implemented yet")
    def testgeovalues_atomcoords(self, data) -> None:
        """Are atomcoords consistent with geovalues?"""
        count_geovalues = len(data.geovalues)
        count_coords = len(data.atomcoords) - self.extracoords
        msg = f"len(atomcoords) is {int(count_coords)} but len(geovalues) is {int(count_geovalues)}"
        assert count_geovalues == count_coords, msg

    @skipForParser("MOPAC", "Not implemented.")
    @skipForParser("xTB", "not implemented yet")
    def testgeovalues_scfvalues(self, data) -> None:
        """Are scfvalues consistent with geovalues?"""
        count_scfvalues = len(data.scfvalues) - self.extrascfs
        count_geovalues = len(data.geovalues)
        assert count_scfvalues == count_geovalues

    @skipForParser("MOPAC", "Not implemented.")
    @skipForParser("xTB", "not implemented yet")
    def testgeotargets(self, data) -> None:
        """Do the geo targets have the right dimensions?"""
        dim_geotargets = data.geotargets.shape
        dim_geovalues = (len(data.geovalues[0]),)
        assert dim_geotargets == dim_geovalues

    @skipForParser("MOPAC", "Not implemented.")
    @skipForParser("xTB", "not implemented yet")
    def testoptdone(self, data) -> None:
        """Has the geometry converged and set optdone to True?"""
        assert isinstance(data.optdone, list)
        assert len(data.optdone) == 1
        assert all(isinstance(val, int) for val in data.optdone)

    @skipForParser("MOPAC", "Not implemented.")
    @skipForParser("xTB", "not implemented yet")
    def testgeoconverged(self, data) -> None:
        """Has the geometry converged and set optdone to True?"""
        assert numpy.all(numpy.abs(data.geovalues[-1]) <= data.geotargets)

    @skipForParser("ADF", "Not implemented.")
    @skipForParser("DALTON", "Not implemented.")
    @skipForParser("GAMESS", "Not implemented.")
    @skipForParser("GAMESSUK", "Not implemented.")
    @skipForParser("Jaguar", "Not implemented.")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Molpro", "Not implemented.")
    @skipForParser("MOPAC", "Not implemented.")
    @skipForParser("NWChem", "Not implemented.")
    @skipForParser("ORCA", "Not implemented.")
    @skipForParser("QChem", "Not implemented.")
    @skipForParser("xTB", "not implemented yet")
    def testoptstatus(self, data) -> None:
        """Is optstatus consistent with geovalues and reasonable?"""
        assert len(data.optstatus) == len(data.geovalues)
        assert data.optstatus[0] == data.OPT_NEW
        for i in range(1, len(data.optstatus) - 1):
            assert data.optstatus[i] == data.OPT_UNKNOWN
        assert data.optstatus[-1] == data.OPT_DONE

    @skipForParser("ADF", "Not implemented yet")
    @skipForParser("DALTON", "Not implemented yet")
    @skipForParser("FChk", "Rotational constants are never written to fchk files")
    @skipForParser("GAMESS", "Not implemented yet")
    @skipForParser("GAMESSUK", "Not implemented yet")
    @skipForParser("Jaguar", "Not implemented yet")
    @skipForParser("Molcas", "Not implemented yet")
    @skipForParser("Molpro", "Not implemented yet")
    @skipForLogfile("MOPAC/basicMOPAC2016", "Not present in this file")
    @skipForParser("NWChem", "Not implemented yet")
    @skipForParser("ORCA", "Not implemented yet")
    @skipForParser("Psi4", "Not implemented yet")
    @skipForParser("QChem", "Not implemented yet")
    @skipForParser("Turbomole", "Not implemented yet")
    @skipForParser("xTB", "not implemented yet")
    def testrotconsts(self, data) -> None:
        """Each geometry leads to a row in the rotational constants entry."""
        assert data.rotconsts.shape == (len(data.atomcoords), 3)

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    def testmoenergies(self, data) -> None:
        """Are only the final MOs parsed?"""
        assert len(data.moenergies) == 1
        if hasattr(data, "mocoeffs"):
            assert len(data.mocoeffs) == 1

    @skipForParser("ADF", "Not implemented.")
    @skipForParser("DALTON", "Not implemented.")
    @skipForParser("GAMESS", "Not implemented.")
    @skipForParser("GAMESSUK", "Not implemented.")
    @skipForParser("Jaguar", "Not implemented.")
    @skipForParser("MOPAC", "Not implemented.")
    @skipForParser("NWChem", "Not implemented.")
    @skipForParser("xTB", "not implemented yet")
    def testgradsdim(self, data) -> None:
        """Do the grads have the right dimensions?"""
        assert data.grads.shape == (len(data.geovalues), data.natom, 3)

    @skipForLogfile("Molcas/basicOpenMolcas18.0/dvb_gopt.out", "not implemented yet")
    def testmetadata_success(self, data) -> None:
        """Does metadata have expected keys and values?"""
        assert "success" in data.metadata
        assert data.metadata["success"]


class ADFGeoOptTest(GenericGeoOptTest):
    """Customized geometry optimization unittest"""

    extracoords = 1
    extrascfs = 1

    b3lyp_energy = -5.144905330719283
    b3lyp_tolerance = 0.03674932379085202


class DALTONGeoOptTest(GenericGeoOptTest):
    """Customzed geometry optimization unittest"""

    # DALTON will normally print the geometry several extra times as the "final geometry"
    # when an optimziation converges. We don't parse those coordinates, but the parser
    # does catch the geometry printed in the final static property calculation when
    # that is done for the final geometry (presumably always).
    extracoords = 1

    # Although DALTON generally has three criteria for convergence, it normally only
    # requires two of them to end a geometry optimization. This is printed in the output
    # and can probably be tweaked in the input, but we don't parsed that in cclib.
    def testgeoconverged(self, data) -> None:
        """Has the geometry converged and set optdone to True?"""
        convergence = numpy.abs(data.geovalues[-1]) <= data.geotargets
        assert sum(convergence) >= 2


class GaussianGeoOptTest(GenericGeoOptTest):
    """Customized geometry optimization unittest"""

    def testgradsorientation(self, data) -> None:
        """Are the orientations for grads and atomcoords are same?"""
        # since z-coordinates of atomcoords are all 0 for dvb, z-values of grads should be all 0
        assert numpy.all(numpy.abs(data.grads[:, :, 2]) < 1e-14)


class MolcasGeoOptTest(GenericGeoOptTest):
    """Customized geometry optimization unittest"""

    # Molcas prints the input coordinates and performs the &scf job
    # once before entering the optimization part where the coordinates and
    # scf section are printed for each iteration. Hence we have an extra set of
    # coordinates and extra set of SCF attributes (scfenergies, scftargets & scfvalues).
    extracoords = 1
    extrascfs = 1


class MolproGeoOptTest(GenericGeoOptTest):
    """Customized geometry optimization unittest"""

    # Note that these extra coordinates and energies will be available only
    # if the appropriate output is parsed, and Molpro often saves the initial
    # SCF run and subsequent geometry optimization to separate files, which
    # both need to be given to the cclib parser (as a list).
    extracoords = 1
    extrascfs = 2

    # Here is what the manual has to say about convergence:
    # The standard MOLPRO convergency criterion requires the maximum component of the gradient
    # to be less then $3 \cdot 10^{-4}$ [a.u.] and the maximum energy change to be less than
    # $1 \cdot 10^{-6}$ [H] or the maximum component of the gradient to be less then
    # $3 \cdot 10^{-4}$ [a.u.] and the maximum component of the step to be less then
    # $3 \cdot 10^{-4}$ [a.u.].
    #
    # It is also possible to use the convergency criterion of (...)
    #
    # Source: https://www.molpro.net/info/2012.1/doc/manual/node592.html
    def testgeoconverged(self, data) -> None:
        """Has the geometry converged and set optdone to True?"""
        target_e, target_g, target_s = data.geotargets
        value_e, value_g, value_s = data.geovalues[-1]
        converged = (value_e < target_e and value_g < target_g) or (
            value_g < target_g and value_s < target_s
        )
        assert converged


class MOPACGeoOptTest(GenericGeoOptTest):
    """Customized geometry optimization unittest for MOPAC."""

    # The geometry optimization unit test logfile uses a PM7 Hamiltonian.
    b3lyp_energy = 0.08166294233283113
    b3lyp_tolerance = 1.6e-5


class NWChemGeoOptTest(GenericGeoOptTest):
    """Customized restricted single point HF unittest"""

    # NWChem typically prints the coordinates in the input module, at the
    # beginning of each geometry optimization step, and then again after
    # the optimziation is finished, so the first and last coordinates
    # are repeated. On the other hand, each optimization step often
    # involves a line search which we don't parse (see parse code for details).
    extracoords = 2
    extrascfs = 0


class OrcaGeoOptTest(GenericGeoOptTest):
    """Customized geometry optimization unittest"""

    extracoords = 1
    extrascfs = 1

    # Besides all the geovalues being below their tolerances, ORCA also considers
    # an optimization finished in some extra cases. These are:
    #   1) everything converged except the energy (within 25 x tolerance)
    #   2) gradient is overachieved and displacement is reasonable (3 x tolerance)
    #   3) displacement is overachieved and gradient is reasonable (3 x tolerance)
    #   4) energy, gradients and angles are converged (displacements not considered)
    # All these exceptions are signaleld in the output with some comments, and here
    # we include the first three exceptions for the pruposes of the unit test.
    def testgeoconverged(self, data) -> None:
        """Has the geometry converged and set optdone to True?"""

        targets = data.geotargets
        values = numpy.abs(data.geovalues[-1])

        target_e = targets[0]
        target_g = targets[1:3]
        target_x = targets[3:]
        value_e = values[0]
        value_g = values[1:3]
        value_x = values[3:]

        conv_all = all(values < targets)
        conv_e = value_e < 25 * target_e and all(value_g < target_g) and all(value_x < target_x)
        conv_g = (
            value_e < target_e and all(value_g < target_g / 3.0) and all(value_x < target_x * 3.0)
        )
        conv_x = (
            value_e < target_e and all(value_g < target_g * 3.0) and all(value_x < target_x / 3.0)
        )
        converged = conv_all or conv_e or conv_g or conv_x
        assert converged


class Psi4GeoOptTest(GenericGeoOptTest):
    """Customized geometry optimization unittest"""

    # Psi has a number of different convergence strategies to choose from, as described here:
    #     http://sirius.chem.vt.edu/psi4manual/latest/optking.html
    # and the default is to check that the max. force is converged and if the max energy change
    # or dispalcement is converged. This is in fact what is tested below.
    def testgeoconverged(self, data) -> None:
        """Has the geometry converged and set optdone to True?"""

        targets = data.geotargets
        values = numpy.abs(data.geovalues[-1])

        # Since the other criteria are not used and are not printed in this case, they should
        # be parsed as numpy.inf, for which we can check.
        assert numpy.isinf(targets[2])
        assert numpy.isinf(targets[4])

        conv = values[1] < targets[1] and (values[0] < targets[0] or values[3] < targets[3])
        assert conv


class XTBGeoOptTest(GenericGeoOptTest):
    """Customized restricted single point unittest"""

    def testscfenergy(self, data) -> None:
        """Is the SCF energy within the target?"""
        assert abs(data.scfenergies[-1] - -26.438242468348) < 1.0e-6


class TurbomoleKeepGeoOptTest(GenericGeoOptTest):
    """Customized geometry optimization unittest"""

    # In Turbomole, each optimisation step is written to its own file,
    # (job.1, job.2 ... job.last) and consists of three (or more)
    # submodule steps:
    # - grad: Calculation of gradient.
    # - statpt: Update of coordinates.
    # - dscf: Calculation of SCF from new geometry.
    # - [rimp2/ricc2/ccsdf12]: Calculation of higher order energies.
    #
    # In addition, this cycle is started by an initial dscf step stored
    # in the job.0 file. Each submodule will print the current atom coords
    # during startup, so there will always be one more atom coords than
    # there are geom steps (we make sure in the parser not to parse
    # identical coords more than once).
    #
    # However, the number of SCF energies will depend on whether the
    # jobex script was called with the -keep flag or not. Without -keep
    # (the default), jobex will delete all intermediate files (including
    # job.0), so we can only parse job.final. In this case, only one
    # SCF energy will be available (the final energy). Alternatively,
    # with -keep we can parse the energy at each step, including the
    # initial energy, and so len(scfenergies) == len(atomcoords) (one
    # greater than the number of opt steps).
    #
    # The test data was called with jobex -keep.
    extracoords = 1
    extrascfs = 1


class TurbomoleGeoOptTest(GenericGeoOptTest):
    """Customized geometry optimization unittest"""

    # The test data was not called with jobex -keep.
    extracoords = 1
    extrascfs = 0

    def testoptstatus(self, data) -> None:
        """Is optstatus consistent with geovalues and reasonable?"""
        assert len(data.optstatus) == len(data.geovalues)
        # We only have the final energy available, so there's no point looking for OPT_NEW.
        for i in range(1, len(data.optstatus) - 1):
            assert data.optstatus[i] == data.OPT_UNKNOWN
        assert data.optstatus[-1] == data.OPT_DONE
