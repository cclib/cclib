final: prev: {
  python3 = prev.python3.override (old: {
    packageOverrides = prev.lib.composeExtensions (old.packageOverrides or (_: _: { })) (pfinal: pprev: {
      cclib = pfinal.callPackage ./cclib.nix {
        inherit (final.qchem.python3.pkgs) psi4 pyquante;
      };

      # Fixes
      # 2 tests currently fail with database lookups
      biopython = pprev.biopython.overrideAttrs (old: {
        doCheck = false;
        doInstallCheck = false;
      });
    });
  });
}
