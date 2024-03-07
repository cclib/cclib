{
  description = "Parsers and algorithms for computational chemistry logfiles";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";

    qchem.url = "github:nix-qchem/nixos-qchem";

    flake-utils.url = "github:numtide/flake-utils";

    cachixPreCommit.url = "github:cachix/pre-commit-hooks.nix";
  };

  nixConfig = {
    bash-prompt = "\\[\\e[0;1;38;5;214m\\]cclib\\[\\e[0;1m\\]:\\[\\e[0;1;38;5;39m\\]\\w\\[\\e[0;1m\\]$ \\[\\e[0m\\]";
  };

  outputs = { self, nixpkgs, qchem, flake-utils, cachixPreCommit }:
    let
      overlay = import ./nix/overlay.nix;
    in
    flake-utils.lib.eachDefaultSystem
      (system:
        let
          pkgs = import nixpkgs {
            inherit system;
            overlays = [
              qchem.overlays.default
              overlay
            ];
          };
        in
        {
          packages.default = pkgs.python3.pkgs.cclib;

          devShells.default =
            let
              pyEnv = pkgs.python3.withPackages (p: with p; [
                pip
                pytest-cov
                coverage
                sphinx
                sphinx-rtd-theme
                identify
                cfgv
                pre-commit-hooks
              ] ++ p.cclib.nativeBuildInputs
              ++ p.cclib.propagatedBuildInputs
              );
            in
            pkgs.mkShell {
              buildInputs = with pkgs; [
                pyEnv
                black
                pre-commit
                ruff
                actionlint
                isort
              ];

              #inherit (self.checks.${system}.pre-commit-check) shellHook;
            };

          formatter = pkgs.nixpkgs-fmt;

          checks.pre-commit-check = cachixPreCommit.lib.${system}.run {
            src = ./.;

            hooks =
              let
                excludes = [
                  "\.bpa"
                  "\.bpaspin"
                  "\.cube"
                  "\.in"
                  "\.inp"
                  "\.molden"
                  "\.out"
                  "\.h5"
                  "\.png"
                  "^data/ADF"
                  "^data/DALTON"
                  "^data/FChk"
                  "^data/GAMESS"
                  "^data/Gaussian"
                  "^data/Jaguar"
                  "^data/Molcas"
                  "^data/Molpro"
                  "^data/MOPAC"
                  "^data/NBO"
                  "^data/NWChem"
                  "^data/ORCA"
                  "^data/Psi4"
                  "^data/QChem"
                  "^data/regression"
                  "^data/Turbomole"
                  "^data/XTB"
                ];
                stdHooks = [
                  "trailing-whitespace-fixer"
                  "end-of-file-fixer"
                  "fix-byte-order-marker"
                  "mixed-line-ending"
                  "check-merge-conflict"
                  "check-added-large-files"
                ];
              in
              {
                nixpkgs-fmt.enable = true;
                actionlint.enable = true;
                isort.enable = true;
                ruff.enable = true;
                check-json = {
                  enable = true;
                  entry = "${pkgs.python3.pkgs.pre-commit-hooks}/bin/check-json";
                  files = "\\.json";
                };
                check-yaml = {
                  enable = true;
                  entry = "${pkgs.python3.pkgs.pre-commit-hooks}/bin/check-yaml";
                  files = "\\.yml";
                };
                check-toml = {
                  enable = true;
                  entry = "${pkgs.python3.pkgs.pre-commit-hooks}/bin/check-toml";
                  files = "\\.toml";
                };
              } // builtins.listToAttrs (builtins.map
                (hook: {
                  name = hook;
                  value = {
                    enable = true;
                    entry = "${pkgs.python3.pkgs.pre-commit-hooks}/bin/${hook}";
                    inherit excludes;
                  };
                })
                stdHooks);
          };
        }) // {
      overlays.default = overlay;
    };
}
