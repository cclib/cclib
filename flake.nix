{
  description = "Parsers and algorithms for computational chemistry logfiles";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";

    qchem.url = "github:nix-qchem/nixos-qchem";

    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, qchem, flake-utils }:
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
              ] ++ p.cclib.nativeBuildInputs
              ++ p.cclib.propagatedBuildInputs
              );
            in
            pkgs.mkShell {
              buildInputs = with pkgs; [
                pyEnv
                black
              ];
            };

          formatter = pkgs.nixpkgs-fmt;
        }) // {
      overlays.default = overlay;
    };
}
