{
  description = "CloningKit: Cloning experiment design library for wet-lab work";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    bioinformatics-toolkits = {
      url = "github:SBEE-Lab/bioinformatics-toolkits";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    treefmt-nix = {
      url = "github:numtide/treefmt-nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs =
    inputs@{
      self,
      nixpkgs,
      ...
    }:
    let
      inherit (nixpkgs) lib;

      systems = [
        "x86_64-linux"
        "aarch64-linux"
        "aarch64-darwin"
      ];
      eachSystem = lib.genAttrs systems;

      overlay = final: prev: {
        pythonPackagesExtensions = prev.pythonPackagesExtensions ++ [
          (_pyfinal: _pyprev: { inherit (final.bioinformatics-toolkits) pydna; })
        ];
      };

      pkgsFor = eachSystem (
        system:
        import nixpkgs {
          inherit system;
          overlays = [
            inputs.bioinformatics-toolkits.overlays.shared-nixpkgs
            overlay
          ];
        }
      );

      pythonEnv = eachSystem (
        system:
        pkgsFor.${system}.python3.withPackages (
          ps: with ps; [
            mypy
            primer3
            pydna
            python-codon-tables
            pytest
          ]
        )
      );

      devShells = eachSystem (
        system:
        let
          pkgs = pkgsFor.${system};
        in
        {
          default = pkgs.mkShellNoCC {
            packages = [
              pkgs.ruff
              pkgs.bash
              pythonEnv.${system}
              formatter.${system}
            ];
          };
        }
      );

      treefmtEval = eachSystem (
        system:
        inputs.treefmt-nix.lib.evalModule pkgsFor.${system} (
          import ./nix/treefmt.nix {
            inherit lib;
            pkgs = pkgsFor.${system};
            pythonEnv = pythonEnv.${system};
          }
        )
      );

      formatter = eachSystem (system: treefmtEval.${system}.config.build.wrapper);

    in
    {
      inherit devShells formatter;
      overlays.default = overlay;
      checks = eachSystem (system: {
        formatting = treefmtEval.${system}.config.build.check self;
        devshell-default = devShells.${system}.default;
      });
    };
}
