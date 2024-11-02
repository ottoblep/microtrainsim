{
  description = "microtrainsim flake";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.05";
    nix-matlab = {
      url = "gitlab:doronbehar/nix-matlab";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs = all@{ self, nixpkgs, nix-matlab, ... }:
    let
      pkgs = import nixpkgs {
        system = "x86_64-linux";
        config.allowUnfree = true;
      };
    in
    {
      # Utilized by `nix develop`
      devShells.x86_64-linux = {
        default = pkgs.mkShell {
          buildInputs = (with nix-matlab.packages.x86_64-linux; [
            matlab
            matlab-mlint
            matlab-mex
          ]);

          shellHook = nix-matlab.shellHooksCommon;

          packages = with pkgs; with pkgs.python311Packages; [
            (python3.withPackages (python-pkgs: [
              pandas
              numpy
              scipy
              seaborn
              plotly
              networkx
              shapely
              geopy
            ]))
          ];
        };

        cpp = pkgs.mkShell {
            packages = with pkgs; [
              # C++ Compiler is already part of stdenv
              boost
              catch2
              cmake
            ];
          };
      };
    };
}
