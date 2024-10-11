{
  description = "Nix Dev Shells";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.05";
  outputs = all@{ self, nixpkgs, ... }:
    let
      pkgs = import nixpkgs {
        system = "x86_64-linux";
        config.allowUnfree = true;
      };
    in
    {
      # Utilized by `nix develop`
      devShell.x86_64-linux =
        pkgs.mkShell {
          packages = with pkgs; with pkgs.python311Packages; [
            octaveFull
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
    };
}
