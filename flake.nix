{
  description = "DSA110 MultiBeam Heimdall nix flake";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";

    psrdada = {
      url = "github:kiranshila/psrdada.nix";
      inputs = {
        nixpkgs.follows = "nixpkgs";
        flake-utils.follows = "flake-utils";
      };
    };
  };

  outputs = {
    nixpkgs,
    flake-utils,
    psrdada,
    ...
  }:
    flake-utils.lib.eachDefaultSystem (
      system: let
        pkgs = import nixpkgs {
          inherit system;
          config.allowUnfree = true;
        };
        cudaStdenv = pkgs.cudaPackages.backendStdenv;
        nativeBuildInputs = with pkgs; [
          cmake
        ];
        buildInputs = with pkgs; [
          # Dependencies
          cudaPackages.cudatoolkit
          psrdada.packages.${system}.default
          boost
          # Linting
          alejandra
        ];
        heimdall = cudaStdenv.mkDerivation rec {
          pname = "heimdall";
          src = ./.;
          version = "0.1";
          inherit nativeBuildInputs buildInputs;
        };
      in rec {
        defaultPackage = heimdall;
      }
    );
}
