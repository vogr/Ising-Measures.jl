let
  pkgs = import <nixpkgs> {};
  /* myjulia = pkgs.callPackage ./nix/julia/0.6.nix {
     gmp = pkgs.gmp6;
    openblas = pkgs.openblasCompat;
    inherit (pkgs.darwin.apple_sdk.frameworks) CoreServices ApplicationServices;
    llvm = pkgs.llvm_39;
    }; */
  myjulia = pkgs.julia;
  julia-fhs = pkgs.callPackage ./nix/julia-fhs.nix { myjulia=myjulia; };
  stdenv = pkgs.stdenv;
in
  stdenv.mkDerivation {
    name = "env";
    buildInputs = [ julia-fhs ];
    shellHook="julia-fhs;exit";
  }
