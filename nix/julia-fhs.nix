{ callPackage, buildFHSUserEnv, julia, cairo, fontconfig, freetype,
myjulia,
qt5, 
qt48Full,
xorg, libGL, zlib}:
let
in
  buildFHSUserEnv {
    name = "julia-fhs";
    targetPkgs = pkgs: with pkgs; [
      myjulia
      cairo fontconfig freetype
      qt5.full
      qt48Full
      libGL
      xorg.libX11 xorg.libXt
      xorg.libXext xorg.libXfixes xorg.libXi xorg.libXrandr xorg.libXrender
      xorg.libXtst xorg.libxcb xorg.xcbutilkeysyms zlib
    ];
    runScript = "julia main.jl";
  }
