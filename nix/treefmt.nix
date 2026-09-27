{
  lib,
  pkgs,
  pythonEnv,
  ...
}:
let
  mypy-check = pkgs.writeShellApplication {
    name = "mypy-check";
    runtimeInputs = [ pythonEnv ];
    text = ''
      mypy cloningkit tests
    '';
  };
in
{
  projectRootFile = "flake.nix";

  programs.deadnix.enable = true;
  programs.nixfmt.enable = true;

  programs.mdformat.enable = true;

  programs.shellcheck.enable = true;
  programs.shfmt.enable = true;

  programs.taplo.enable = true;

  programs.ruff-check.enable = true;
  programs.ruff-format.enable = true;

  settings.formatter.deadnix.priority = 1;
  settings.formatter.nixfmt.priority = 2;

  settings.formatter.shellcheck.priority = 1;
  settings.formatter.shfmt.priority = 2;

  settings.formatter.ruff-check.priority = 1;
  settings.formatter.ruff-format.priority = 2;

  # Runs after ruff so it never reports against code ruff is about to rewrite.
  settings.formatter.mypy-check = {
    command = lib.getExe mypy-check;
    includes = [ "*.py" ];
    priority = 3;
  };
}
