{ pkgs ? import <nixpkgs> {} }:
with pkgs;


let


gmp-dev-meta = (pkgs.gmp.meta // {outputsToInstall =["out" "dev"]; });

in


{

inherit pkgs;


gmp-dev = pkgs.gmp // {meta = gmp-dev-meta;};

# Run :
# source /applis/site/nix.sh
# nix-env -f . -iA gmp-dev

}

