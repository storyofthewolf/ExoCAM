#!/bin/sh
# install_sourcemods.sh -- assemble an ExoCAM config's SourceMods into a case.
#
# Usage:  install_sourcemods.sh CONFIG DEST
#   CONFIG  a config under cesm1.2.1/configs/ that has a COMPONENTS file
#           (cam_aqua_fv, cam_land_fv, cam_mixed_fv)
#   DEST    the case's SourceMods directory, e.g. $CASEROOT/SourceMods
#
# Two layers, copied in this order:
#   1. cesm1.2.1/SourceMods/src.<comp>/       shared by every FV config
#   2. cesm1.2.1/configs/CONFIG/SourceMods/   the config's own files
# Only the components listed in configs/CONFIG/COMPONENTS are copied, so a
# case receives exactly the files it would have had from the per-config
# layout. A file present in both layers is an error, not an override: each
# file has exactly one home.
#
# Existing files in DEST are overwritten; nothing is deleted.

set -eu

if [ $# -ne 2 ]; then
  echo "usage: $0 CONFIG DEST" >&2
  exit 2
fi
config=$1
dest=$2

root=$(cd "$(dirname "$0")" && pwd)
shared=$root/SourceMods
cfgdir=$root/configs/$config
complist=$cfgdir/COMPONENTS

if [ ! -f "$complist" ]; then
  echo "ERROR: $complist not found; '$config' is not a shared-tree config." >&2
  exit 1
fi
if [ ! -d "$dest" ]; then
  echo "ERROR: destination $dest does not exist (expected \$CASEROOT/SourceMods)." >&2
  exit 1
fi

comps=$(grep -v '^[[:space:]]*#' "$complist" | tr -s ' \t\n' ' ')

# Refuse a file that lives in both layers.
for comp in $comps; do
  [ -d "$cfgdir/SourceMods/src.$comp" ] || continue
  for f in "$cfgdir/SourceMods/src.$comp"/*; do
    [ -e "$f" ] || continue
    if [ -e "$shared/src.$comp/$(basename "$f")" ]; then
      echo "ERROR: src.$comp/$(basename "$f") is in both the shared tree and $config." >&2
      exit 1
    fi
  done
done

# A config overlay may not carry a component its COMPONENTS file omits.
if [ -d "$cfgdir/SourceMods" ]; then
  for d in "$cfgdir/SourceMods"/src.*; do
    [ -d "$d" ] || continue
    c=${d##*/src.}
    case " $comps " in
      *" $c "*) ;;
      *) echo "ERROR: $config has SourceMods/src.$c but COMPONENTS does not list '$c'." >&2; exit 1 ;;
    esac
  done
fi

n=0
for comp in $comps; do
  mkdir -p "$dest/src.$comp"
  for layer in "$shared/src.$comp" "$cfgdir/SourceMods/src.$comp"; do
    [ -d "$layer" ] || continue
    for f in "$layer"/*; do
      [ -f "$f" ] || continue
      cp -p "$f" "$dest/src.$comp/"
      n=$((n + 1))
    done
  done
done

# exoplanet_mod.F90 is generated, not copied: the template plus the config's
# user_nl_exo plus the case's own user_nl_exo ($CASEROOT/user_nl_exo, next to
# DEST), when that exists. Rerun gen_exoplanet_mod.py after editing it.
if [ -f "$cfgdir/user_nl_exo" ]; then
  gen=$root/exoplanet/gen_exoplanet_mod.py
  mkdir -p "$dest/src.share"
  casenl=$(cd "$dest/.." && pwd)/user_nl_exo
  if [ -f "$casenl" ]; then
    python3 "$gen" "$config" --user-nl "$casenl" --out "$dest/src.share/exoplanet_mod.F90"
  else
    python3 "$gen" "$config" --out "$dest/src.share/exoplanet_mod.F90"
  fi
  n=$((n + 1))
fi
echo "install_sourcemods: $config -> $dest ($n files; components: $comps)"
