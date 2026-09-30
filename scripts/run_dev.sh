#!/bin/bash
# Runs the panteraGA script from this clone using only the tools and data of an
# installed panteraga conda environment (clean PATH, no user R library).
#
#   scripts/run_dev.sh -g genomes.txt -b mylib -o out ...
#
# PANTERA_ENV selects the environment (default: ~/miniconda/envs/pantera).
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_PREFIX="${PANTERA_ENV:-$HOME/miniconda/envs/pantera}"

if [ ! -x "$ENV_PREFIX/bin/Rscript" ]; then
  echo "Rscript not found in $ENV_PREFIX/bin (set PANTERA_ENV)" >&2
  exit 1
fi

# Data (libs/, model/) from the installed package, whatever its version
shopt -s nullglob
homes=("$ENV_PREFIX"/share/panteraga-*)
if [ ${#homes[@]} -ne 1 ]; then
  echo "Expected exactly one panteraga-* folder in $ENV_PREFIX/share, found ${#homes[@]}" >&2
  exit 1
fi

exec env -i HOME="$HOME" LANG=en_US.UTF-8 \
  PATH="$ENV_PREFIX/bin:/usr/bin:/bin:/usr/sbin:/sbin" \
  R_LIBS_USER=/nonexistent \
  PANTERA_HOME="${homes[0]}" \
  "$ENV_PREFIX/bin/Rscript" "$REPO/panteraGA" "$@"
