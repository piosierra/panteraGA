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

# HMM libraries (model/hmm): from the clone when unpacked there (pantera_hmm_v1.tar.gz -> model/hmm), otherwise
# from the installed package, whatever its version. The classifier and model/label_map.tsv always come from the
# clone (next to the script).
shopt -s nullglob
if [ -d "$REPO/model/hmm" ]; then
  homes=("$REPO")
else
  homes=("$ENV_PREFIX"/share/panteraga-*)
fi
if [ ${#homes[@]} -ne 1 ]; then
  echo "Expected exactly one panteraga-* folder in $ENV_PREFIX/share, found ${#homes[@]}" >&2
  exit 1
fi

exec env -i HOME="$HOME" LANG=en_US.UTF-8 \
  PATH="$ENV_PREFIX/bin:/usr/bin:/bin:/usr/sbin:/sbin" \
  R_LIBS_USER=/nonexistent \
  PANTERA_HOME="${homes[0]}" \
  "$ENV_PREFIX/bin/Rscript" "$REPO/panteraGA" "$@"
