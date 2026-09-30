#!/bin/bash
# Builds both recipes and runs their tests on linux-64 in Docker, as bioconda CI
# does. Run from anywhere inside the clone.
#
#   scripts/linux_test.sh          # panteraga source = git archive of HEAD (before tagging)
#   scripts/linux_test.sh --tag    # panteraga source = the recipe as is (after tagging)
#
# The HEAD archive has exactly what the tag tarball will have, so uncommitted
# changes are NOT tested: commit first.
# Built packages go to $PANTERA_BLD (default ~/conda-bld-linux) and are reused.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BLD="${PANTERA_BLD:-$HOME/conda-bld-linux}"
USE_TAG=0
[ "${1:-}" = "--tag" ] && USE_TAG=1
mkdir -p "$BLD"

if [ "$USE_TAG" = 0 ]; then
  git -C "$REPO" archive --format=tar.gz --prefix=panteraGA/ -o "$BLD/panteraGA-HEAD.tar.gz" HEAD
fi

docker run --rm --platform linux/amd64 \
  -e USE_TAG="$USE_TAG" \
  -v "$REPO/recipes":/recipes:ro \
  -v "$BLD":/bld \
  condaforge/miniforge3 bash -euo pipefail -c '
    conda install -y -q -n base conda-build
    CH="-c file:///bld -c conda-forge -c bioconda"

    conda build /recipes/alntools $CH --output-folder /bld

    cp -r /recipes/panteraga /tmp/panteraga
    if [ "$USE_TAG" = 0 ]; then
      # First source: the HEAD archive instead of the tag tarball (and no sha256)
      sed -i -e "s#- url: https://github.com/piosierra/panteraGA/.*#- url: file:///bld/panteraGA-HEAD.tar.gz#" \
             -e "0,/sha256:/{/sha256:/d}" /tmp/panteraga/meta.yaml
    fi
    conda build /tmp/panteraga $CH --output-folder /bld

    conda create -y -q -n pantera $CH panteraga
    conda run -n pantera panteraGA -h
    echo "Linux build and tests OK"
  '
