#!/usr/bin/env bash
# Test of the HMM classifier (scripts/pantera_classify.R) with exact expected output.
# Input test/hmm/input.fa.gz: 48 sequences (36 consensi of the VGP library aPleWal, up to 6 for each evidence x
# confidence class; 12 decoys: mRNAs, Rfam genes, shuffled windows). Expected output: test/hmm/expected_classification.tsv
# (made by the Python caller of the classifier project, step 12). Run twice: default chunks, and chunks of 7
# sequences (the result must not depend on the chunk size).
# usage: scripts/test_classifier.sh <HMM library folder (pantera_hmm_v2)> [threads (default 4)]
# Needs on the PATH: Rscript (data.table, getopt), hmmsearch, nhmmer, miniprot, cmscan, cmfetch, cmalign, trf,
# dustmasker, makeblastdb, blastn (and diamond if the library has the host-protein part).
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
LIB=$1; T=${2:-4}
[ -s "$LIB/prot/pantera_v5.cal.hmm.gz" ] || { echo "not an HMM library folder: $LIB"; exit 1; }
missing=""
for t in Rscript hmmsearch nhmmer miniprot cmscan cmfetch cmalign trf dustmasker makeblastdb blastn; do
  command -v "$t" > /dev/null || missing="$missing $t"
done
[ -z "$missing" ] || { echo "test_classifier: not on the PATH:$missing"; exit 1; }
echo "test_classifier: Rscript $(command -v Rscript)"
W=$(mktemp -d)
gzip -dc "$REPO/test/hmm/input.fa.gz" > "$W/input.fa"
fail=0
for C in "" 7; do
  O="$W/out${C:+_c$C}"
  cmd=(Rscript "$REPO/scripts/pantera_classify.R" -i "$W/input.fa" -o "$O" -L "$LIB" -T "$T" ${C:+-C $C})
  rc=0; "${cmd[@]}" > "$O.log" 2>&1 || rc=$?
  if [ "$rc" -ne 0 ]; then
    tail -20 "$O.log"
    echo "test_classifier: FAIL, exit status $rc of: ${cmd[*]}"
    echo "test_classifier: full log (kept): $O.log"
    exit 1
  fi
  if cmp -s "$O/classification.tsv" "$REPO/test/hmm/expected_classification.tsv"; then
    echo "test_classifier${C:+ (chunks of $C)}: PASS ($(($(wc -l < "$O/classification.tsv") - 1)) sequences)"
  else
    diff "$REPO/test/hmm/expected_classification.tsv" "$O/classification.tsv" | head -10
    echo "test_classifier${C:+ (chunks of $C)}: FAIL"; fail=1
  fi
done
rm -rf "$W"
exit $fail
