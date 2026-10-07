#!/bin/bash
# Compares the libraries of two panteraGA runs on the same data.
#
#   scripts/compare_runs.sh <run A folder> <run B folder> <library name (-b)>
#
# Reports processing time (after the FastGA alignments), loop 2 time when the
# timing files exist, element and pass counts, classes, which consensus sequences
# are shared, elements whose class changed, and, for each element found in only
# one run, its closest element in the other run (blastn). blastn must be on the
# PATH (e.g. PATH=~/miniconda/envs/pantera/bin:$PATH).
set -euo pipefail

if [ $# -ne 3 ]; then
  echo "Usage: $0 <run A folder> <run B folder> <library name>" >&2
  exit 1
fi
A=$1; B=$2; LIB=$3
FA="$A/$LIB-pantera-final.fa"; FB="$B/$LIB-pantera-final.fa"
for f in "$FA" "$FB"; do [ -s "$f" ] || { echo "Missing $f" >&2; exit 1; }; done
command -v blastn >/dev/null || { echo "blastn not found on the PATH" >&2; exit 1; }
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT

secs() {  # yy-mm-dd:HH:MM:SS -> seconds since midnight (runs within one day)
  echo "$1" | awk -F: '{print $2*3600 + $3*60 + $4}'
}
echo "== Processing time (from reading the alignments to the end)"
for d in "$A" "$B"; do
  s=$(grep -m1 "Procesing" "$d/pantera.log" | cut -c1-17)
  e=$(tail -1 "$d/pantera.log" | cut -c1-17)
  echo "  $d: $(( $(secs "$e") - $(secs "$s") )) s"
done

if [ -s "$A/loop2_mafft_timing.tsv" ] && [ -s "$B/loop2_mafft_timing.tsv" ]; then
  echo "== Loop 2 mafft time (sum over alignments)"
  for d in "$A" "$B"; do
    awk -F'\t' -v d="$d" 'NR==1{for(i=1;i<=NF;i++) if($i=="seconds") c=i; next} {s+=$c; n++} END{printf "  %s: %.0f s in %d alignments\n", d, s, n}' "$d/loop2_mafft_timing.tsv"
  done
fi

echo "== Elements"
for d in "$A" "$B"; do
  echo "  $d: $(grep -c '>' "$d/$LIB-pantera-final.fa") elements, $(grep -c '>' "$d/$LIB-pantera-final-pass.fa") pass"
done

echo "== Classes ($A | $B)"
grep '>' "$FA" | sed 's/.*#//' | sort | uniq -c | awk '{print $2"\t"$1}' | sort > "$T/ca"
grep '>' "$FB" | sed 's/.*#//' | sort | uniq -c | awk '{print $2"\t"$1}' | sort > "$T/cb"
join -t $'\t' -a1 -a2 -e 0 -o 0,1.2,2.2 "$T/ca" "$T/cb" |
  awk -F'\t' '{printf "  %-24s %4d %4d%s\n", $1, $2, $3, ($2 != $3 ? "   <-" : "")}'

# sequence <TAB> name <TAB> length, sorted by sequence
tab() { awk '/^>/{n=substr($0,2); next}{print $0"\t"n"\t"length($0)}' "$1" | sort -t $'\t' -k1,1; }
tab "$FA" > "$T/a"; tab "$FB" > "$T/b"
shared=$(join -t $'\t' "$T/a" "$T/b" | wc -l | tr -d ' ')
echo "== Consensus sequences: $shared identical in both"

echo "== Identical sequence, different class"
join -t $'\t' "$T/a" "$T/b" |
  awk -F'\t' '{split($2,x,"#"); split($4,y,"#"); if (x[2] != y[2]) print "  "$2" -> "$4}' | grep . || echo "  none"

only() {  # elements of run 1 not in run 2, with their closest element in run 2
  local name=$1 a=$2 b=$3 fb=$4
  join -t $'\t' -v1 "$a" "$b" | cut -f2,3,1 > "$T/u"
  local n; n=$(wc -l < "$T/u" | tr -d ' ')
  echo "== Only in $name: $n"
  [ "$n" -gt 0 ] || return 0
  awk -F'\t' '{print ">"$2"\n"$1}' "$T/u" > "$T/u.fa"
  # Best match per element by bit score; coverage from the query positions aligned
  blastn -query "$T/u.fa" -subject "$fb" -outfmt "6 qseqid sseqid pident qstart qend bitscore" -max_hsps 1 2>/dev/null |
    sort -t $'\t' -k1,1 -k6,6gr | awk -F'\t' '!seen[$1]++ {print $1"\t"$2"\t"$3"\t"($5 - $4 + 1)}' > "$T/hits"
  awk -F'\t' '{print $2"\t"$3}' "$T/u" | sort -t $'\t' -k1,1 > "$T/ul"
  sort -t $'\t' -k1,1 "$T/hits" > "$T/hs"
  join -t $'\t' -a1 "$T/ul" "$T/hs" |
    awk -F'\t' '{ if (NF >= 5) printf "  %-28s %6d bp  closest: %-28s %5.1f%% identity over %5.0f%% of its length\n", $1, $2, $3, $4, 100*$5/$2;
                  else printf "  %-28s %6d bp  closest: none\n", $1, $2 }'
}
only "$A" "$T/a" "$T/b" "$FB"
only "$B" "$T/b" "$T/a" "$FA"
