# Changelog

## [Unreleased]

### Changed
- In the alignment files (`alignments/*.maf`), each copy shows its flanks in lowercase and the insertion itself in uppercase. The library and the stats are not affected.

### Fixed
- Cluster members could be silently lost when reading cd-hit's cluster files, if the lines R sampled to guess the file's columns held only single-sequence clusters (shown as the R warning "Stopped early ... Expected 3 fields but found 4"). Clusters could then fall below `-m` and be discarded. The files are now read line by line. This affected all earlier versions.
- `-s/--min_size` and `-l/--max_size` now apply to the insertion itself. They were applied to the insertion plus its 200 bp of flanks, so `-s 100` let through insertions from 50 bp (the minimum passed to svfind). The default `-s` is now 50, which keeps the same segments as before; values below 50 are also passed to svfind. Consensus sequences must still be at least 100 bp long (or `-s`, if larger), as before.
- The flanks removed before clustering are the exact ones reported by svfind (its `X` line) instead of a fixed 100 bp per side, which also left 1 flank base on each side.
- `-n/--max_ns` had no effect: Ns were counted in uppercase on sequences that are still lowercase at that point. (With current FastGA versions segments contain no Ns, so results do not change.)
- The "Largest/Smallest insertion" log lines report the insertion length (previously "segment", with the round 2 values inconsistent).

## [1.3.1]

### Fixed
- When the genome list gave files with a folder (e.g. `genomes/sp1.fa` or an absolute path), clustering failed in most length windows and those segments were silently lost; typically all polymorphisms shorter than ~1.5 kb. Genome lists with plain file names were not affected.
### Changed
- The `name` column of `*-pantera-final.stats.tsv` no longer starts with `>`, so it matches the element names directly.

## [1.3.0]

First version distributed through bioconda (`conda install -c conda-forge -c bioconda panteraga`).
The XGBoost model is now downloaded from Zenodo (doi:10.5281/zenodo.22990589) and is no longer stored in git.

### Option names (breaking)
Long option names were made consistent; short flags are unchanged.

| Old | New |
|---|---|
| `--identity` | `--identity1` |
| `--min_cl` | `--min_cluster` |
| `--cl_size` | `--max_cluster` |
| `--mingen` | `--min_copies` |
| `--Ns` | `--max_ns` |
| `--cons_Ns` | `--max_cons_ns` |
| `--pAs` | `--min_polya` |
| `--anno_per` | `--anno_coverage` |
| `--anno_div` | `--anno_identity` |
| `--keep` | `--keep_alignments` |

- `-f/--flanking` is automatic when not given; `-f 100` now means a fixed 100 bp (it used to mean automatic).
- New `-V/--version`.
- `-h` lists every option with its default.

### Packaging
- R packages are no longer installed at runtime; all dependencies come from conda. `qualV` removed.
- The install location is resolved through symlinks; `PANTERA_HOME` overrides it.
- Startup checks for all required tools and data files, including detection of a Git LFS pointer in place of the model.
- `-h` exits with status 0; errors exit with status 1 and are reported on stderr.
- Portable shell redirection (`>/dev/null 2>&1`).

### Fixed
- `-y/--identity2` was ignored (round 2 clustered at `-i`).
- `-p` (polyA length) was ignored.
- `-n` is now a double (was declared integer).
- `-T` now also limits BLAST; mafft runs single-threaded inside the parallel workers.
- The flanking filter was computed but never applied; `-f` other than 100 crashed.
- LTR reclassification used the TIR table.
- Pass logic is now explicitly "fail unless rescued".
- Crashes on edge cases: no ORF hits, no self-BLAST hits, empty svfind output, `kneedle` failures, clusters without a conserved block.
- Per-consensus data (cluster size, TSD, genome counts) could be swapped between elements of equal length.
- Element names now match the final class (`<class>_<n>-<lib>#<class>`), numbered by length within class; alignment files renamed to match and reverse-complemented when the consensus is flipped (`_R_` prefix toggled).
- The null base composition for TSD statistics was always 0.25; it is now computed from the flanks.

### Changed
- Flanking filter: discards the fraction `-q` of segments with the shortest flank (the previous rule removed ~95% of segments).
- Consensus ties: transition ties (A/G, C/T) resolve to A/C; other ties become N. Consensus threshold 0.3.
- Edge search: inward shifts cannot exceed the TSD length.
- TSD support is tested for significance (binomial, corrected for the lengths and edge positions tried).
- Edge/TSD search rewritten (13–20x faster, identical results).
- Terminal TIR finder (anchored local alignment with a shuffled null, p <= 0.01) replaces the `short_tir` heuristic.
- Unknown elements with TG...CA termini and a 4–6 bp TSD (confidence > 0.8) are reclassified as LTR.
- DNA elements with a TIR within 8 bp of both ends pass.
- LINEs pass with a 3' tail (polyA or short tandem repeat) or 3'-anchored relatives.

### Added
- Options `-q/--flank_quantile` (default 0.05) and `-c/--max_cons_ns` (default 0.02).
- Stats columns `homology_class`, `homology_prob`, `tail3_motif`, `tail3_len`, `tail5_motif`, `tail5_len`, `tr3_relatives`, `fragment_of`.

## [1.1.0]
- Previous release (conda-free installation, model in Git LFS).
