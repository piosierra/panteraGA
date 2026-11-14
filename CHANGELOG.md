# Changelog

## [Unreleased]

### Changed (classification)
- The classification by homology uses the HMM classifier (`scripts/pantera_classify.R`, libraries `pantera_hmm_v1`) instead of blastp against RepeatMasker's peptide library and the XGBoost model. Each consensus gets one call from protein domain HMMs, TIR end models and end partners, Helitron end models or SINE models, with a confidence (`high`, `low`, `weak`, `none`). The call's label becomes the class through `model/label_map.tsv` (RepeatMasker names, refined by the sub label when RepeatMasker has a name for it, e.g. `DNA/TcMar-Pogo`, `DNA/hAT-hobo`); without a call the class is `Unknown`. Protein sub labels whose label has a single sub-label model are not reported (`class_cfg$min_sub_models`). The structure rules are unchanged.
- Stats file: the columns `homology_class` and `homology_prob` are replaced by `homology_label`, `homology_confidence`, `homology_evidence`, `homology_sub`, `homology_domains`, `homology_alternatives`, `homology_flags` and `homology_notes` (README, "Classification").
- Results, with the same alignments. Test data: 6 of 17 elements change class, none changes `pass`: 5 elements called `DNA` by their TIRs now have a superfamily (CMC-Transib, P, TcMar-Pogo, TcMar-Tc1, hAT-hobo), and one hAT becomes hAT-hobo. Five D. melanogaster genomes (168 elements, compared with RepBase 31.10): 31 change class and 3 change `pass`. Gains: 6 `DNA` elements get their superfamily (P, MULE-MuDR); 8 of the 10 elements called TcMar-Tc1 (the other 2 stay TcMar-Tc1) and 1 called `DNA` are TcMar-Pogo (RepBase POGO, POGON1), and 3 hAT are hAT-hobo; 7 Jockey fragments that were Unknown are LINE/I-Jockey (still rejected). Losses: 3 Jockey consensi with ORF1 only and 1 shorter one get no call, because no HMM covers LINE ORF1: 2 that were accepted as LINE are rejected, and one is accepted as "Unknown: large cluster". 2 BURDOCK consensi (GAG hit below threshold) stay accepted, as `LTR` instead of `LTR/Gypsy`. Accepted elements with the right class: 94 of 96 (was 96 of 97); right subclass: 92 of 96 (was 90 of 97). Classification took about as long as before (171 s against 152 s for 195 consensi, 10 threads).
- No longer used: getorf (EMBOSS), the R package xgboost, `libs/RepeatPeps.lib*` and `model/{typesnames,featurenames,xgbmodel.ubj}`. New tools needed: HMMER, miniprot, Infernal, TRF (and DIAMOND for the optional host-protein part). The recipe and the repository files change together with the Zenodo record of `pantera_hmm_v1` (next patch).
- `-d` also keeps the classifier's per-chunk files (`hmm_classifier/`).
- `scripts/run_dev.sh` uses the clone's `model/hmm/` when it is there.

### Changed
- Output: `<lib>-pantera-final.fa` is now the library of accepted elements only, and the rejected ones are written to `<lib>-pantera-rejected.fa`. `<lib>-pantera-final-pass.fa` is no longer written (its content is now `-final.fa`). The stats file still lists every element, accepted or rejected. The genome annotation (`annotations/`) uses only the accepted elements. Element names and numbering are unchanged.
- `scripts/compare_runs.sh` compares the accepted libraries, for runs of this and earlier versions.

### Added
- `scripts/pantera_classify.R`: the HMM classifier of TE sequences (protein domain HMMs, TIR and Helitron end models, SINE models, end partners, and flags for tandem, low-complexity, ncRNA and host-protein sequences), an R/data.table port of the unified caller of the panteraGA HMM classifier project (step 12). Used by panteraGA for the classification (see "Changed (classification)"); it also runs alone (`Rscript scripts/pantera_classify.R -i library.fa -o out -L pantera_hmm_v1 -T 8`) and writes `classification.tsv`, one row per sequence. Its output is identical, byte for byte, to that of the Python caller on the 48 test sequences, the VGP library aPleWal (999 consensi), 1,200 decoys and 12 more libraries. Results do not depend on the number of threads or the chunk size. The libraries (`pantera_hmm_v1.tar.gz`, 47 MB; optional host-protein part `pantera_hmm_v1_host.tar.gz`, 60 MB) are not in the repository. `scripts/test_classifier.sh <library folder>` checks the exact expected output on `test/hmm/` (input and expected classification). It needs HMMER, miniprot, Infernal, TRF and BLAST (and DIAMOND for the host-protein flag) on the `PATH`; they are not yet in the recipe.
- `-t/--tandem_mask`: the tandem arrays of each genome, found with FasTAN, are soft masked before aligning (FastGA `-M`), so FastGA starts no alignment in satellites and centromeres. Only arrays of at least `-w/--tandem_min_len` bp (default 1000; 0 masks all) are masked: short arrays add little alignment time, but are found inside some TEs, and masking them lost a young Gypsy with tandem repeats in its LTRs. With the two haplotypes of gorilla chromosome 22 (T2T), the alignment took 21 s instead of 632 s, masking 363 arrays of 1 kb or more (95% of the arrays' length); masking all 8,952 arrays took 20 s but kept fewer of the segments found without masking (176 of 201, against 186 masking arrays of 2 kb or more). Coordinates and sequences are those of the original genomes. The genome databases (`.1gdb`), FasTAN's arrays (`_tan.1ano`) and the masks (`_tan<w>.1ano`) are written to `tandem/` in the output folder, and kept with `-k` (the alignments refer to them) or `-d`. Needs FasTAN (new dependency, bioconda `fastan`); FAtoGDB and GIXmake come with FastGA.
- Stats file: column `reason`, after `pass`: the rule that accepted or rejected the element (the last rule that decided, since rules run in order and later ones overwrite earlier ones).
- Stats file: column `TSD_pvalue`, after `TSD_motif`: the p-value of the reported TSD (agreement between copies beyond chance, corrected for the TSD lengths tried). Empty when no TSD is reported. The consensus names in intermediate files carry it as `_tsdp`.

### Changed (rules)
- Contained elements are rejected as the last step, so no rule can accept them again: an element covered over at least 98% of its length, at 98% identity or more, by one hit of a longer accepted LINE, LTR, PLE or RC element (not DNA). The reason names the container. This replaces the previous step, which applied only to elements contained in LINEs (at more than 98%) and could be overridden by later rules.
- LINEs that are mostly a repeat of themselves (longest self-match, `maxTR`, over 40% of their length) are rejected, as LTR-structure elements already were. This replaces the rule that accepted every LINE with `maxTR` over 50 bp, after all other filters (including Ns and containment).
- Removed: the rule accepting any element with TIRs within 3 bp of both ends, whatever its class (DNA and Unknown elements with TIRs are accepted by their own rules), and the rule making Unknowns with TIRs and an 8 bp TSD DNA (the next rule made every Unknown with TIRs DNA anyway).

### Fixed
- FastGA, GIXmake and FasTAN are run up to 3 times when they fail. FastGA's index step (GIXmake) sometimes crashes with a segmentation fault (seen with bioconda's FastGA 1.5.20260729, about once in a few runs of the test genomes) and stopped the whole run; run again, it works. Before a new attempt, the files left by the failed one are removed: the partial alignment, and the GDB, mask and index files FastGA had made next to the genomes (files that were there before the run are kept). Each retry is logged.
- `-t` with no tandem array at least `-w` bp long in a genome: GIXmake could not read the empty mask and the run stopped. That genome is now indexed without a mask.
- `-d/--debug` kept almost no intermediate files: the segments (`pantera_lib_0.fa`), the loop 1 output (`segments_candidates.fa`), `loop_1/`, `loop_2/` and the consensi were deleted at the end of every run. They are now kept with `-d`.
- `-i/--identity1` and `-y/--identity2` must be between 0.8 and 1. cd-hit-est refuses lower identities with the word length panteraGA uses, and the run stopped later with a misleading error ("cannot open the connection"). Values outside the range now stop the run at the start, with a message naming the option.
- A TSD length found in a single copy no longer counts as a TSD. A long chance match in one copy (e.g. 11 bp in 1 of 5 copies) passed the significance test, and its length was reported with confidence 0 and no motif. The minimum number of copies sharing a TSD length is now 2 (`tsd_cfg$min_support`, also used for the motif). In the test data, two SINEs lose such 10-11 bp "TSDs"; nothing else changes.

## [1.4.1]

Libraries, stats, alignments and annotations are identical to 1.4.0, and runs are faster: in a two-genome test (Eledone cirrhosa, 8 threads), processing after the FastGA alignments took about 30% less time. Two rare cases that stopped the run with an error now complete (see Fixed). The code is reorganised internally.

### Fixed
- No more R warnings "'/etc/localtime' is not a symlink" / "is not identical to any known timezone file" when running in a container: if `TZ` is unset and `/etc/localtime` is a regular file, panteraGA uses it directly (`TZ=:/etc/localtime`), keeping the host's local time.
- Loop 2 no longer stops the run when a cluster's alignment has no conserved block, or when a cluster holds a single sequence (possible only with `-m 1`): the cluster gives no consensus. Both cases used an R `next` outside a loop and would have stopped the run with an error; neither happened in the test runs.
- A file named `all_segments.fa` in the directory the run was started from (not the output folder) made the run stop with an error.
- The log overcounted removed duplicates ("Removed dups"): an element duplicated in several pairs was counted once per pair. Only the logged number was wrong.
- Log typos: "Processing N alignments" (was "Procesing"), "successfully". `scripts/compare_runs.sh` accepts both spellings.

### Changed
- The FastGA alignments are named after the two genomes they align: `<genome1>_<genome2>.1aln` (file names without folders or extensions), instead of `pantera<n>.1aln`. Alignments named `pantera<n>.1aln` by earlier versions are still reused (renamed) when found in the output folder. Genome files whose names are the same without folders and extensions (e.g. `a.fa.gz` and `a.fna`) are now refused.
- Internal thresholds that are not command-line options (loop 1 and 2 windows, consensus, edge and TSD search, classification probability, structure rules) are grouped in named settings near the top of the script, each with a short description. Results are unchanged.
- Loop 2's per-cluster work (`process_cluster_2`, about 450 lines) is split into separate functions: alignment, column consensus, raw edges, the TSD caller, the edge search, the TSD offset refinement and writing the alignment. Results are unchanged.
- panteraGA no longer changes the working directory (`setwd`): every file is written with its full path in the output folder, and file paths passed to external tools are quoted. Functions no longer assign global variables (`<<-`): the genome names and the segment flank lengths are set once in the main part of the script. Results are unchanged.
- The terminal TIR test is about 8 times faster: the 100 shuffled copies it compares each element with are aligned all at once, and only their scores are computed. Results are unchanged (same shuffles, same scores, same p-values).
- BLAST (classification, structure checks, annotation) splits the work between threads by query (`-mt_mode 1`). By default BLAST splits by database, and with one small database most threads stayed idle (annotation used 1.35 of 8 threads). Same hits, in the same order. Requires BLAST+ 2.12 or later (recipe updated).

### Removed
- Unused function `cdhit1` and commented-out code.

## [1.4.0]

### Added
- Checks before aligning: every genome in the list must exist and appear only once. The run stops if FastGA produces no alignment for a pair.
- With `-d`, loop 2 writes timing files per window and per alignment (`loop2_window_timing.tsv`, `loop2_mafft_timing.tsv`); `scripts/loop2_timing_summary.R` summarises them.
- Option `-I/--align_iterations`: maximum mafft refinement iterations for small clusters (< 100 sequences, < 3 kb), which mafft `--auto` aligns with L-INS-i and up to 1000 iterations. The default (1000) keeps those alignments unchanged. `-I 2` is somewhat faster (about 10% of the processing time in a test with ~600 elements) but changes some consensi slightly, and a few short elements are gained or lost.
- `-r/--seed` (default 1): runs are reproducible. The TIR test compares each candidate with shuffled copies, and without a fixed seed borderline TIRs could be detected in one run and not in the next. Each element is tested with random draws seeded from its own sequence and `--seed`, so the same data and seed give identical results, whatever `-T`, and an element's result does not depend on the other elements.

### Changed
- The run stops with an error when any parallel task fails (an R error, or a worker that ends without a result, for example killed for using too much memory), naming the step and the first failing window or cluster. Failed tasks used to be skipped silently, losing their segments or clusters.
- Loop 2 clusters each length window, then aligns all clusters as separate parallel tasks, largest first. Previously each window aligned its own clusters one after another, which left a few windows running alone at the end. Same results; loop 2 took about 40% less time on 8 genomes with 8 threads, and now scales with `-T`.
- mafft is called directly instead of through the R package ips, which is no longer a dependency (it brought 38 other R packages). Alignments are identical.
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
