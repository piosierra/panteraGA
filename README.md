[![pantera](images/panteraGA.png?raw=true "pantera")](images/panteraGA.png?raw=true)

**Identification of transposable element families from whole genome alignments using FastGA**

panteraGA is the new version of [pantera](https://github.com/piosierra/pantera). It aligns two or more genomes with [FastGA](https://github.com/thegenemyers/FASTGA) and builds a library of transposable elements (TEs) from the polymorphic segments between them. It generates TE libraries in minutes for most species and can handle genomes over 10 Gb in a few hours.

Rather than an all-vs-all alignment, genomes in the list are aligned pairwise in a ring (1 vs 2, 2 vs 3, ..., last vs 1), and the polymorphic segments from all the alignments are processed together.

## Installation

### Conda (recommended)

```bash
conda create -n pantera -c conda-forge -c bioconda panteraga
conda activate pantera
```

This installs panteraGA with all its dependencies (FastGA, FasTAN, alntools, BLAST, MAFFT, CD-HIT, HMMER, miniprot, Infernal, TRF and the R packages) and the HMM libraries of the classifier.

### Containers

Every bioconda package is also available as a container from [BioContainers](https://biocontainers.pro). Find the current tag at [quay.io/repository/biocontainers/panteraga](https://quay.io/repository/biocontainers/panteraga?tab=tags), then:

```bash
singularity pull panteraga.sif docker://quay.io/biocontainers/panteraga:1.4.0--hdfd78af_0
singularity exec panteraga.sif panteraGA -V
singularity exec panteraGA.sif panteraGA -g genomes.txt -b mylib -o out
```

### From source

Cloning the repository is only needed for development. The HMM libraries of the classifier are not stored in git: download `pantera_hmm_v1.tar.gz` (47 MB) <!-- TODO before 2.0.0: Zenodo record of pantera_hmm_v1 --> and unpack it as `model/hmm/` (`tar -xzf pantera_hmm_v1.tar.gz && mv pantera_hmm_v1 model/hmm`). The optional `pantera_hmm_v1_host.tar.gz` (60 MB, Swiss-Prot host proteins, needs DIAMOND) adds the `host_protein` flag; unpack it next to the first one, before moving the folder. The dependencies listed in [`recipes/panteraga/meta.yaml`](recipes/panteraga/meta.yaml) must be on the `PATH`.

panteraGA reads the classifier (`scripts/pantera_classify.R`) and its label map (`model/label_map.tsv`) from its own folder, and the HMM libraries from `model/hmm/` next to the script, or in the folder given by the `PANTERA_HOME` environment variable. `scripts/test_classifier.sh <library folder>` checks the classifier against an exact expected output.

## Usage

```bash
panteraGA -g genomes.txt -b mylib -o out -T 16
```

- `genomes.txt`: a plain text file with the paths to the genomes (FASTA), one per line, at least two.
- `mylib`: an identifier appended to the names of the elements in the library.

### Options

| Option | Default | Description |
|---|---|---|
| `-g, --genomes` | required | File with the paths to the genomes, one per line (min 2) |
| `-b, --lib_name` | required | Identifier appended to element names |
| `-o, --output_folder` | `pantera_output` | Output folder |
| `-T, --threads` | 8 (or the number of cores, if fewer) | Threads |
| `-s, --min_size` | 50 | Min length of a polymorphic insertion, without flanks |
| `-l, --max_size` | 30000 | Max length of a polymorphic insertion, without flanks |
| `-i, --identity1` | 0.90 | Identity for clustering, first round (0.8 to 1) |
| `-y, --identity2` | 0.85 | Identity for clustering, second round (0.8 to 1) |
| `-m, --min_cluster` | 3 | Min sequences to build a consensus |
| `-u, --max_cluster` | 300 | Max sequences per cluster, first round (the second round uses 1000) |
| `-e, --min_copies` | 2 | Min copies in at least one genome for an element to be accepted |
| `-n, --max_ns` | 0.001 | Max fraction of Ns in a segment |
| `-c, --max_cons_ns` | 0.02 | Max fraction of Ns in a final consensus |
| `-p, --min_polya` | 10 | Min length of a 3' tail (polyA or short tandem repeat, unit 1-6 bp) |
| `-f, --flanking` | automatic | Min length of the flanks of a segment. If not given, the fraction `-q` of segments with the shortest flanks is discarded instead |
| `-q, --flank_quantile` | 0.05 | Without `-f`: fraction of segments with the shortest flanks to discard |
| `-a, --anno_coverage` | 0.80 | Min fraction of an element covered by a hit, for genome annotation |
| `-z, --anno_identity` | 0.80 | Min identity of a hit, and min length ratio hit/element, for genome annotation |
| `-I, --align_iterations` | 1000 | Max mafft refinement iterations (L-INS-i) for clusters of fewer than 100 sequences shorter than 3 kb; 1000 is what mafft `--auto` uses. `-I 2` is somewhat faster (about 10% of the processing time in a test with ~600 elements) but changes some consensi: in that test about 20% differed slightly and a few short elements were gained or lost |
| `-r, --seed` | 1 | Random seed for the TIR test. The same data and seed give identical results, whatever `-T` |
| `-N, --no_tandem_mask` | | Do not mask tandem arrays before aligning. By default the tandem arrays of each genome, found with [FasTAN](https://github.com/thegenemyers/FASTAN), are soft masked before aligning (`FastGA -M`): much faster with genomes rich in satellites, with almost the same segments (see [Tandem masking](#tandem-masking)). The genome databases and masks are written to `tandem/` in the output folder, and kept with `-k` or `-d`. `-t, --tandem_mask` is still accepted and does nothing |
| `-w, --tandem_min_len` | 1000 | Min length of the tandem arrays masked, in bp. Shorter arrays add little alignment time, but are found inside some TEs (e.g. in LTRs), whose copies masking can lose. `0` masks them all |
| `-k, --keep_alignments` | | Keep the FastGA alignments, named after the two genomes (`<genome1>_<genome2>.1aln`). A run reuses the alignments it finds in its output folder, including ones named `pantera<n>.1aln` by earlier versions |
| `-d, --debug` | | Keep intermediate files (all segments in `pantera_lib_0.fa`, the loop 1 output in `segments_candidates.fa`, `loop_1/`, `loop_2/`, the consensi before classification, the classifier's per-chunk tool outputs in `hmm_classifier/`), and write the loop 2 timing files (`loop2_window_timing.tsv`, `loop2_mafft_timing.tsv`; summarise them with `scripts/loop2_timing_summary.R`) |
| `-v, --verbose` | | Show log messages on screen |
| `-V, --version` | | Show the version and exit |
| `-h, --help` | | Show the options and exit |

### Output

In the output folder:

| File | Content |
|---|---|
| `mylib-pantera-final.fa` | The library: the TE consensus sequences accepted by the structural checks |
| `mylib-pantera-rejected.fa` | The consensus sequences rejected by the checks |
| `mylib-pantera-final.stats.tsv` | Per element, accepted or rejected: length, class, `pass` and the `reason` (the rule that accepted or rejected it), cluster size (`cluster_size`, copies) and `loci` (copies of the same insertion in several genomes count once), TSD, terminal repeats, tails, copies per genome... |
| `alignments/` | The alignment (`.maf`) behind each consensus, named as the element. Each copy shows its flanks in lowercase and the insertion in uppercase |
| `annotations/` | A BED file per genome with the copies of each accepted element |
| `pantera.log` | Log of the run |

Unless `-N` is given, `tandem/` holds the genome databases (`.1gdb`), the tandem arrays found by FasTAN (`<genome>_tan.1ano`) and the ones masked, at least `-w` bp long (`<genome>_tan<w>.1ano`), when kept with `-k` or `-d`. The alignments refer to the databases, so keep them together.

Elements are named `<class>_<n>-<lib_name>#<class>`, numbered by length within each class (accepted and rejected elements share the numbering).

## Classification

Each consensus is classified by homology with the HMM classifier (`scripts/pantera_classify.R`, libraries `pantera_hmm_v1`), then checked and, for some classes, reclassified by its structure (terminal repeats, TSDs, tails, ORFs; the `reason` column of the stats says which rule decided). The classifier combines, in this order of precedence:

1. protein domains: HMMs of TE protein domains (`hmmsearch --cut_ga`) on miniprot translations of curated TE proteins and on six-frame ORFs;
2. TIR end models (DNA transposons) and end partners (ends shared with an autonomous DNA transposon of the same library or of a curated set);
3. Helitron end models;
4. SINE head and body models.

The classifier's label (pantera classification) becomes the class written in the library through `model/label_map.tsv`, which gives RepeatMasker names (e.g. `LINE/Jockey` -> `LINE/I-Jockey`, `DIRS/DIRS` -> `LTR/DIRS`). When a sub label is reported and RepeatMasker has a name for it, the name is refined (e.g. `DNA/Tc1-Mariner` with sub label Pogo -> `DNA/TcMar-Pogo`). Without a call the class is `Unknown`. The classifier can also run alone on any library: `Rscript scripts/pantera_classify.R -i library.fa -o out -L model/hmm -T 8` (output `out/classification.tsv`).

Homology columns of the stats file:

| Column | Content |
|---|---|
| `homology_label` | Label of the classifier (pantera classification, `TYPE/subtype`), or `none` |
| `homology_confidence` | `high`, `low`, `weak` or `none` (below) |
| `homology_evidence` | What made the call: `protein`, `TIR` (TIR end models), `partner` (shared ends), `TIR+partner`, `Helitron` (end models), `SINE` |
| `homology_sub` | Sub label (lineage below the subtype, e.g. Tekay, Pogo); not reported for protein sub labels whose label has a single sub-label model (`class_cfg$min_sub_models`) |
| `homology_domains` | Domain classes of the protein call (e.g. `GAG,INT,RH,RT`; abbreviations below) |
| `homology_alternatives` | For `weak` calls: other labels that share the only domain found |
| `homology_flags` | Warnings: `tandem` (TRF arrays cover >= 50%), `low_complexity` (dustmasker >= 30%), `ncRNA:<families>` (Rfam hits >= 30%, or any rRNA), `host_protein:<Swiss-Prot entry>` (only with the optional host part), `chimera:<what>` (strong evidence for two TYPEs) |
| `homology_notes` | `single_domain:<class>` for weak calls; `close_runner_up:<label>` when the call leads the next label by < 50 bits; `partner_conflict:<label>` when TIR ends and partner disagree |
| `homology_ends` | DNA transposon class of the element's ends (TIR end models, or end partner) when they match both ends, whatever the final call; empty otherwise |

Confidence: `high` and `low` are the protein rule's confidence (`low`: one domain class below 100 bits, or RH, YR, CHD or HEL alone), TIR + partner of the same superfamily is `high`, TIR alone keeps its own, partner alone is `low`. `weak`: protein calls whose only domain class is shared with other TYPEs (LTR integrase, RNase H or low reverse transcriptase alone, DIRS tyrosine recombinase alone, PLE GIY-YIG alone), and Helitron-end calls; structure must confirm them. In cross-validation on curated elements, protein calls were right in 99.98% (high), 99.7% (low) and 72% (weak) of cases.

Domain classes: AP aspartic protease, ATPase, CHD chromodomain, ENDO endonuclease, ENV envelope, GAG capsid, GIY GIY-YIG nuclease, HEL helicase, INT integrase, POLB DNA polymerase B, REP Rep nuclease, RH RNase H, RT reverse transcriptase, TPase DDE transposase, YR tyrosine recombinase.

## Element edges, TSDs and partners

Copies from different genomes whose flanks next to the insertion (60 bp each side) share 16-mers on both sides are the same insertion found in several genomes (one locus, an orthologous copy). When a cluster has at least `-m` loci, one copy per locus is aligned: orthologous copies have conserved flanks, which would look like a segmental duplication and bias the TSD. Unknown elements found at more than 5 loci are accepted ("Unknown: many loci"; `rules_cfg$unknown_min_loci`).

The edges of each element are where most of its copies' svfind segments start and end: svfind finds each insertion against the other genome, without its target site duplication (TSD). They are moved by at most one base to where the copies agree best on a TSD. When no copy has svfind boundaries, the edges come from the conservation of the alignment and a search for the TSD.

The TSD is called from both sources. svfind sees the target of each insertion in the empty allele, so the bases of the empty allele found at both ends of the insertion are its TSD (confirmed by the copy's own flanks). The flank caller compares the bases next to the edges in every copy. When at least two copies have an svfind TSD, the element has one, of the flank caller's length if it is within the copies' range, otherwise of their median. Without svfind support, the flank caller's TSD is reported only if the copies' own flanks match better than flanks of different copies do (cross-copy null). TSDs up to 20 bp are tested.

An Unknown element whose two termini are those of an element classified by homology takes that element's class: as a non-autonomous derivative (`-NA` added to the class, e.g. `DNA/TcMar-NA`) if it has the partner's edges (allowing up to 3 A, T or N bases of microsatellite target at either end) and terminal repeat at high identity, as a low-confidence call (`NA_low`) otherwise. Solo LTRs and internal deletions of LINEs take the partner's class and are rejected.

Elements are oriented by the strand of their protein call; without one, as the longer element (with a protein call) whose 3' end they match, if any (solo LTRs, 5'-truncated LINE copies); otherwise by their longest ORF. The stats column `orientation` says which: `protein`, `container` (turned to the longer element's strand), `tail` (a SINE turned so that its tail is 3') or `ORF` (including fragments already on their container's strand).

Tails (stats `tail3_motif`, `tail3_len`, `tail5_motif`, `tail5_len`; `polyA` is TRUE when either is an A tail) are short tandem repeats (unit 1-6 bp, polyA = unit 1) of at least `-p` bp ending within 10 bp of an end. A LINE needs a long ORF and a tail at its 3' end (or no terminal repeat). SINEs come only from the classifier's SINE models; a SINE with its tail only at the 5' end is turned (`orientation` `tail`). Elements contained in a longer accepted LINE, LTR, PLE or RC element (one hit over 98% of the element at 98% identity, tails excluded), and Unknowns that are 3' fragments of an accepted LINE or PLE (90% at 90%), take the container's class and are rejected ("contained in"). Elements with no homology call that the classifier flags both `tandem` and `low_complexity` are microsatellites: class `Simple_repeat`, rejected ("microsatellite"). DNA elements with a superfamily call (protein, TIR end models or end partner, high or low confidence) are accepted without the structural TIR test. An element whose ends are those of a known DNA transposon (`homology_ends`) but whose call is another class (e.g. a reverse transcriptase inside) is a non-autonomous DNA transposon of that superfamily: class `DNA/<superfamily>-NA`, with the inner class named in `reason`. Non-autonomous elements found by shared termini (`partner_type` `NA_low`) are accepted too. The `-NA` LTR elements with terminal LTRs of their own are then partners for solo-LTR calls, so that solo LTRs of non-autonomous LTR elements are rejected as such.

TSD and partner columns of the stats file:

| Column | Content |
|---|---|
| `TSD_length`, `TSD_confidence`, `TSD_motif`, `TSD_pvalue` | TSD length (bp, 0 = none), fraction of copies supporting it, IUPAC motif, p-value (flank caller and cross-copy null; empty for svfind lengths) |
| `TSD_svfind` | Most frequent svfind empty-site overlap of the copies (bp) |
| `TSD_svfind_support` | Fraction of copies with that overlap |
| `TSD_source` | What gave the TSD: `flanks` (flank caller, within the svfind range), `svfind` (median svfind TSD), `null` (flank caller against the cross-copy null), `none` |
| `TSD_target_tandem` | Fraction of copies inserted into a short tandem repeat (both flanks repeat with a period up to 6 bp); from 0.5 the TSD is not used to reclassify |
| `TSD_superfamily` | DNA elements only: superfamily suggested by the TSD alone (low confidence; the class is not changed): `DNA/TcMar` (TA; any TcMar-*), `DNA/PIF-Harbinger|DNA/PIF-ISL2EU` (AT), `DNA/hAT` (8 bp; any hAT-*), `DNA/PiggyBac|DNA/Kolobok` (TTAA); empty for other TSDs and for targets in short tandem repeats |
| `partner` | Element of the run whose termini this one shares |
| `partner_type` | `NA_high` (non-autonomous, `-NA` class), `NA_low` (same, low termini identity), `solo_LTR`, `internal_deletion` |
| `partner_ident` | Identity over the partner's terminal repeat (TIR or LTR; 60 bp without one), the lower of the two termini |

## Tandem masking

Long tandem arrays (satellites, centromeres) make FastGA extend huge alignments through them, and can dominate the alignment time. By default (`-N` turns it off), FasTAN finds the tandem arrays of each genome and those at least `-w` bp long (1 kb by default) are soft masked in the genome index, so FastGA starts no alignment inside them (alignments can still extend through them). Coordinates and sequences are unchanged: segments, alignments and annotations refer to the original genomes.

Short arrays are left unmasked because they add little alignment time, but are found inside some TEs: masking them all lost, in a test, a young Gypsy element with 51 bp tandem repeats in its LTRs.

In a test with the two haplotypes of gorilla chromosome 22 (T2T, 59 Mb, 36% in tandem arrays, 2 threads), FastGA took 632 s without masking, 20 s masking all the arrays (8,952 in one haplotype) and 21 s masking those of 1 kb or more (363, holding 95% of the arrays' length). svfind found 201 segments without masking; masking all arrays kept 176 of them, and masking arrays of 2 kb or more 186. The segments lost were at tandem arrays (likely array length differences). In a 25 Mb genome with 6% in tandem arrays, masking made no difference in time.

## Citation

SOON?

The HMM libraries of the classifier: <!-- TODO before 2.0.0: Zenodo DOI of pantera_hmm_v1 --> (up to 1.4.1, the XGBoost classification model: doi:[10.5281/zenodo.22990589](https://doi.org/10.5281/zenodo.22990589))

## Acknowledgements

Thanks to Arian Smit and Robert Hubley (Dfam, RepeatMasker, RepeatModeler) for allowing us to use their curated peptide library in releases up to 1.4.1.
Thanks to Claude for helping with the development.

## Licence

MIT, see [LICENSE](LICENSE).













