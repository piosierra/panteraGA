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

This installs panteraGA with all its dependencies (FastGA, FasTAN, alntools, BLAST, MAFFT, CD-HIT, EMBOSS and the R packages) and the classification model.

### Containers

Every bioconda package is also available as a container from [BioContainers](https://biocontainers.pro). Find the current tag at [quay.io/repository/biocontainers/panteraga](https://quay.io/repository/biocontainers/panteraga?tab=tags), then:

```bash
singularity pull panteraga.sif docker://quay.io/biocontainers/panteraga:1.4.0--hdfd78af_0
singularity exec panteraga.sif panteraGA -V
singularity exec panteraGA.sif panteraGA -g genomes.txt -b mylib -o out
```

### From source

Cloning the repository is only needed for development. The classification model is not stored in git: download `model.tar.gz` from [Zenodo](https://zenodo.org/records/22990589) and unpack `xgbmodel.ubj` into `model/`. The dependencies listed in [`recipes/panteraga/meta.yaml`](recipes/panteraga/meta.yaml) must be on the `PATH`.

panteraGA looks for `libs/` and `model/` next to the script, or in the folder given by the `PANTERA_HOME` environment variable.

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
| `-p, --min_polya` | 10 | Min length of a polyA tail |
| `-f, --flanking` | automatic | Min length of the flanks of a segment. If not given, the fraction `-q` of segments with the shortest flanks is discarded instead |
| `-q, --flank_quantile` | 0.05 | Without `-f`: fraction of segments with the shortest flanks to discard |
| `-a, --anno_coverage` | 0.80 | Min fraction of an element covered by a hit, for genome annotation |
| `-z, --anno_identity` | 0.80 | Min identity of a hit, and min length ratio hit/element, for genome annotation |
| `-I, --align_iterations` | 1000 | Max mafft refinement iterations (L-INS-i) for clusters of fewer than 100 sequences shorter than 3 kb; 1000 is what mafft `--auto` uses. `-I 2` is somewhat faster (about 10% of the processing time in a test with ~600 elements) but changes some consensi: in that test about 20% differed slightly and a few short elements were gained or lost |
| `-r, --seed` | 1 | Random seed for the TIR test. The same data and seed give identical results, whatever `-T` |
| `-t, --tandem_mask` | | Soft mask the tandem arrays of each genome, found with [FasTAN](https://github.com/thegenemyers/FASTAN), before aligning (`FastGA -M`). Much faster with genomes rich in satellites, with almost the same segments (see [Tandem masking](#tandem-masking)). The genome databases and masks are written to `tandem/` in the output folder, and kept with `-k` or `-d` |
| `-w, --tandem_min_len` | 1000 | With `-t`, min length of the tandem arrays masked, in bp. Shorter arrays add little alignment time, but are found inside some TEs (e.g. in LTRs), whose copies masking can lose. `0` masks them all |
| `-k, --keep_alignments` | | Keep the FastGA alignments, named after the two genomes (`<genome1>_<genome2>.1aln`). A run reuses the alignments it finds in its output folder, including ones named `pantera<n>.1aln` by earlier versions |
| `-d, --debug` | | Keep intermediate files, and write the loop 2 timing files (`loop2_window_timing.tsv`, `loop2_mafft_timing.tsv`; summarise them with `scripts/loop2_timing_summary.R`) |
| `-v, --verbose` | | Show log messages on screen |
| `-V, --version` | | Show the version and exit |
| `-h, --help` | | Show the options and exit |

### Output

In the output folder:

| File | Content |
|---|---|
| `mylib-pantera-final.fa` | The library: the TE consensus sequences accepted by the structural checks |
| `mylib-pantera-rejected.fa` | The consensus sequences rejected by the checks |
| `mylib-pantera-final.stats.tsv` | Per element, accepted or rejected: length, class, `pass` and the `reason` (the rule that accepted or rejected it), cluster size, TSD, terminal repeats, tails, copies per genome... |
| `alignments/` | The alignment (`.maf`) behind each consensus, named as the element. Each copy shows its flanks in lowercase and the insertion in uppercase |
| `annotations/` | A BED file per genome with the copies of each accepted element |
| `pantera.log` | Log of the run |

With `-t`, `tandem/` holds the genome databases (`.1gdb`), the tandem arrays found by FasTAN (`<genome>_tan.1ano`) and the ones masked, at least `-w` bp long (`<genome>_tan<w>.1ano`), when kept with `-k` or `-d`. The alignments refer to the databases, so keep them together.

Elements are named `<class>_<n>-<lib_name>#<class>`, numbered by length within each class (accepted and rejected elements share the numbering).

## Tandem masking

Long tandem arrays (satellites, centromeres) make FastGA extend huge alignments through them, and can dominate the alignment time. With `-t`, FasTAN finds the tandem arrays of each genome and those at least `-w` bp long (1 kb by default) are soft masked in the genome index, so FastGA starts no alignment inside them (alignments can still extend through them). Coordinates and sequences are unchanged: segments, alignments and annotations refer to the original genomes.

Short arrays are left unmasked because they add little alignment time, but are found inside some TEs: masking them all lost, in a test, a young Gypsy element with 51 bp tandem repeats in its LTRs.

In a test with the two haplotypes of gorilla chromosome 22 (T2T, 59 Mb, 36% in tandem arrays, 2 threads), FastGA took 632 s without masking, 20 s masking all the arrays (8,952 in one haplotype) and 21 s masking those of 1 kb or more (363, holding 95% of the arrays' length). svfind found 201 segments without masking; masking all arrays kept 176 of them, and masking arrays of 2 kb or more 186. The segments lost were at tandem arrays (likely array length differences). In a 25 Mb genome with 6% in tandem arrays, masking made no difference in time.

## Citation

SOON?

The classification model: doi:[10.5281/zenodo.22990589](https://doi.org/10.5281/zenodo.22990589)

## Acknowledgements

Thanks to Arian Smit and Robert Hubley (Dfam, RepeatMasker, RepeatModeler) for allowing us to use their curated peptide library in this release.
Thanks to Claude for helping with the development.

## Licence

MIT, see [LICENSE](LICENSE).













