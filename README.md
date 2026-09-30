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

This installs panteraGA with all its dependencies (FastGA, alntools, BLAST, MAFFT, CD-HIT, EMBOSS and the R packages) and the classification model.

### Containers

Every bioconda package is also available as a container from [BioContainers](https://biocontainers.pro). Find the current tag at [quay.io/repository/biocontainers/panteraga](https://quay.io/repository/biocontainers/panteraga?tab=tags), then:

```bash
singularity pull panteraGA.sif docker://quay.io/biocontainers/panteraga:<tag>
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
| `-g, --genomes` | required | File with the list of genomes |
| `-b, --lib_name` | required | Identifier appended to element names |
| `-o, --output_folder` | `pantera_output` | Output folder |
| `-T, --threads` | 8 (or the number of cores, if fewer) | Threads |
| `-s, --min_size` | 100 | Min length of polymorphic segments |
| `-l, --max_size` | 30000 | Max length of polymorphic segments |
| `-i, --identity` | 0.90 | Identity for clustering, first round |
| `-y, --identity2` | 0.85 | Identity for clustering, second round |
| `-m, --min_cl` | 3 | Min sequences to build a consensus |
| `-e, --mingen` | 2 | Min copies in at least one genome for an element to pass |
| `-n, --Ns` | 0.001 | Max fraction of Ns in a segment |
| `-c, --cons_Ns` | 0.02 | Max fraction of Ns in a final consensus |
| `-p, --pAs` | 10 | Min length of a polyA tail |
| `-u, --cl_size` | 300 | Max sequences per cluster |
| `-f, --flanking` | 100 (automatic) | Min length of the flanks of a segment; 100 uses `-q` instead |
| `-q, --flank_quantile` | 0.05 | Fraction of segments with the shortest flanks to discard |
| `-a, --anno_per` | 0.80 | Min coverage of an element for genome annotation |
| `-z, --anno_div` | 0.80 | Min identity for genome annotation |
| `-k, --keep` | | Keep the FastGA alignments (`.1aln`) |
| `-d, --debug` | | Keep intermediate files |
| `-v, --verbose` | | Show log messages on screen |
| `-h, --help` | | Show the options |

### Output

In the output folder:

| File | Content |
|---|---|
| `mylib-pantera-final.fa` | All TE consensus sequences |
| `mylib-pantera-final-pass.fa` | Only the elements that pass the structural checks |
| `mylib-pantera-final.stats.tsv` | Per element: length, class, pass, cluster size, TSD, terminal repeats, tails, copies per genome... |
| `alignments/` | The alignment (`.maf`) behind each consensus, named as the element |
| `annotations/` | A BED file per genome with the copies of each element |
| `pantera.log` | Log of the run |

Elements are named `<class>_<n>-<lib_name>#<class>`, numbered by length within each class.

## Citation

<!-- SOON? -->

The classification model: doi:[10.5281/zenodo.22990589](https://doi.org/10.5281/zenodo.22990589)

## Acknowledgements

Thanks to Arian Smit and Robert Hubley (Dfam, RepeatMasker, RepeatModeler) for allowing us to use their curated peptide library in this release.

## Licence

MIT, see [LICENSE](LICENSE).













