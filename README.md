![pantera](images/panteraGA.png?raw=true "pantera")
# Identification of transposable element families from whole genome alignments using FastGA

[**IMPORTANT**]
Only the full releases include the necesary files. If you just clone the repository the model for classification will be missing and you will receive an error.

This new version of [**pantera**](https://github.com/piosierra/pantera) aligns one or more genomes with [**FastGA**](https://github.com/thegenemyers/FASTGA) to generate TE libraries in just minutes for most species, and can handle large, over 10Gb genomes, in a few hours. Rather than performing an all vs all alignment, in the case of more than two genomes, pairwise aligments will be generated on the list of genomes, and all the polymorphic segments processed together afterwards. 

### The easy way
Build a singularity image from the Docker container:

`singularity pull panteraGA.sif docker://piosierra/pantera_ga:latest`

Run pantera:

`singularity exec panteraGA.sif panteraGA -g listofgenomes -b identifier -o outputfolder`

- `listofgenomes` is a plain text file with the paths to the genomes you want to use to build the library (min 2)
- `identifier` is the name to append to element names in the library



### Akwnoledgements
Thanks to Arian Smit and Robert Hubley (Dfam, RepeatMasker, RepeatModeler) for allowing us make use of their curated peptides library in this release.















