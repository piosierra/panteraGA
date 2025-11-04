![pantera](images/panteraGA.png?raw=true "pantera")
# Identification of transposable element families from whole genome alignments using FastGA

This new version of [**pantera**](https://github.com/piosierra/pantera) extracts the polymorphisms found in a whole genome alignment 
and uses them to build a TE library. It greatly benefits from the speed of the blazing fast whole genome aligner [**FastGA**](https://github.com/thegenemyers/FASTGA) to generate TE libraries in just minutes for most species, and can handle large, over 10Gb genomes, in a few hours.

### 0- Requirements

In addition to FastGA to build the alignments, **pantera** requires that several utilities are available on the path:

[alntools](https://github.com/richarddurbin/alntools)

[mafft](https://mafft.cbrc.jp/alignment/software/) 

[blast](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata)

[cd-hit](https://github.com/weizhongli/cdhit) 

[emboss](https://emboss.sourceforge.net/download/#Stable/) 

**pantera** has been tested in Linux with R 4.5.1

**pantera** needs several R packages in your system, if they are not found they will be installed during the first run, which will make it slower, subsequent runs will not require that and will be faster.


### 1- Installing
Simply download the last release and run panteraGA to install any missing R libraries.
```
./panteraGA
```

### 1- Aligning two genomes
```
FastGA -1:alignment_file genome1.fa genome2.fa
```

### 2- Obtain the library from the genome alignment
```
./panteraGA -1 alignment_file 
```
For more options just use:
```
./panteraGA -h
```

### 3 Content of the output folder

- `Myspecies-pantera-pass.fa` Final filtered library (recomended to use this one)
- `Myspecies-pantera-final.fa` Final unfiltered library.
- `Myspecies-pantera-final.stats` Stats for the families found.
- `Myspecies-pantera.log` Runtime log.

### Akwnoledgements
Thanks to Arian Smit and Robert Hubley (Dfam, RepeatMasker, RepeatModeler) for allowing us make use of their curated peptides library in this release.















