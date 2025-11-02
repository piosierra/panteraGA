![pantera](images/panteraGA.png?raw=true "pantera")
# Identification of transposable element families from whole genome alignments using FastGA

This new version of [**pantera**](https://github.com/piosierra/pantera) extracts the polymorphisms found in a whole genome alignment 
to build with them a TE library. It benefits of the amazing speed of [**FastGA**](https://github.com/thegenemyers/FASTGA) to generate TE libraries in just minutes for most species, and can handle large, over 10Gb genomes, in a few hours.

### 0- Requirements

**pantera** requires that several utilities are available on the path:

[ONEcode](https://github.com/thegenemyers/ONEcode)

[alntools](https://github.com/richarddurbin/alntools)

[mafft](https://mafft.cbrc.jp/alignment/software/) 

[blast](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata)

[cd-hit](https://github.com/weizhongli/cdhit) 

[emboss](https://emboss.sourceforge.net/download/#Stable/) 

**pantera** has been tested in Linux with R 4.5.1

**pantera** needs several R packages in your system, if they are not found the first time it runs it will be slower as it will proceed to install them.


### 1- Installing
Simply clone the repository. Then configure **panteraGA** with the path to your RepeatMasker installation:
```
./panteraGA install path_to_RepeatMasker_folder
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
```
Myspecies-pantera-final.fa  Myspecies-pantera-final.stats  LmiMyspeciesc-pantera-pass.fa  pantera.log
```
- `pantera-pass.fa` Final filtered library (recomended to use this one)
- `pantera-final.fa` Final unfiltered library.
- `pantera-final.stats` Stats for the families found.
- `pantera.log` Runtime log.













