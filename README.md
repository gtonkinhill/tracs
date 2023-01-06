# Tracm

Trac'm provides robust estimates of pairwise transmission distances from single isolate, multi-strain and metagenomic samples. It uses an empirical Bayes approach to account for variable sequence coverage and aligns to multiple reference genomes to estimate a lower bound for both the SNP distance and the number of intermediate hosts seperating two samples.

## Installation

### Conda

Coming soon...

### Manual

Tracm is a python package and can be installed easily using pip. 

```
pip3 install git+https://github.com/gtonkinhill/tracm
```

This is all that is needed for the pairwise distance and clustering commands. To generate alignments you will also need to install [minimap2](https://github.com/lh3/minimap2) and [htsbox](https://github.com/lh3/htsbox).
