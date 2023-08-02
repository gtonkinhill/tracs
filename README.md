# Trac'm

<p align="center">
<img src="https://github.com/gtonkinhill/tracm/blob/main/docs/_figures/tracm_logo.png" alt="alt text" width="500">
</p>

Trac'm provides robust estimates of pairwise transmission distances from single isolate, multi-strain and metagenomic samples. It uses an empirical Bayes approach to account for variable sequence coverage and aligns to multiple reference genomes to estimate a lower bound for both the SNP distance and the number of intermediate hosts seperating two samples.

## Documentation

Trac'm is currently under development and frequent backwards incompatable changes may be made.

Documentation for Trac'm can be found [here](https://gtonkinhill.github.io/tracm)

## Installation

### Conda

A proper Conda install is coming soon.

In the mean time, Trac'm can be installed using conda by running

```
conda create -n tracm python=3.10
conda activate tracm
conda install -c bioconda samtools htsbox minimap2
pip3 install git+https://github.com/gtonkinhill/tracm
```

### Manual

Tracm is a python package and can be installed easily using pip. 

```
pip3 install git+https://github.com/gtonkinhill/tracm
```

This is all that is needed for the pairwise distance and clustering commands. To generate alignments you will also need to install [samtools](http://www.htslib.org/), [minimap2](https://github.com/lh3/minimap2) and [htsbox](https://github.com/lh3/htsbox).
