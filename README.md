# TRACS

[![tracs-CI](https://github.com/gtonkinhill/tracs/actions/workflows/tracs_test.yml/badge.svg)](https://github.com/gtonkinhill/tracs/actions/workflows/tracs_test.yml)

<p align="center">
<img src="https://github.com/gtonkinhill/tracs/blob/main/docs/_figures/tracs_logo.png" alt="alt text" width="500">
</p>

TRACS provides robust estimates of pairwise transmission distances from single isolate, multi-strain and metagenomic samples. It uses an empirical Bayes approach to account for variable sequence coverage and aligns to multiple reference genomes to estimate a **lower bound** for both the SNP distance and the number of intermediate hosts seperating two samples.

**Note: TRACS is not intended to estimate very large SNP distances**

## Documentation

TRACS is currently under development and frequent backwards incompatable changes may be made.

Documentation for TRACS can be found [here](https://gtonkinhill.github.io/tracs)

## Installation

### Conda

A proper Conda install is coming soon.

In the mean time, TRACS can be installed using conda by running

```
conda create -n tracs python=3.10
conda activate tracs
conda install -c bioconda samtools htsbox minimap2
pip3 install git+https://github.com/gtonkinhill/tracs
```

### Manual

TRACS is a python package and can be installed easily using pip. 

```
pip3 install git+https://github.com/gtonkinhill/tracs
```

This is all that is needed for the pairwise distance and clustering commands. To generate alignments you will also need to install [samtools](http://www.htslib.org/), [minimap2](https://github.com/lh3/minimap2) and [htsbox](https://github.com/lh3/htsbox).
