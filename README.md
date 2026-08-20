# TRACS

[![tracs-CI](https://github.com/gtonkinhill/tracs/actions/workflows/tracs_test.yml/badge.svg)](https://github.com/gtonkinhill/tracs/actions/workflows/tracs_test.yml)

<p align="center">
<img src="https://github.com/gtonkinhill/tracs/blob/main/docs/_figures/tracs_logo.png" alt="alt text" width="500">
</p>

TRACS provides robust estimates of pairwise transmission distances from single isolate, multi-strain and metagenomic samples. It uses an empirical Bayes approach to account for variable sequence coverage and aligns to multiple reference genomes to estimate a **lower bound** for both the SNP distance and the number of intermediate hosts separating two samples.

**Note: TRACS is not intended to estimate very large SNP distances**

## Documentation

TRACS is currently under development and frequent backwards incompatible changes may be made.

Documentation for TRACS can be found [here](https://gtonkinhill.github.io/tracs)

## Installation

### Conda

TRACS can be installed using conda by running

```
conda install bioconda::tracs
```

**Note:** Conda and container builds are not guaranteed to be portable across all CPU generations due to microarchitecture differences. If you run into a `SIGILL` or `Illegal instruction (core dumped)` error when using the bioconda package or biocontainers image, it is likely because your CPU does not support the instructions the binary was built with. Building from source using `pip` can solve this problem.

### Manual

TRACS is a python package and can be installed easily using pip. 

```
pip3 install git+https://github.com/gtonkinhill/tracs
```

By default, building from source targets a baseline architectural compatibility (`x86-64-v2`) for portability. If you want to optimize the build for your current CPU to obtain maximum performance, you can set the `TRACS_MARCH` environment variable before installing:

```
TRACS_MARCH=native pip3 install git+https://github.com/gtonkinhill/tracs
```

This is all that is needed for the pairwise distance and clustering commands. To generate alignments you will also need to install [samtools](http://www.htslib.org/), [minimap2](https://github.com/lh3/minimap2) and [htsbox](https://github.com/lh3/htsbox).
