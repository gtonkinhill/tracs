# Installation

## Conda

A proper install is coming soon.

In the mean time, Trac'm can be installed using conda by running

```
conda create -n tracm python=3.10
conda activate tracm
conda install -c bioconda samtools htsbox minimap2
pip3 install git+https://github.com/gtonkinhill/tracm
```

## Manual

Tracm is a python package and can be installed easily using pip. 

```
pip3 install git+https://github.com/gtonkinhill/tracm
```

This is all that is needed for the pairwise distance and clustering commands. To generate alignments you will also need to install [minimap2](https://github.com/lh3/minimap2) and [htsbox](https://github.com/lh3/htsbox).

