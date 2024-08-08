# Installation

## Conda

A proper Conda install is coming soon.

In the meantime, TRACS can be installed using conda by running

```
conda create -n tracs python=3.10
conda activate tracs
conda install -c bioconda samtools htsbox minimap2
pip3 install git+https://github.com/gtonkinhill/tracs
```

## Manual

TRACS is a python package and can be installed easily using pip. 

```
pip3 install git+https://github.com/gtonkinhill/tracs
```

This is all that is needed for the pairwise distance and clustering commands. To generate alignments you will also need to install [minimap2](https://github.com/lh3/minimap2) and [htsbox](https://github.com/lh3/htsbox).

