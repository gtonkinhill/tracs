# Trac'm

Trac'm provides robust estimates of pairwise transmission distances from single isolate, multi-strain and metagenomic samples. It uses an empirical Bayes approach to account for variable sequence coverage and aligns to multiple reference genomes to estimate a lower bound for both the SNP distance and the number of intermediate hosts separating two samples.

## Getting started

The Trac'm package consists of several modules that can be used independently. 

<center><img src="_figures/tracm_flow.drawio.svg" width="700"></center>

If you already have a Multiple Sequence Alignment (MSA) you can jump straight to the [distance](distance.md) command.

To start from raw single isolate or metagenomic sequencing data, you need to generate alignments by running the [align](alignment.md) command on each sample separately.

These can then be combined into multiple sequence alignments using the [combine](combine.md) command before running the [distance](distance.md) function.

Finally, putative transmission clusters can be generated using single linkage hierarchical clustering via the [cluster](cluster.md) command.

When dealing with a small number of samples, the [pipe](pipe.md) command can be used to automate the running of each step of the Trac'm pipeline.

## Examples

### Isolate

### Metagenomic

### Multi-strain