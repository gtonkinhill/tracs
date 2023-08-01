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

All files needed to run these examples are provided [here](github.com).

### Isolate

Isolate data can either be generated using the align command or using alternative pipelines such as [Snippy](https://github.com/tseemann/snippy). 

#### SARS-CoV-2

We consider a multiple sequence alignment of SARS-CoV-2 genomes from Tonkin-Hill, Martinconera et al., *Elife* 2021. Here, we are assuming a transmission generation time of 5 days (5/356 = 73) and a clock rate of 1e-3 per base per year (1e-3 * 29903 = 29.03). The sample dates are proved as a csv formatted file with one sample per line. 

As we do not expect much recombination within these sample we do not use the `--filter` parameter.

```
tracm distance --msa MA_combined_consensus_replicates_filt_dates.fa  --meta combined_consensus_replicates_filt_dates.csv -o test.out --trans_rate 73 --clock_rate 29.03
```

### Metagenomic



### Multi-strain