# Trac'm

Trac'm provides robust estimates of pairwise transmission distances from single isolate, multi-strain and metagenomic samples. It uses an empirical Bayes approach to account for variable sequence coverage and aligns to multiple reference genomes to estimate a lower bound for both the SNP distance and the number of intermediate hosts separating two samples.

## Getting started

The Trac'm package consists of several modules that can be used independently. 

<center><img src="../_figures/tracm_flow.drawio.svg" width="700"></center>

If you already have a Multiple Sequence Alignment (MSA) you can jump straight to the [distance](distance.md) command.

To start from raw single isolate or metagenomic sequencing data, you need to generate alignments by running the [align](alignment.md) command on each sample separately.

These can then be combined into multiple sequence alignments using the [combine](combine.md) command before running the [distance](distance.md) function.

Finally, putative transmission clusters can be generated using single linkage hierarchical clustering via the [cluster](cluster.md) command.

When dealing with a small number of samples, the [pipe](pipe.md) command can be used to automate the running of each step of the Trac'm pipeline.

## Examples

All files needed to run these examples are provided [here](https://zenodo.org/record/8202050).

### Isolate

Isolate data can either be generated using the [align](alignment.md) command or using alternative pipelines such as [Snippy](https://github.com/tseemann/snippy). 

#### SARS-CoV-2

Let's start by considering a multiple sequence alignment of SARS-CoV-2 genomes from Tonkin-Hill, Martinconera et al., *Elife* 2021. 

To calculate SNP distances we can run Trac'm using the following command:

```
tracm distance --msa MA_combined_consensus_replicates_filt_dates.fa  --meta combined_consensus_replicates_filt_dates.csv -o transmission_distances.csv
```

In order to estimate the the transmission distance or number of intermediate hosts, we need to provide a transmission generation time and clock rate.

Here, we assume a transmission generation time of 5 days (5/356 = 73) and a clock rate of 1e-3 per base per year (1e-3 * 29903 = 29.03). The sample dates are proved as a csv formatted file with one sample per line

```
tracm distance --msa MA_combined_consensus_replicates_filt_dates.fa  --meta combined_consensus_replicates_filt_dates.csv -o transmission_distances.csv --trans_rate 73 --clock_rate 29.03
```

#### Cluster

To cluster the inferred distances we can run the [cluster](cluster.md) command.

Clustering can be performed using a number of different distance metrics including SNP and the estimated number of intermediate hosts according to the Trac'm model. In the case we use the latter and choose a threshold of 5 intermediate hosts to separate clusters.

```
tracm cluster -d transmission_distances.csv -D expectedK -c 5 -o clusters.csv
```

### Metagenomic

For the matagenomics example, we consider a pair of simulated human gut microbiome samples. The simulated proportions for each species and the respective SNP distances is provided in the `simulated_proportions.csv` file.

#### Align

We first need to align the sequencing reads to reference genomes. The set of references can either be a custom database generated using the [build-db](database.md) command or a GTDB/sourmash database. 

If only the Sourmash database is supplied, Trac'm automatically downloads the reference genomes corresponding to the species observed within the given sample. This eliminates the inconvenience associated with downloading the comprehensive reference genome database from GTDB. However, when analysing many samples, it can be more advantageous to preemptively download the complete set of genomes from GTDB. 

Run the [align](alignment.md) command separately on each sample

```
tracm align -i sim_d5_ref_GCF_018292165.1_ASM1829216v1_genomic_A*.fastq.gz -o sampleA --prefix sampleA --keep-all -t 20 --database gtdb-rs214-reps.k51.sbt.zip
```

```
tracm align -i sim_d5_ref_GCF_018292165.1_ASM1829216v1_genomic_B*.fastq.gz -o sampleB --prefix sampleB --keep-all -t 20 --database gtdb-rs214-reps.k51.sbt.zip
```

#### Combine

We can now combine the alignments by species into multiple sequence alignments (MSAs).

```
tracm combine -i sampleA sampleB -o combined_alignments
```

#### Distance

Finally, we can estimate SNP distances across all species common between the two samples. We use the `--filter` command to account for shared homology between genomes and other sources of noise

```
tracm distance --msa ./combined_alignments/*.fasta.gz -o transmission_distances.csv --filter
```


### Multi-strain

Here we consider two labratory mixtures of pneumococcal strains originally published in Knight et al., *Microbial Genomics* 2021. The first contains strains of serotype 1 and 19F. The second contains the same two serotypes in addition to serotypes 4 and 18C. 

We use a custom database of pneumococcal reference genomes provided as part of this example (pneumoexampleDB.zip) which was built using the [build-db](database.md) command.

#### Align

First we align the samples in a similar way to the metagenomics example.

```
tracm align -i subset_SRR9998185_WGS_of_strep_pneumoniae_Serotype_mixture_1_19F_*.fastq.gz -o sampleA --prefix sampleA --keep-all -t 20 --database pneumoexampleDB.zip
```

```
tracm align -i subset_SRR9998201_WGS_of_strep_pneumoniae_Serotype_mixture_1_4_18C_19F_*.fastq.gz -o sampleB --prefix sampleB --keep-all -t 20 --database pneumoexampleDB.zip
```

#### Combine

We can now combine the results into individual MSAs for each strain that is shared between the two samples

```
tracm combine -i sample* -o combined_alignments
```

#### Distance

Finally, we can calculate the transmission distance between the isolates. Here, we are assuming a transmission generation time of 2 months (2/12 = 6) and a clock rate of 5.3 mutations/genome/year. An artificial sample date is provided to indicate that these samples were taken on the same day

```
tracm distance --msa ./combined_alignments/*.fasta.gz -o transmission_distances.csv --trans_rate 6 --clock_rate 5.3 --meta dates.csv --filter
```
