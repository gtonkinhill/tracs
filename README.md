# Mtran

Mtran provides robust estimates of pairwise transmission distances from single isolate, multi-strain and metagenomic samples. It uses an empirical Bayes approach to account for variable sequence coverage and aligns to multiple reference genomes to estimate a lower bound for both the SNP distance and the number of intermediate hosts seperating two samples.

## Installation

### Conda

Coming soon...

### Manual

Mtran is a python package and can be installed easily using pip. 

```
pip3 install git+https://github.com/gtonkinhill/mtran
```

This is all that is needed for the pairwise distance command. To run the alignment command you will also need to install [minimap2](https://github.com/lh3/minimap2) and [htsbox](https://github.com/lh3/htsbox).

## Alignment

```
Uses sourmash to identify reference genomes within a read set and then aligns
reads to each reference using minimap2

options:
  -h, --help            show this help message and exit
  -t N_CPU, --threads N_CPU
                        number of threads to use (default=1)
  --quiet               turns off some console output
  --version             show program's version number and exit

Input/output:
  -i INPUT_FILES [INPUT_FILES ...], --input INPUT_FILES [INPUT_FILES ...]
                        path to query signature
  --database DATABASE   path to database signatures
  -o OUTPUT_DIR, --output OUTPUT_DIR
                        location of an output directory
  -p PREFIX, --prefix PREFIX
                        prefix to describe the input sample read files

Alignment options:
  --minimap_preset MINIMAP_PRESET
                        minimap preset to use - one of 'sr' (default), 'map-
                        ont' or 'map-pb'

Pileup options:
  -Q MIN_BASE_QUAL, --min_base_qual MIN_BASE_QUAL
                        minimum base quality (default=0)
  -q MIN_MAP_QUAL, --min_map_qual MIN_MAP_QUAL
                        minimum mapping quality (default=0)
  -l MIN_QUERY_LEN, --min_query_len MIN_QUERY_LEN
                        minimum query length (default=0)

Posterior count estimates:
  --threshold EXPECTED_FREQ_THRESHOLD
                        Minimum posterior read frequency threshold. The
                        default is set that a variant at a location is
                        discounted if it is not found with a coverage of ~100x
  --both-strands        turns on the requirement that a variant is supported
                        by both strands
  --filter-all          turns on filtering of variants with support below the
                        posterior frequency threshold
```

## Pairwise Transmission Distance

usage: dist [-h] --msa MSA_FILES [MSA_FILES ...] --dates METADATA -o
            OUTPUT_FILE [-D SNP_THRESHOLD] [--clock_rate CLOCK_RATE]
            [--trans_rate TRANS_RATE] [-K TRANS_THRESHOLD] [-t N_CPU]
            [--quiet] [--version]

Estimates pairwise SNP and transmission distances between each pair of samples
aligned to the same reference genome.

options:
  -h, --help            show this help message and exit
  -t N_CPU, --threads N_CPU
                        number of threads to use (default=1)
  --quiet               turns off some console output
  --version             show program's version number and exit

Input/output:
  --msa MSA_FILES [MSA_FILES ...]
                        Input fasta files formatted by the align and merge
                        functions
  --dates METADATA      Location of metadata in csv format. The first column
                        must include the sequence names and the second column
                        must include sampling dates.
  -o OUTPUT_FILE, --output OUTPUT_FILE
                        name of the output file to store the pairwise distance
                        estimates.

SNP distance options:
  -D SNP_THRESHOLD, --snp_threshold SNP_THRESHOLD
                        Only output those transmission pairs with a SNP
                        distance <= D

Transmission distance options:
  --clock_rate CLOCK_RATE
                        clock rate as defined in the transcluster paper
                        (SNPs/genome/year) default=1e-3 * 29903
  --trans_rate TRANS_RATE
                        transmission rate as defined in the transcluster paper
                        (transmissions/year) default=73
  -K TRANS_THRESHOLD, --trans_threshold TRANS_THRESHOLD
                        Only outputs those pairs where the most likely number
                        of intermediate hosts <= K

## Clustering

Coming soon...

