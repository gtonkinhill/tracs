# Alignment

The **align** module takes raw sequence data as input and produces alignments to multiple reference genomes in FASTA format. IUPAC ambiguity codes are used to represent polymorphic sites in the alignment. This module should be run on each sample separately.

If a database containing multiple references is provided, the algorithm will first select those represented in the sample using [sourmash](https://sourmash.readthedocs.io/en/latest/index.html). 

It will then align all reads to each reference in turn and report the alleles present taking into account sequencing depth.

The major options that should be considered are 

|     **Column**    |                           **Description**                           |
|:-----------------:|:-------------------------------------------------------------------:|
|    `--min-cov`    |     Specifies the minimum sequencing coverage to consider a site    |
|   `--error-perc`  |          Specifies the minimum frequency to retain variants (will be overridden at the reporting stage by `--keep-all`)         |
| `--either-strand` |       Allows for variants that are present on only one strand       |
|    `--keep-all`   | Whether to keep all observed variants regardless of their frequency |


### Example

#### Single reference 

To align to a single reference genome the `--refseqs` option can be set as

```
tracm align -i read_file_1.fq read_file_2.fq --refseqs reference_genome.fasta -o output_folder -p sample_name
```

#### Metagenomic data with GTDB database

To analyse a metagenomic dataset using the GTDB database as a reference you can either have Trac'm automatically download the required genomes using the precombiled sourmash databases availble [here](https://sourmash.readthedocs.io/en/latest/databases.html)

Trac'm expects the SBT format files ending in `.sbt.zip`

```
tracm align -i read_file_1.fq read_file_2.fq --database gtdb-rs207.genomic-reps.dna.k51.sbt.zip -o output_folder -p sample_name
```

Alternatively, you can also specify the location of the predownloaded GTDB fasta files. This avoids redownloading them which is particularly useful when analysing large numbers of samples on a compute cluster. These folders can be downloaded from [here](https://data.gtdb.ecogenomic.org/). You will need the folder containing all representative genomes, usually labelled something like `gtdb_genomes_reps_r207.tar.gz`. This then needs to be unzipped prior to running Trac'm. It is the users responsibility to ensure the correct pairing of sourmash and GTDB reference folders is used.

```
tracm align -i read_file_1.fq read_file_2.fq --database gtdb-rs207.genomic-reps.dna.k51.sbt.zip --refseqs gtdb_genomes_reps_r207 -o output_folder -p sample_name
```

#### Species level and custom databases

Finally, it possible to run Trac'm using a custom database generated using the [database]() command. Pre-compiled version of the for S. pneumoniae is available [here]().

```
tracm align -i read_file_1.fq read_file_2.fq --database reference_db.zip -o output_folder -p sample_name
```

### Options

```
usage: tracm align [-h] -i INPUT_FILES [INPUT_FILES ...] [--database DATABASE]
                   [--refseqs REFSEQS] -o OUTPUT_DIR [-p PREFIX]
                   [--minimap_preset MINIMAP_PRESET] [-Q MIN_BASE_QUAL]
                   [-q MIN_MAP_QUAL] [-l MIN_QUERY_LEN] [-V MAX_DIV]
                   [--trim TRIM] [--min-cov MIN_COV]
                   [--error-perc ERROR_THRESHOLD] [--either-strand]
                   [--keep-all] [-t N_CPU]
                   [--loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

Uses sourmash to identify reference genomes within a read set and then aligns
reads to each reference using minimap2

options:
  -h, --help            show this help message and exit
  -t N_CPU, --threads N_CPU
                        number of threads to use (default=1)
  --loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Set the logging threshold.

Input/output:
  -i INPUT_FILES [INPUT_FILES ...], --input INPUT_FILES [INPUT_FILES ...]
                        path to query signature
  --database DATABASE   path to database signatures
  --refseqs REFSEQS     path to reference fasta files
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
  -V MAX_DIV, --max_div MAX_DIV
                        ignore queries with per-base divergence > max_div
                        (default=1)
  --trim TRIM           ignore bases within TRIM-bp from either end of a read
                        (default=0)

Posterior count estimates:
  --min-cov MIN_COV     Minimum read coverage (default=5).
  --error-perc ERROR_THRESHOLD
                        Threshold to exclude likely erroneous variants
  --either-strand       turns off the requirement that a variant is supported
                        by both strands
  --keep-all            turns on filtering of variants with support below the
                        posterior frequency threshold
```