# Pipe

The pipe command can be used to run the full TRACS pipeline on a set of samples. It is most useful when considering a small number of samples. It is usually better to run large datasets in stages, which makes it easier to run multiple `align` commands in parallel.

The `pipe` command takes a text file as input with the sample name followed by the sequencing data in the subsequent columns. These can either be paired-end read files, a single combined read file, or an assembly. An example of the format is given below.

```
prefix  read1 read2
sampleA sampleA_R1.fastq.gz sampleA_R2.fastq.gz
sampleB sampleA_combined.fastq.gz
sampleC sampleC.fasta
```


### Example

```
tracs pipe -i input_file_list.txt -o output_directory -t 15 --database gtdb-rs207.genomic-reps.dna.k51.sbt.zip --keep-all --filter --cluster_distance filter -c 10
```

### Options

```
usage: tracs pipe [-h] -i INPUT_FILE --database DATABASE -o OUTPUT_DIR
                      [--meta METADATA] [--minimap_preset MINIMAP_PRESET]
                      [-Q MIN_BASE_QUAL] [-q MIN_MAP_QUAL] [-l MIN_QUERY_LEN]
                      [-V MAX_DIV] [--trim TRIM] [--min-cov MIN_COV]
                      [--error-perc ERROR_THRESHOLD] [--either-strand]
                      [--keep-all] [-D SNP_THRESHOLD] [--filter]
                      [--clock_rate CLOCK_RATE] [--trans_rate TRANS_RATE]
                      [-K TRANS_THRESHOLD] [--precision PRECISION]
                      [-c THRESHOLD]
                      [--cluster_distance {SNP,direct,expectedK}] [-t N_CPU]
                      [--loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

A script to run the full TRACS pipeline.

options:
  -h, --help            show this help message and exit
  -t N_CPU, --threads N_CPU
                        number of threads to use (default=1)
  --loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Set the logging threshold.

Input/output:
  -i INPUT_FILE, --input INPUT_FILE
                        path to text file containing input file paths
  --database DATABASE   path to database signatures
  -o OUTPUT_DIR, --output OUTPUT_DIR
                        location of an output directory
  --meta METADATA       Location of metadata in csv format. The first column
                        must include the sequence names and the second column
                        must include sampling dates.

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
                        Threshold to exclude likely erroneous variants prior
                        to fitting Dirichlet multinomial model
  --either-strand       turns off the requirement that a variant is supported
                        by both strands
  --keep-all            turns on filtering of variants with support below the
                        posterior frequency threshold

SNP distance options:
  -D SNP_THRESHOLD, --snp_threshold SNP_THRESHOLD
                        Only output those transmission pairs with a SNP
                        distance <= D
  --filter              Filter out regions with unusually high SNP distances
                        often caused by HGT

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
  --precision PRECISION
                        The precision used to calculate E(K) (default=0.01).

Cluster options:
  -c THRESHOLD, --cluster_threshold THRESHOLD
                        Distance threshold. Samples will be grouped together
                        if the distance between them is below this threshold.
                        (default=10)
  --cluster_distance {snp,filter,direct,expectedK}
                        The type of transmission distance to use. Can be one
                        of 'snp', 'filter', 'direct', 'expectedK'
```