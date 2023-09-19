# Transmission distance

The **distance** command rapidly calculates pairwise SNP and transmission distance estimates. It accepts a FASTA formatted multiple sequence alignment as input and is aware of IUPAC ambiguity codes. Transmission distance is estimated using a modified version of the [TransCluster](https://doi.org/10.1093/molbev/msy242) algorithm and will only be calculated if both the clock and transmission rate parameters are set.

Trac'm includes a filtering strategy to account for problematic regions in alignments. These are usually the result of shared homology between species and strains within mixed samples or by hard-to-align regions of the genome. This can be enabled using the `--filter` option.

It is recommended that a maximum SNP distance is set to speed up the algorithm and reduce the size of the output file. This is calculated before any filtering is done and thus it is recommended that it be an order of magnitude larger than the minimum SNP distance of interest.

## Example

### Simple SNP distance

To run the **distance** command on an MSA file called "msa.fasta" with filtering and a minimum SNP distance of 1000 you can use the following command.

```
tracm distance -i msa.fasta -o transmission_distances.csv --snp_threshold 1000 --filter
```

### Transmission distance (TransCluster)

Using estimates of the transmission and clock rate of SARS-CoV-2 (in transmissions/year and SNPs/genome/year respectively), the distance pipeline can be run as

```
tracm distance -i msa.fasta --clock_rate 29.03 --trans_rate 73 --meta dates.csv
```

Here, we are assuming a transmission generation time of 5 days (5/356 = 73) and a clock rate of 1e-3 per base per year (1e-3 * 29903 = 29.03). The sample dates are proved as a csv formatted file with one sample per line for example

```
sample,date
seq1,2021-01-01
seq2,2021-12-31
```

## Output

The **distance** command generates a csv with the following columns

|       **Column**      |                                          **Description**                                          |
|:---------------------:|:-------------------------------------------------------------------------------------------------:|
|        sampleA        |                                            sample name                                            |
|        sampleB        |                                            sample name                                            |
|    date difference    |                           time between samples (if `--meta` is provided)                          |
|      SNP distance     |                    classic SNP distance (accounting for IUPAC ambiguity codes)                    |
| transmission distance |                 the probability of a direct transmission (if `--clock_rate` and `--trans_rate` are provided)                |
|       expected K      |                the expected number of intermediate hosts (if `--clock_rate` and `--trans_rate` are provided)                |
| filtered SNP distance |                   SNP distance after removing regions with unusual SNP patterns                   |
|    sites considered   | the number of informative sites (not ambiguous Ns) - can be used to calculate percentage identity |

## Options

```
usage: tracm distance [-h] --msa MSA_FILES [MSA_FILES ...] [--msa-db MSA_DB]
                      [--meta METADATA] -o OUTPUT_FILE [-D SNP_THRESHOLD]
                      [--filter] [--clock_rate CLOCK_RATE]
                      [--trans_rate TRANS_RATE] [-K TRANS_THRESHOLD]
                      [--precision PRECISION] [-t N_CPU]
                      [--loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

Estimates pairwise SNP and transmission distances between each pair of samples
aligned to the same reference genome.

options:
  -h, --help            show this help message and exit
  -t N_CPU, --threads N_CPU
                        number of threads to use (default=1)
  --loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Set the logging threshold.

Input/output:
  --msa MSA_FILES [MSA_FILES ...]
                        Input fasta files formatted by the align and merge
                        functions
  --msa-db MSA_DB       A database MSA used to compare each sequence to. By
                        default this is not uses and all pairwise comparisons
                        within each MSA are considered.
  --meta METADATA       Location of metadata in csv format. The first column
                        must include the sequence names and the second column
                        must include sampling dates.
  -o OUTPUT_FILE, --output OUTPUT_FILE
                        name of the output file to store the pairwise distance
                        estimates.

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
```
