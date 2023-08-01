# Databases

## Prebuilt GTDB database

Trac'm can accept prebuilt databases provided on the GTDB and Sourmash websites. If only the Sourmash database is supplied, Trac'm automatically downloads the reference genomes corresponding to the species observed within the given sample from RefSeq. This eliminates the inconvenience associated with downloading the comprehensive reference genome database from GTDB which can be very large. However, when analysing multiple samples, it can be more advantageous to preemptively download the complete set of genomes from GTDB.

The Sourmash file can be provided to Trac'm as is. When using the GTDB references it is necessary to untar the file first by running `tar -xf`


| Database 	| Sourmash DB 	| GTDB reference genomes 	|
|:---:	|:---:	|:---:	|
| GTDB R08-RS214 genomic representatives (85k) 	| [gtdb-rs214-reps.k51.sbt.zip](https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k51.sbt.zip) (4.4Gb) 	| [gtdb_genomes_reps_r214.tar.gz](https://data.gtdb.ecogenomic.org/releases/release214/214.1/genomic_files_reps/gtdb_genomes_reps_r214.tar.gz) (74.90G) 	|
| GTDB R07-RS207 genomic representatives (66k) 	| [gtdb-rs207.genomic-reps.dna.k51.sbt.zip]( https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs207/gtdb-rs207.genomic-reps.dna.k51.sbt.zip) (4.4Gb) 	| [gtdb_genomes_reps_r207.tar.gz](https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_reps/gtdb_genomes_reps_r207.tar.gz) (60.91G) 	|


### Example

The [align](align.md) command using just the Sourmash RS214 database can be run as

```
tracm align -i example_sample_R1.fastq example_sample_R1.fastq -o output_folder --database gtdb-rs214-reps.k51.sbt.zip
```

To avoid repeatedly downloading the same reference genomes, the path to the untar'd GTDB database can be provided as 

```
tracm align -i example_sample_R1.fastq example_sample_R1.fastq -o output_folder --database gtdb-rs214-reps.k51.sbt.zip --refseqs  
```


## Creating a database

To investigate multiple strains within a species or if only a subset of species are of interest it is possible to create a custom reference database using the `build-db` command.

### Example

```
tracm build-db -i refA.fasta refB.fasta refC.fast -o custom_db.zip -t 10
```

Alternatively, you can pass many fasta files located in the same directory as

```
tracm build-db -i /path/to/references/*.fasta -o custom_db.zip -t 10
```

### Options

```
Builds a database for tracm

options:
  -h, --help            show this help message and exit
  -i INPUT_FILES [INPUT_FILES ...], --input INPUT_FILES [INPUT_FILES ...]
                        path to genome fasta files (one per reference genome)
  -o DBNAME, --output DBNAME
                        name of the database file
  --ksize KSIZE         the kmer length used in sourmash (default=51)
  --scale SCALE         the scale used in sourmash (default=100)
  -t N_CPU, --threads N_CPU
                        number of threads to use (default=1)
  --quiet               turns off some console output
```