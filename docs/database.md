# Creating a database

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