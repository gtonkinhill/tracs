# Cluster

The **cluster** module takes pairwise distance estimates in the format produced by the **distance** module. It uses single-linkage hierarchical clustering to identify potential transmission clusters. To cluster using a SNP threshold (after filtering) of 10 the module can be run as

### Example

```
tracm cluster -d transmission_distances.csv -o clusters.csv -D filter -c 10
```

### Output

The resulting clustering is formatted as a comma separated file with the sample name followed by its cluster. For example, the following output indicates that samples A and B belong to the same transmission cluster while C is in a separate cluster.

```
sample,cluster
A,0
B,0
C,1
```

### Options

```
Groups samples into putative transmission clusters using single linkage
clustering

options:
  -h, --help            show this help message and exit
  --quiet               turns off some console output

Input/output:
  -d DISTANCE_FILE, --distances DISTANCE_FILE
                        Pairwise distance estimates obtained from running the
                        'distance' function
  -o OUTPUT_FILE, --output OUTPUT_FILE
                        name of the output file to store the resulting cluster
                        assignments

Cluster options:
  -c THRESHOLD, --threshold THRESHOLD
                        Distance threshold. Samples will be grouped together
                        if the distance between them is below this threshold.
  -D {snp,filter,direct,expectedK}, --distance {snp,filter,direct,expectedK}
                        The type of transmission distance to use. Can be one
                        of 'snp', 'filter', 'direct', 'expectedK'
```

